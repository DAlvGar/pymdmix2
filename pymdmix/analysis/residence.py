"""
Residence Time Analysis
=======================

Calculate which solvent residues visit specified hotspot locations
over the course of a trajectory.

This helps identify:
- Which probe types preferentially bind to hotspots
- How long probes stay in binding sites
- The dynamic behavior of solvent around the protein

Features:
- Sphere-based hotspot tracking
- Residue name filtering
- Parallel processing with multiprocessing
- Detailed occupancy statistics

Examples
--------
>>> from pymdmix.analysis import ResidenceAction
>>>
>>> action = ResidenceAction()
>>> result = action.run(
...     trajectory=traj,
...     hotspot_coords=[(10.0, 20.0, 15.0), (25.0, 30.0, 12.0)],
...     tolerance=3.0,
...     track_residues=["ETA", "WAT"],
...     n_workers=4,  # Use parallel processing
... )
"""

from __future__ import annotations

import json
import logging
from collections.abc import Sequence
from dataclasses import dataclass, field
from multiprocessing import Lock, Manager, Process, Queue
from pathlib import Path
from typing import Any, cast

import numpy as np
from numpy.typing import NDArray

try:
    from scipy.spatial import KDTree

    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False
    KDTree = None

from pymdmix.analysis.base import Action, ActionResult, register_action
from pymdmix.core.trajectory import TrajectoryReader

log = logging.getLogger(__name__)


def _flatten_nested_list(nested_list: list) -> list:
    """Flatten a nested list to a single list."""
    result = []
    for item in nested_list:
        if isinstance(item, list):
            result.extend(_flatten_nested_list(item))
        else:
            result.append(item)
    return result


@dataclass
class ResidenceResult:
    """
    Result of residence analysis for a single hotspot.

    Attributes
    ----------
    hotspot_id : int
        Index of the hotspot
    hotspot_coord : tuple[float, float, float]
        Hotspot coordinates
    frame_residues : dict[int, list[int]]
        Frame number -> list of residue IDs present
    residue_counts : dict[int, int]
        Residue ID -> number of frames present
    total_frames : int
        Total frames analyzed
    """

    hotspot_id: int
    hotspot_coord: tuple[float, float, float]
    frame_residues: dict[int, list[int]] = field(default_factory=dict)
    residue_counts: dict[int, int] = field(default_factory=dict)
    total_frames: int = 0

    @property
    def occupancy(self) -> float:
        """Fraction of frames with any residue present."""
        occupied = sum(1 for resids in self.frame_residues.values() if resids and resids != [0])
        return occupied / self.total_frames if self.total_frames > 0 else 0.0

    def top_residues(self, n: int = 10) -> list[tuple[int, int]]:
        """Get top N most frequent residues."""
        sorted_residues = sorted(self.residue_counts.items(), key=lambda x: x[1], reverse=True)
        return sorted_residues[:n]


class ResidenceWorker(Process):
    """
    Multiprocessing worker for parallel residence calculation.

    Uses a shared Manager.dict for synchronized result storage.
    """

    def __init__(
        self,
        frame_queue: Queue,
        hotspot_coords: list[tuple[float, float, float]],
        tolerance: float,
        non_h_mask: NDArray[np.bool_],
        index_to_resid: dict[int, int],
        resid_to_resname: dict[int, str],
        track_residues: list[str] | None,
        results: Any,  # Manager.dict
        lock: Any,
    ):
        super().__init__()
        self.frame_queue = frame_queue
        self.hotspot_coords = hotspot_coords
        self.tolerance = tolerance
        self.non_h_mask = non_h_mask
        self.index_to_resid = index_to_resid
        self.resid_to_resname = resid_to_resname
        self.track_residues = track_residues
        self.results = results
        self.lock = lock

    def run(self):
        """Process frames from queue until None is received."""
        while True:
            item = self.frame_queue.get()
            if item is None:
                break

            frame_coords, frame_num = item

            # Apply non-hydrogen mask
            masked_coords = frame_coords[self.non_h_mask]

            # Build KDTree for this frame
            tree = KDTree(masked_coords)

            # Query all hotspots
            for hotspot_idx, hotspot_coord in enumerate(self.hotspot_coords):
                # Find atoms within tolerance
                indices_nested = tree.query_ball_point(hotspot_coord, self.tolerance)
                indices = (
                    _flatten_nested_list([indices_nested])
                    if isinstance(indices_nested, list)
                    else [indices_nested]
                )

                result_key = (hotspot_idx, frame_num)

                if indices:
                    # Get unique residue IDs
                    resids = np.unique(
                        [
                            self.index_to_resid.get(idx, 0)
                            for idx in indices
                            if idx in self.index_to_resid
                        ]
                    )
                    resnames = [self.resid_to_resname.get(resid, "") for resid in resids]

                    # Filter by tracked residue names if specified
                    if self.track_residues:
                        filtered_resids = []
                        for resid, resname in zip(resids, resnames):
                            if resname in self.track_residues:
                                filtered_resids.append(resid)
                        if filtered_resids:
                            resids = np.array(filtered_resids)
                        else:
                            resids = np.array([0])  # No tracked residues found

                    self.results[result_key] = resids.tolist()
                else:
                    self.results[result_key] = [0]


@register_action("residence")
class ResidenceAction(Action):
    """
    Calculate residence times of solvent at specified hotspot locations.

    Uses KDTree for efficient spatial queries to find which atoms
    are within a tolerance distance of each hotspot coordinate.

    Parameters (in run())
    ---------------------
    trajectory : TrajectoryReader
        Trajectory to analyze
    hotspot_coords : list[tuple[float, float, float]]
        List of hotspot (x, y, z) coordinates
    tolerance : float
        Distance threshold in Angstroms (default: 3.0)
    atom_indices : NDArray | None
        Indices of atoms to consider (e.g., probe atoms only)
    residue_ids : NDArray | None
        Residue ID for each atom (same length as coordinates)
    residue_names : NDArray | None
        Residue name for each atom
    track_residues : list[str] | None
        Only track these residue types
    non_hydrogen_mask : NDArray[bool] | None
        Mask to exclude hydrogen atoms
    output_dir : Path | None
        Output directory
    output_prefix : str
        Prefix for output files
    n_workers : int
        Number of parallel workers (default: 1, sequential)

    Outputs
    -------
    - {prefix}residence.json : Full residence data
    - {prefix}residence_summary.txt : Human-readable summary
    """

    name = "residence"
    description = "Calculate residence times at hotspot locations"

    def run(
        self,
        trajectory: TrajectoryReader | None = None,
        reference=None,
        output_dir: Path | None = None,
        hotspot_coords: Sequence[tuple[float, float, float]] | None = None,
        tolerance: float = 3.0,
        atom_indices: NDArray | None = None,
        residue_ids: NDArray | None = None,
        residue_names: NDArray | None = None,
        track_residues: list[str] | None = None,
        non_hydrogen_mask: NDArray[np.bool_] | None = None,
        output_prefix: str = "",
        n_workers: int = 1,
        **kwargs,
    ) -> ActionResult:
        """Execute residence analysis."""

        if trajectory is None:
            return ActionResult(
                success=False,
                error="trajectory is required for ResidenceAction.run()",
            )

        if hotspot_coords is None:
            return ActionResult(
                success=False,
                error="hotspot_coords must be provided to ResidenceAction.run()",
            )

        if not HAS_SCIPY:
            return ActionResult(
                success=False,
                error="scipy is required for residence analysis",
            )

        output_dir = Path(output_dir) if output_dir else Path.cwd()
        output_dir.mkdir(parents=True, exist_ok=True)

        hotspot_coords = [
            cast(tuple[float, float, float], tuple(float(v) for v in c)) for c in hotspot_coords
        ]
        n_hotspots = len(hotspot_coords)

        self.log.info(f"Analyzing residence at {n_hotspots} hotspots")
        self.log.info(f"Tolerance: {tolerance} Å")

        if track_residues:
            self.log.info(f"Tracking residue types: {track_residues}")
        else:
            self.log.warning(
                "No track_residues specified. Will track ALL residues including protein!"
            )

        # Build residue mappings
        if non_hydrogen_mask is None:
            # Default: include all atoms
            non_hydrogen_mask = np.ones(
                residue_ids.shape[0] if residue_ids is not None else 0, dtype=bool
            )

        # Build index to residue ID mapping (for masked indices)
        if residue_ids is not None:
            masked_residue_ids = residue_ids[non_hydrogen_mask]
            index_to_resid = dict(enumerate(masked_residue_ids.astype(int)))
        else:
            index_to_resid = {}

        # Build residue ID to name mapping
        if residue_ids is not None and residue_names is not None:
            masked_resids = residue_ids[non_hydrogen_mask].astype(int)
            masked_resnames = residue_names[non_hydrogen_mask]
            resid_to_resname = dict(zip(masked_resids, masked_resnames))
            resid_to_resname[0] = "NO_RESIDENCE"  # Special marker for unoccupied
        else:
            resid_to_resname = {0: "NO_RESIDENCE"}

        # Run analysis
        if n_workers > 1:
            self.log.info(f"Using parallel processing with {n_workers} workers")
            frame_results, n_frames = self._run_parallel(
                trajectory,
                hotspot_coords,
                tolerance,
                non_hydrogen_mask,
                index_to_resid,
                resid_to_resname,
                track_residues,
                n_workers,
            )
        else:
            self.log.info("Using sequential processing")
            frame_results, n_frames = self._run_sequential(
                trajectory,
                hotspot_coords,
                tolerance,
                non_hydrogen_mask,
                index_to_resid,
                resid_to_resname,
                track_residues,
            )

        self.log.info(f"Processed {n_frames} frames")

        # Organize results by hotspot
        results: list[ResidenceResult] = [
            ResidenceResult(
                hotspot_id=i,
                hotspot_coord=coord,
                total_frames=n_frames,
            )
            for i, coord in enumerate(hotspot_coords)
        ]

        for key, resids in frame_results.items():
            if isinstance(key, tuple):
                hotspot_idx, frame_num = key
            else:
                # Handle dict returned from parallel processing
                continue

            results[hotspot_idx].frame_residues[frame_num] = resids

            # Update counts
            for resid in resids:
                if resid != 0:  # Don't count "unoccupied" marker
                    results[hotspot_idx].residue_counts[resid] = (
                        results[hotspot_idx].residue_counts.get(resid, 0) + 1
                    )

        # Save results
        output_files = []

        # JSON output
        json_path = output_dir / f"{output_prefix}residence.json"

        def _to_serializable(obj):
            """Convert numpy types to JSON-serializable Python types."""
            if isinstance(obj, np.integer):
                return int(obj)
            elif isinstance(obj, np.floating):
                return float(obj)
            elif isinstance(obj, np.ndarray):
                return obj.tolist()
            elif isinstance(obj, dict):
                # Convert both keys and values
                return {
                    int(k) if isinstance(k, np.integer) else k: _to_serializable(v)
                    for k, v in obj.items()
                }
            elif isinstance(obj, (list, tuple)):
                return [_to_serializable(v) for v in obj]
            return obj

        json_data = _to_serializable(
            {
                "n_frames": n_frames,
                "tolerance": tolerance,
                "n_hotspots": n_hotspots,
                "track_residues": track_residues,
                "hotspots": [
                    {
                        "id": r.hotspot_id,
                        "coord": r.hotspot_coord,
                        "occupancy": r.occupancy,
                        "top_residues": r.top_residues(10),
                        "residue_counts": r.residue_counts,
                    }
                    for r in results
                ],
            }
        )
        with open(json_path, "w") as f:
            json.dump(json_data, f, indent=2)
        output_files.append(json_path)

        # Summary text (similar to old format)
        summary_path = output_dir / f"{output_prefix}residence_summary.txt"
        with open(summary_path, "w") as f:
            f.write("Residence Analysis Summary\n")
            f.write("==========================\n\n")
            f.write(f"Total frames: {n_frames}\n")
            f.write(f"Tolerance: {tolerance} Å\n")
            f.write(f"Hotspots analyzed: {n_hotspots}\n")
            if track_residues:
                f.write(f"Tracked residues: {track_residues}\n")
            f.write("\n")

            # Build resname to resid map for output
            f.write("# RESNAME-RESID MAP\n")
            all_resids: set[int] = set()
            for r in results:
                all_resids.update(r.residue_counts.keys())
            resname_to_resids: dict[str, list[int]] = {}
            for resid in all_resids:
                resname = resid_to_resname.get(resid, "UNK")
                if resname not in resname_to_resids:
                    resname_to_resids[resname] = []
                resname_to_resids[resname].append(resid)
            for resname, resids in sorted(resname_to_resids.items()):
                f.write(f"# {resname} = {','.join(map(str, sorted(resids)))}\n")
            f.write("\n")

            for result in results:
                f.write(f"Hotspot {result.hotspot_id}: {result.hotspot_coord}\n")
                f.write(f"  Occupancy: {result.occupancy:.1%}\n")
                f.write("  Top residues:\n")
                for resid, count in result.top_residues(5):
                    pct = count / n_frames * 100
                    resname = resid_to_resname.get(resid, "UNK")
                    f.write(f"    {resname} (ID {resid}): {count} frames ({pct:.1f}%)\n")
                f.write("\n")

        output_files.append(summary_path)

        # ASCII data file (old format compatibility)
        data_path = output_dir / f"{output_prefix}residence_data.txt"
        with open(data_path, "w") as f:
            f.write(f"# Residence results for {n_hotspots} hotspots\n")
            if track_residues:
                f.write(f"# Tracked residues: {track_residues}\n")
            f.write("# RESNAME-RESID MAP\n")
            for resname, resids in sorted(resname_to_resids.items()):
                f.write(f"# {resname} = {','.join(map(str, sorted(resids)))}\n")
            f.write("# DATA (frame_num, hotspot_0_resids, hotspot_1_resids, ...)\n")

            for frame_num in sorted(
                set(fn for _, fn in frame_results.keys() if isinstance(_, int))
            ):
                row = [str(frame_num)]
                for hotspot_idx in range(n_hotspots):
                    resids = frame_results.get((hotspot_idx, frame_num), [0])
                    row.append(",".join(map(str, resids)))
                f.write("\t".join(row) + "\n")

        output_files.append(data_path)

        self.log.info(f"Wrote {json_path}")
        self.log.info(f"Wrote {summary_path}")
        self.log.info(f"Wrote {data_path}")

        return ActionResult(
            success=True,
            output_files=output_files,
            metadata={
                "n_frames": n_frames,
                "n_hotspots": n_hotspots,
                "tolerance": tolerance,
                "hotspot_occupancies": [r.occupancy for r in results],
                "parallel_workers": n_workers,
            },
        )

    def _run_sequential(
        self,
        trajectory: TrajectoryReader,
        hotspot_coords: list[tuple[float, float, float]],
        tolerance: float,
        non_h_mask: NDArray[np.bool_],
        index_to_resid: dict[int, int],
        resid_to_resname: dict[int, str],
        track_residues: list[str] | None,
    ) -> tuple[dict[tuple[int, int], list[int]], int]:
        """Run residence analysis sequentially."""
        results = {}
        n_frames = 0

        for frame in trajectory:
            n_frames += 1

            # Apply non-hydrogen mask
            masked_coords = frame.coordinates[non_h_mask]

            # Build KDTree for this frame
            tree = KDTree(masked_coords)

            for hotspot_idx, hotspot_coord in enumerate(hotspot_coords):
                # Find atoms within tolerance
                indices_nested = tree.query_ball_point(hotspot_coord, tolerance)
                indices = (
                    _flatten_nested_list([indices_nested])
                    if isinstance(indices_nested, list)
                    else [indices_nested]
                )

                result_key = (hotspot_idx, n_frames)

                if indices:
                    # Get unique residue IDs
                    resids = list(
                        set(index_to_resid.get(idx, 0) for idx in indices if idx in index_to_resid)
                    )
                    resnames = [resid_to_resname.get(resid, "") for resid in resids]

                    # Filter by tracked residue names
                    if track_residues:
                        filtered_resids = [
                            resid
                            for resid, resname in zip(resids, resnames)
                            if resname in track_residues
                        ]
                        resids = filtered_resids if filtered_resids else [0]

                    results[result_key] = list(np.unique(resids))
                else:
                    results[result_key] = [0]

        return results, n_frames

    def _run_parallel(
        self,
        trajectory: TrajectoryReader,
        hotspot_coords: list[tuple[float, float, float]],
        tolerance: float,
        non_h_mask: NDArray[np.bool_],
        index_to_resid: dict[int, int],
        resid_to_resname: dict[int, str],
        track_residues: list[str] | None,
        n_workers: int,
    ) -> tuple[dict, int]:
        """Run residence analysis in parallel using multiprocessing."""
        manager = Manager()
        results = manager.dict()
        lock = Lock()
        frame_queue: Any = Queue(maxsize=n_workers * 2)

        # Create workers
        workers = []
        for _ in range(n_workers):
            worker = ResidenceWorker(
                frame_queue=frame_queue,
                hotspot_coords=hotspot_coords,
                tolerance=tolerance,
                non_h_mask=non_h_mask,
                index_to_resid=index_to_resid,
                resid_to_resname=resid_to_resname,
                track_residues=track_residues,
                results=results,
                lock=lock,
            )
            workers.append(worker)
            worker.start()

        # Feed frames to workers
        n_frames = 0
        for frame in trajectory:
            n_frames += 1
            frame_queue.put((frame.coordinates.copy(), n_frames))

        # Signal workers to stop
        for _ in workers:
            frame_queue.put(None)

        # Wait for workers to finish
        for worker in workers:
            worker.join()

        # Convert manager.dict to regular dict
        return dict(results), n_frames

    def validate(self, trajectory, **kwargs) -> list[str]:
        """Validate inputs."""
        errors = super().validate(trajectory, **kwargs)

        if not HAS_SCIPY:
            errors.append("scipy is required for residence analysis")

        hotspot_coords = kwargs.get("hotspot_coords")
        if not hotspot_coords:
            errors.append("hotspot_coords is required")

        return cast(list[str], errors)


def calculate_residence(
    trajectory: TrajectoryReader,
    hotspot_coords: Sequence[tuple[float, float, float]],
    atom_indices: NDArray | None = None,
    residue_ids: NDArray | None = None,
    residue_names: NDArray | None = None,
    tolerance: float = 3.0,
    track_residues: list[str] | None = None,
    non_hydrogen_mask: NDArray[np.bool_] | None = None,
    n_workers: int = 1,
) -> ActionResult:
    """
    Convenience function to calculate residence times.

    Parameters
    ----------
    trajectory : TrajectoryReader
        Trajectory to analyze
    hotspot_coords : sequence
        Hotspot coordinates
    atom_indices : NDArray | None
        Atom indices to track
    residue_ids : NDArray | None
        Residue ID for each atom
    residue_names : NDArray | None
        Residue name for each atom
    tolerance : float
        Distance threshold in Angstroms
    track_residues : list[str] | None
        Residue types to track
    non_hydrogen_mask : NDArray[bool] | None
        Mask for non-hydrogen atoms
    n_workers : int
        Number of parallel workers

    Returns
    -------
    ActionResult
        Residence analysis action result
    """
    action = ResidenceAction()
    result = action.run(
        trajectory=trajectory,
        hotspot_coords=hotspot_coords,
        atom_indices=atom_indices,
        residue_ids=residue_ids,
        residue_names=residue_names,
        tolerance=tolerance,
        track_residues=track_residues,
        non_hydrogen_mask=non_hydrogen_mask,
        n_workers=n_workers,
    )

    if not result.success:
        raise RuntimeError(result.error)

    return result
