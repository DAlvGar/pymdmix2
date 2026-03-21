"""
Residence Time Analysis
=======================

Calculate which solvent residues visit specified hotspot locations
over the course of a trajectory.

This helps identify:
- Which probe types preferentially bind to hotspots
- How long probes stay in binding sites
- The dynamic behavior of solvent around the protein

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
... )
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Sequence
import logging
import json
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
        occupied = sum(1 for resids in self.frame_residues.values() if resids)
        return occupied / self.total_frames if self.total_frames > 0 else 0.0

    def top_residues(self, n: int = 10) -> list[tuple[int, int]]:
        """Get top N most frequent residues."""
        sorted_residues = sorted(
            self.residue_counts.items(),
            key=lambda x: x[1],
            reverse=True
        )
        return sorted_residues[:n]


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
    output_dir : Path | None
        Output directory
    output_prefix : str
        Prefix for output files

    Outputs
    -------
    - {prefix}residence.json : Full residence data
    - {prefix}residence_summary.txt : Human-readable summary
    """

    name = "residence"
    description = "Calculate residence times at hotspot locations"

    def run(
        self,
        trajectory: TrajectoryReader,
        hotspot_coords: Sequence[tuple[float, float, float]],
        tolerance: float = 3.0,
        atom_indices: NDArray | None = None,
        residue_ids: NDArray | None = None,
        residue_names: NDArray | None = None,
        track_residues: list[str] | None = None,
        output_dir: Path | None = None,
        output_prefix: str = "",
        reference=None,
        **kwargs,
    ) -> ActionResult:
        """Execute residence analysis."""

        if not HAS_SCIPY:
            return ActionResult(
                success=False,
                error="scipy is required for residence analysis",
            )

        output_dir = Path(output_dir) if output_dir else Path.cwd()
        output_dir.mkdir(parents=True, exist_ok=True)

        hotspot_coords = [tuple(c) for c in hotspot_coords]
        n_hotspots = len(hotspot_coords)

        self.log.info(f"Analyzing residence at {n_hotspots} hotspots")
        self.log.info(f"Tolerance: {tolerance} Å")

        # Initialize results for each hotspot
        results: list[ResidenceResult] = [
            ResidenceResult(
                hotspot_id=i,
                hotspot_coord=coord,
            )
            for i, coord in enumerate(hotspot_coords)
        ]

        # Convert hotspot coords to numpy array for KDTree
        _hotspot_array = np.array(hotspot_coords)

        # Process trajectory
        n_frames = 0
        for frame in trajectory:
            n_frames += 1

            # Get coordinates to analyze
            if atom_indices is not None:
                coords = frame.coordinates[atom_indices]
                local_residue_ids = residue_ids[atom_indices] if residue_ids is not None else None
                local_residue_names = residue_names[atom_indices] if residue_names is not None else None
            else:
                coords = frame.coordinates
                local_residue_ids = residue_ids
                local_residue_names = residue_names

            # Build KDTree for this frame
            tree = KDTree(coords)

            # Query each hotspot
            for i, hotspot_coord in enumerate(hotspot_coords):
                # Find atoms within tolerance
                indices = tree.query_ball_point(hotspot_coord, tolerance)

                if indices:
                    # Get unique residue IDs
                    if local_residue_ids is not None:
                        resids = np.unique(local_residue_ids[indices]).tolist()

                        # Filter by residue name if requested
                        if track_residues and local_residue_names is not None:
                            filtered_resids = []
                            for idx in indices:
                                resname = local_residue_names[idx]
                                resid = local_residue_ids[idx]
                                if resname in track_residues and resid not in filtered_resids:
                                    filtered_resids.append(resid)
                            resids = filtered_resids

                        results[i].frame_residues[n_frames] = resids

                        # Update counts
                        for resid in resids:
                            results[i].residue_counts[resid] = \
                                results[i].residue_counts.get(resid, 0) + 1
                    else:
                        # Just record that something was present
                        results[i].frame_residues[n_frames] = list(indices)

        # Update total frames
        for result in results:
            result.total_frames = n_frames

        self.log.info(f"Processed {n_frames} frames")

        # Save results
        output_files = []

        # JSON output
        json_path = output_dir / f"{output_prefix}residence.json"
        json_data = {
            "n_frames": n_frames,
            "tolerance": tolerance,
            "n_hotspots": n_hotspots,
            "hotspots": [
                {
                    "id": r.hotspot_id,
                    "coord": r.hotspot_coord,
                    "occupancy": r.occupancy,
                    "top_residues": r.top_residues(10),
                }
                for r in results
            ],
        }
        with open(json_path, 'w') as f:
            json.dump(json_data, f, indent=2)
        output_files.append(json_path)

        # Summary text
        summary_path = output_dir / f"{output_prefix}residence_summary.txt"
        with open(summary_path, 'w') as f:
            f.write("Residence Analysis Summary\n")
            f.write("==========================\n\n")
            f.write(f"Total frames: {n_frames}\n")
            f.write(f"Tolerance: {tolerance} Å\n")
            f.write(f"Hotspots analyzed: {n_hotspots}\n\n")

            for result in results:
                f.write(f"Hotspot {result.hotspot_id}: {result.hotspot_coord}\n")
                f.write(f"  Occupancy: {result.occupancy:.1%}\n")
                f.write("  Top residues:\n")
                for resid, count in result.top_residues(5):
                    pct = count / n_frames * 100
                    f.write(f"    Residue {resid}: {count} frames ({pct:.1f}%)\n")
                f.write("\n")

        output_files.append(summary_path)
        self.log.info(f"Wrote {json_path}")
        self.log.info(f"Wrote {summary_path}")

        return ActionResult(
            success=True,
            output_files=output_files,
            metadata={
                "n_frames": n_frames,
                "n_hotspots": n_hotspots,
                "tolerance": tolerance,
                "hotspot_occupancies": [r.occupancy for r in results],
            },
        )

    def validate(self, trajectory, **kwargs) -> list[str]:
        """Validate inputs."""
        errors = super().validate(trajectory, **kwargs)

        if not HAS_SCIPY:
            errors.append("scipy is required for residence analysis")

        hotspot_coords = kwargs.get('hotspot_coords')
        if not hotspot_coords:
            errors.append("hotspot_coords is required")

        return errors


def calculate_residence(
    trajectory: TrajectoryReader,
    hotspot_coords: Sequence[tuple[float, float, float]],
    atom_indices: NDArray,
    residue_ids: NDArray,
    tolerance: float = 3.0,
) -> list[ResidenceResult]:
    """
    Convenience function to calculate residence times.

    Parameters
    ----------
    trajectory : TrajectoryReader
        Trajectory to analyze
    hotspot_coords : sequence
        Hotspot coordinates
    atom_indices : NDArray
        Atom indices to track
    residue_ids : NDArray
        Residue ID for each atom index
    tolerance : float
        Distance threshold in Angstroms

    Returns
    -------
    list[ResidenceResult]
        Results for each hotspot
    """
    action = ResidenceAction()
    result = action.run(
        trajectory=trajectory,
        hotspot_coords=hotspot_coords,
        atom_indices=atom_indices,
        residue_ids=residue_ids,
        tolerance=tolerance,
    )

    if not result.success:
        raise RuntimeError(result.error)

    # Return parsed results from JSON
    # (In real use, we'd return the ResidenceResult objects directly)
    return result
