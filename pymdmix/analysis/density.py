"""
Density Grid Calculation
========================

Calculate probe atom density grids from MD trajectories.
This is the core analysis for MDMix binding site detection.

The algorithm:
1. Create 3D grid encompassing the reference structure
2. For each frame, count probe atoms in each grid cell
3. Normalize by number of frames to get density
4. Optionally convert to free energy

Modes of operation:
- Standard probe density (DensityAction)
- Protein/solute density (DensityProteinAction)
- All heavy atoms density (DensityAllHAAction)
- cpptraj-based density (CpptrajDensityAction)

Features:
- Center of mass (COM) tracking with includeCOM/onlyCOM
- Subregion calculation for focused analysis
- Parallel processing with multiprocessing

Examples
--------
>>> from pymdmix.analysis import DensityAction
>>> from pymdmix.core import open_trajectory, load_structure
>>>
>>> traj = open_trajectory("system.prmtop", "aligned.nc")
>>> ref = load_structure("reference.pdb")
>>>
>>> action = DensityAction()
>>> result = action.run(
...     trajectory=traj,
...     reference=ref,
...     probe_selections={"OH": "resname ETA and name O1"},
...     spacing=0.5,
...     include_com=True,  # Also track center of mass
...     n_workers=4,  # Use parallel processing
... )
"""

from __future__ import annotations

import logging
import os
import tempfile
from dataclasses import dataclass
from multiprocessing import Lock, Process, Queue
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

from pymdmix.analysis.base import Action, ActionResult, register_action
from pymdmix.core.grid import Grid
from pymdmix.core.trajectory import TrajectoryReader

log = logging.getLogger(__name__)


@dataclass
class ProbeConfig:
    """Configuration for a probe to track."""
    name: str
    selection: str | None = None  # MDAnalysis selection string
    atom_indices: NDArray[np.int64] | None = None  # Direct indices
    is_com: bool = False  # Is this a center of mass probe?
    residue_indices: list[NDArray[np.int64]] | None = None  # For COM: indices per residue

    def __post_init__(self):
        if self.selection is None and self.atom_indices is None and not self.is_com:
            raise ValueError("Either selection, atom_indices, or is_com must be provided")


@dataclass 
class Subregion:
    """Defines a subregion for focused density calculation."""
    min_coord: tuple[float, float, float]
    max_coord: tuple[float, float, float]

    def __post_init__(self):
        for i in range(3):
            if self.min_coord[i] >= self.max_coord[i]:
                raise ValueError(f"min_coord[{i}] must be < max_coord[{i}]")

    def filter_coords(self, coords: NDArray) -> NDArray:
        """Filter coordinates to only those within the subregion."""
        if coords.size == 0:
            return coords
        x, y, z = coords.T
        valid_x = (x >= self.min_coord[0]) & (x < self.max_coord[0])
        valid_y = (y >= self.min_coord[1]) & (y < self.max_coord[1])
        valid_z = (z >= self.min_coord[2]) & (z < self.max_coord[2])
        return coords[valid_x & valid_y & valid_z]


class DensityWorker(Process):
    """
    Multiprocessing worker for parallel density calculation.
    
    Reads frames from a queue and updates shared memory-mapped count grids.
    """

    def __init__(
        self,
        frame_queue: Queue,
        probe_configs: list[ProbeConfig],
        count_grids: dict[str, np.memmap],
        grid_origin: NDArray,
        grid_shape: tuple[int, int, int],
        grid_spacing: float,
        subregion: Subregion | None = None,
    ):
        super().__init__()
        self.frame_queue = frame_queue
        self.probe_configs = probe_configs
        self.count_grids = count_grids
        self.grid_origin = grid_origin
        self.grid_shape = grid_shape
        self.grid_spacing = grid_spacing
        self.subregion = subregion

    def _coord_to_index(self, coord: NDArray) -> tuple[int, int, int] | None:
        """Convert coordinate to grid index."""
        idx = np.floor((coord - self.grid_origin) / self.grid_spacing).astype(int)
        if np.all(idx >= 0) and np.all(idx < self.grid_shape):
            return tuple(idx)
        return None

    def run(self):
        """Process frames from queue until None is received."""
        while True:
            item = self.frame_queue.get()
            if item is None:
                break

            frame_coords, frame_num = item

            for probe in self.probe_configs:
                if probe.is_com and probe.residue_indices is not None:
                    # Calculate center of mass for each residue
                    coords_list = []
                    for res_indices in probe.residue_indices:
                        res_coords = frame_coords[res_indices]
                        com = res_coords.mean(axis=0)
                        coords_list.append(com)
                    coords = np.array(coords_list) if coords_list else np.empty((0, 3))
                elif probe.atom_indices is not None:
                    coords = frame_coords[probe.atom_indices]
                else:
                    continue

                # Apply subregion filter
                if self.subregion is not None:
                    coords = self.subregion.filter_coords(coords)

                if coords.size == 0:
                    continue

                # Add counts to grid
                for coord in coords:
                    idx = self._coord_to_index(coord)
                    if idx is not None:
                        self.count_grids[probe.name][idx] += 1


class ProteinDensityWorker(Process):
    """
    Worker for protein/solute density calculation with locking.
    
    Uses a lock because all workers write to the same grid.
    """

    def __init__(
        self,
        frame_queue: Queue,
        lock: Lock,
        atom_mask: NDArray[np.bool_],
        count_grid: np.memmap,
        grid_origin: NDArray,
        grid_shape: tuple[int, int, int],
        grid_spacing: float,
        subregion: Subregion | None = None,
    ):
        super().__init__()
        self.frame_queue = frame_queue
        self.lock = lock
        self.atom_mask = atom_mask
        self.count_grid = count_grid
        self.grid_origin = grid_origin
        self.grid_shape = grid_shape
        self.grid_spacing = grid_spacing
        self.subregion = subregion

    def _coord_to_index(self, coord: NDArray) -> tuple[int, int, int] | None:
        """Convert coordinate to grid index."""
        idx = np.floor((coord - self.grid_origin) / self.grid_spacing).astype(int)
        if np.all(idx >= 0) and np.all(idx < self.grid_shape):
            return tuple(idx)
        return None

    def run(self):
        """Process frames from queue until None is received."""
        while True:
            item = self.frame_queue.get()
            if item is None:
                break

            frame_coords, frame_num = item
            coords = frame_coords[self.atom_mask]

            if self.subregion is not None:
                coords = self.subregion.filter_coords(coords)

            if coords.size == 0:
                continue

            # Add counts with lock
            with self.lock:
                for coord in coords:
                    idx = self._coord_to_index(coord)
                    if idx is not None:
                        self.count_grid[idx] += 1


class AllHADensityWorker(Process):
    """
    Worker for all heavy atoms density calculation.
    
    Tracks each heavy atom type and COM for each residue type.
    """

    def __init__(
        self,
        frame_queue: Queue,
        ha_info: dict[str, dict[str, list[int]]],
        count_grids: dict[str, np.memmap],
        grid_origin: NDArray,
        grid_shape: tuple[int, int, int],
        grid_spacing: float,
        subregion: Subregion | None = None,
    ):
        super().__init__()
        self.frame_queue = frame_queue
        self.ha_info = ha_info  # {resname: {atomname: [indices]}}
        self.count_grids = count_grids
        self.grid_origin = grid_origin
        self.grid_shape = grid_shape
        self.grid_spacing = grid_spacing
        self.subregion = subregion

    def _coord_to_index(self, coord: NDArray) -> tuple[int, int, int] | None:
        """Convert coordinate to grid index."""
        idx = np.floor((coord - self.grid_origin) / self.grid_spacing).astype(int)
        if np.all(idx >= 0) and np.all(idx < self.grid_shape):
            return tuple(idx)
        return None

    def run(self):
        """Process frames from queue until None is received."""
        while True:
            item = self.frame_queue.get()
            if item is None:
                break

            frame_coords, frame_num = item

            for resname, atoms in self.ha_info.items():
                # Process each heavy atom type
                for atomname, indices in atoms.items():
                    probe_name = f"{resname}_{atomname}"
                    if probe_name not in self.count_grids:
                        continue

                    coords = frame_coords[indices]

                    if self.subregion is not None:
                        coords = self.subregion.filter_coords(coords)

                    for coord in coords:
                        idx = self._coord_to_index(coord)
                        if idx is not None:
                            self.count_grids[probe_name][idx] += 1

                # Calculate COM for residue
                com_probe = f"{resname}_COM"
                if com_probe in self.count_grids:
                    # Get all atom indices for this residue type
                    all_indices = []
                    for indices in atoms.values():
                        all_indices.extend(indices)
                    if all_indices:
                        all_coords = frame_coords[all_indices]
                        com = all_coords.mean(axis=0)
                        if self.subregion is None or len(self.subregion.filter_coords(com.reshape(1, 3))) > 0:
                            idx = self._coord_to_index(com)
                            if idx is not None:
                                self.count_grids[com_probe][idx] += 1


def _create_memmap_grid(shape: tuple[int, int, int]) -> np.memmap:
    """Create a temporary memory-mapped array for grid counts."""
    tmp = tempfile.mktemp(prefix='pymdmix_mmap_')
    return np.memmap(tmp, mode='w+', dtype='uint32', shape=shape)


def _calculate_grid_params(
    ref_coords: NDArray,
    spacing: float,
    padding: float,
    subregion: Subregion | None = None,
) -> tuple[NDArray, tuple[int, int, int]]:
    """
    Calculate grid origin and shape from reference coordinates.
    
    Returns
    -------
    origin : NDArray
        Grid origin (3,)
    shape : tuple[int, int, int]
        Grid dimensions
    """
    if subregion is not None:
        min_coord = np.array(subregion.min_coord)
        max_coord = np.array(subregion.max_coord)
    else:
        min_coord = ref_coords.min(axis=0) - padding
        max_coord = ref_coords.max(axis=0) + padding

    # Snap origin to nearest half-integer
    origin = np.trunc(min_coord)
    mask = np.abs(min_coord - origin) >= 0.5
    origin[~mask] -= 0.5
    origin[mask] += 0.5

    dimensions = max_coord - origin
    shape = tuple(np.ceil(dimensions / spacing).astype(int))

    return origin, shape


@register_action("density")
class DensityAction(Action):
    """
    Calculate probe density grids from trajectory.

    For each probe, creates a 3D grid counting atom occurrences,
    then normalizes to density (counts per frame).

    Parameters (in run())
    ---------------------
    trajectory : TrajectoryReader
        Aligned trajectory to analyze
    reference : parmed.Structure
        Reference structure for grid bounds
    probe_selections : dict[str, str]
        Probe name -> MDAnalysis selection string
    probe_indices : dict[str, NDArray] | None
        Alternative: probe name -> atom indices
    spacing : float
        Grid spacing in Angstroms (default: 0.5)
    padding : float
        Padding around reference structure (default: 5.0)
    output_dir : Path | None
        Output directory (default: current)
    output_prefix : str
        Prefix for output files (default: "")
    compute_free_energy : bool
        Also compute free energy grids (default: False)
    temperature : float
        Temperature for free energy (default: 300.0 K)
    include_com : bool
        Include center of mass probes (default: False)
    only_com : bool
        Only calculate COM probes (default: False)
    com_residue_indices : dict[str, list[NDArray]] | None
        For COM: probe_name -> list of atom indices per residue
    subregion : tuple[tuple[float,float,float], tuple[float,float,float]] | None
        ((min_x, min_y, min_z), (max_x, max_y, max_z)) for focused calculation
    n_workers : int
        Number of parallel workers (default: 1, sequential)

    Outputs
    -------
    - {prefix}{probe}_density.dx : Density grid for each probe
    - {prefix}{probe}_dg.dx : Free energy grid (if requested)
    """

    name = "density"
    description = "Calculate probe density grids from aligned trajectory"

    def run(
        self,
        trajectory: TrajectoryReader,
        reference=None,
        probe_selections: dict[str, str] | None = None,
        probe_indices: dict[str, NDArray] | None = None,
        spacing: float = 0.5,
        padding: float = 5.0,
        output_dir: Path | None = None,
        output_prefix: str = "",
        compute_free_energy: bool = False,
        temperature: float = 300.0,
        include_com: bool = False,
        only_com: bool = False,
        com_residue_indices: dict[str, list[NDArray]] | None = None,
        subregion: tuple[tuple[float, float, float], tuple[float, float, float]] | None = None,
        n_workers: int = 1,
        **kwargs,
    ) -> ActionResult:
        """Execute density calculation."""

        output_dir = Path(output_dir) if output_dir else Path.cwd()
        output_dir.mkdir(parents=True, exist_ok=True)

        # Parse subregion
        sub = None
        if subregion is not None:
            sub = Subregion(min_coord=subregion[0], max_coord=subregion[1])
            self.log.info(f"Using subregion: {subregion}")

        # Build probe configurations
        probes = self._build_probe_configs(
            trajectory, probe_selections, probe_indices,
            include_com=include_com, only_com=only_com,
            com_residue_indices=com_residue_indices,
        )

        if not probes:
            return ActionResult(
                success=False,
                error="No probes configured. Provide probe_selections or probe_indices.",
            )

        self.log.info(f"Calculating density for {len(probes)} probes")
        for p in probes:
            if p.is_com:
                n_residues = len(p.residue_indices) if p.residue_indices else 0
                self.log.info(f"  - {p.name}: COM of {n_residues} residues")
            else:
                n_atoms = len(p.atom_indices) if p.atom_indices is not None else "?"
                self.log.info(f"  - {p.name}: {n_atoms} atoms")

        # Get reference coordinates for grid bounds
        if reference is not None:
            ref_coords = reference.coordinates
        else:
            self.log.warning("No reference provided, using first frame for grid bounds")
            for frame in trajectory:
                ref_coords = frame.coordinates
                break

        # Calculate grid parameters
        grid_origin, grid_shape = _calculate_grid_params(
            ref_coords, spacing, padding, sub
        )
        self.log.debug(f"Grid origin: {grid_origin}, shape: {grid_shape}")

        # Decide sequential vs parallel
        if n_workers > 1:
            self.log.info(f"Using parallel processing with {n_workers} workers")
            results = self._run_parallel(
                trajectory, probes, grid_origin, grid_shape, spacing, sub, n_workers
            )
        else:
            self.log.info("Using sequential processing")
            results = self._run_sequential(
                trajectory, probes, grid_origin, grid_shape, spacing, sub
            )

        count_grids, n_frames = results

        self.log.info(f"Processed {n_frames} frames")

        # Convert to density and save
        output_files = []

        for probe in probes:
            count_data = count_grids[probe.name]

            # Create Grid object and normalize
            grid = Grid(
                data=count_data.astype(np.float64) / n_frames,
                origin=grid_origin,
                spacing=spacing,
            )

            # Density output
            density_path = output_dir / f"{output_prefix}{probe.name}_density.dx"
            grid.write_dx(density_path)
            output_files.append(density_path)
            self.log.info(f"Wrote {density_path}")

            # Free energy (optional)
            if compute_free_energy:
                dg = grid.to_free_energy(temperature)
                dg_path = output_dir / f"{output_prefix}{probe.name}_dg.dx"
                dg.write_dx(dg_path)
                output_files.append(dg_path)
                self.log.info(f"Wrote {dg_path}")

            # Clean up memmap temp file
            if hasattr(count_data, 'filename') and count_data.filename:
                try:
                    os.remove(count_data.filename)
                except OSError:
                    pass

        return ActionResult(
            success=True,
            output_files=output_files,
            metadata={
                "n_frames": n_frames,
                "n_probes": len(probes),
                "probe_names": [p.name for p in probes],
                "spacing": spacing,
                "grid_shape": grid_shape,
                "parallel_workers": n_workers,
                "subregion": subregion,
            },
        )

    def _run_sequential(
        self,
        trajectory: TrajectoryReader,
        probes: list[ProbeConfig],
        grid_origin: NDArray,
        grid_shape: tuple[int, int, int],
        spacing: float,
        subregion: Subregion | None,
    ) -> tuple[dict[str, NDArray], int]:
        """Run density calculation sequentially."""
        # Create count grids
        count_grids = {probe.name: np.zeros(grid_shape, dtype=np.uint32) for probe in probes}

        def coord_to_index(coord: NDArray) -> tuple[int, int, int] | None:
            idx = np.floor((coord - grid_origin) / spacing).astype(int)
            if np.all(idx >= 0) and np.all(idx < grid_shape):
                return tuple(idx)
            return None

        n_frames = 0
        for frame in trajectory:
            n_frames += 1

            if n_frames % 100 == 0:
                self.log.debug(f"Processing frame {n_frames}")

            for probe in probes:
                if probe.is_com and probe.residue_indices is not None:
                    coords_list = []
                    for res_indices in probe.residue_indices:
                        res_coords = frame.coordinates[res_indices]
                        com = res_coords.mean(axis=0)
                        coords_list.append(com)
                    coords = np.array(coords_list) if coords_list else np.empty((0, 3))
                elif probe.atom_indices is not None:
                    coords = frame.coordinates[probe.atom_indices]
                else:
                    continue

                if subregion is not None:
                    coords = subregion.filter_coords(coords)

                if coords.size == 0:
                    continue

                for coord in coords:
                    idx = coord_to_index(coord)
                    if idx is not None:
                        count_grids[probe.name][idx] += 1

        return count_grids, n_frames

    def _run_parallel(
        self,
        trajectory: TrajectoryReader,
        probes: list[ProbeConfig],
        grid_origin: NDArray,
        grid_shape: tuple[int, int, int],
        spacing: float,
        subregion: Subregion | None,
        n_workers: int,
    ) -> tuple[dict[str, np.memmap], int]:
        """Run density calculation in parallel using multiprocessing."""
        # Create memory-mapped count grids
        count_grids = {probe.name: _create_memmap_grid(grid_shape) for probe in probes}

        # Create frame queue
        frame_queue = Queue(maxsize=n_workers * 2)

        # Create workers
        workers = []
        for _ in range(n_workers):
            worker = DensityWorker(
                frame_queue=frame_queue,
                probe_configs=probes,
                count_grids=count_grids,
                grid_origin=grid_origin,
                grid_shape=grid_shape,
                grid_spacing=spacing,
                subregion=subregion,
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

        return count_grids, n_frames

    def _build_probe_configs(
        self,
        trajectory: TrajectoryReader,
        probe_selections: dict[str, str] | None,
        probe_indices: dict[str, NDArray] | None,
        include_com: bool = False,
        only_com: bool = False,
        com_residue_indices: dict[str, list[NDArray]] | None = None,
    ) -> list[ProbeConfig]:
        """Build probe configurations from selections or indices."""
        probes = []

        if only_com:
            # Only COM probes
            if com_residue_indices:
                for name, res_indices in com_residue_indices.items():
                    probes.append(ProbeConfig(
                        name=name,
                        is_com=True,
                        residue_indices=res_indices,
                    ))
            return probes

        # From selections (requires MDAnalysis trajectory)
        if probe_selections:
            if hasattr(trajectory, 'select_atoms'):
                for name, selection in probe_selections.items():
                    try:
                        indices = trajectory.select_atoms(selection)
                        probes.append(ProbeConfig(
                            name=name,
                            selection=selection,
                            atom_indices=indices,
                        ))
                        self.log.debug(f"Probe {name}: {len(indices)} atoms from '{selection}'")
                    except Exception as e:
                        self.log.warning(f"Failed to select probe {name}: {e}")
            else:
                self.log.warning(
                    "probe_selections requires MDAnalysis trajectory reader. "
                    "Use probe_indices instead."
                )

        # From direct indices
        if probe_indices:
            for name, indices in probe_indices.items():
                probes.append(ProbeConfig(
                    name=name,
                    atom_indices=np.asarray(indices),
                ))

        # Add COM probes if requested
        if include_com and com_residue_indices:
            for name, res_indices in com_residue_indices.items():
                probes.append(ProbeConfig(
                    name=f"{name}_COM",
                    is_com=True,
                    residue_indices=res_indices,
                ))

        return probes

    def validate(self, trajectory, **kwargs) -> list[str]:
        """Validate inputs."""
        errors = super().validate(trajectory, **kwargs)

        probe_selections = kwargs.get('probe_selections')
        probe_indices = kwargs.get('probe_indices')
        only_com = kwargs.get('only_com', False)
        com_residue_indices = kwargs.get('com_residue_indices')

        if not probe_selections and not probe_indices and not only_com:
            errors.append("Either probe_selections, probe_indices, or only_com must be provided")

        if only_com and not com_residue_indices:
            errors.append("only_com requires com_residue_indices")

        subregion = kwargs.get('subregion')
        if subregion:
            if len(subregion) != 2 or len(subregion[0]) != 3 or len(subregion[1]) != 3:
                errors.append("subregion must be ((x0,y0,z0), (x1,y1,z1))")

        return errors


@register_action("density_protein")
class DensityProteinAction(Action):
    """
    Calculate density grid for protein/solute atoms.

    Useful for:
    - Analyzing protein flexibility
    - Checking alignment quality
    - Visualizing protein motion

    Parameters (in run())
    ---------------------
    trajectory : TrajectoryReader
        Aligned trajectory to analyze
    reference : parmed.Structure
        Reference structure (solute part used for grid bounds)
    solute_mask : NDArray[bool]
        Boolean mask for solute atoms
    spacing : float
        Grid spacing in Angstroms (default: 0.5)
    padding : float
        Padding around solute (default: 5.0)
    output_dir : Path | None
        Output directory
    output_prefix : str
        Prefix for output files
    subregion : tuple | None
        Subregion for focused calculation
    n_workers : int
        Number of parallel workers

    Outputs
    -------
    - {prefix}protein_density.dx : Protein density grid
    """

    name = "density_protein"
    description = "Calculate density grid for protein/solute atoms"

    def run(
        self,
        trajectory: TrajectoryReader,
        reference=None,
        solute_mask: NDArray[np.bool_] | None = None,
        solute_indices: NDArray[np.int64] | None = None,
        spacing: float = 0.5,
        padding: float = 5.0,
        output_dir: Path | None = None,
        output_prefix: str = "",
        subregion: tuple[tuple[float, float, float], tuple[float, float, float]] | None = None,
        n_workers: int = 1,
        **kwargs,
    ) -> ActionResult:
        """Execute protein density calculation."""

        output_dir = Path(output_dir) if output_dir else Path.cwd()
        output_dir.mkdir(parents=True, exist_ok=True)

        # Get solute mask
        if solute_mask is not None:
            mask = solute_mask
        elif solute_indices is not None:
            # Will be applied per frame
            mask = solute_indices
        else:
            return ActionResult(
                success=False,
                error="Either solute_mask or solute_indices must be provided",
            )

        # Parse subregion
        sub = None
        if subregion is not None:
            sub = Subregion(min_coord=subregion[0], max_coord=subregion[1])

        # Get reference coordinates (solute only)
        if reference is not None:
            if solute_mask is not None:
                ref_coords = reference.coordinates[solute_mask]
            else:
                ref_coords = reference.coordinates[solute_indices]
        else:
            self.log.warning("No reference provided, using first frame")
            for frame in trajectory:
                if isinstance(mask, np.ndarray) and mask.dtype == bool:
                    ref_coords = frame.coordinates[mask]
                else:
                    ref_coords = frame.coordinates[mask]
                break

        # Calculate grid parameters
        grid_origin, grid_shape = _calculate_grid_params(
            ref_coords, spacing, padding, sub
        )
        self.log.info(f"Grid origin: {grid_origin}, shape: {grid_shape}")

        if n_workers > 1:
            self.log.info(f"Using parallel processing with {n_workers} workers")
            count_grid, n_frames = self._run_parallel(
                trajectory, mask, grid_origin, grid_shape, spacing, sub, n_workers
            )
        else:
            self.log.info("Using sequential processing")
            count_grid, n_frames = self._run_sequential(
                trajectory, mask, grid_origin, grid_shape, spacing, sub
            )

        self.log.info(f"Processed {n_frames} frames")

        # Convert to density
        grid = Grid(
            data=count_grid.astype(np.float64) / n_frames,
            origin=grid_origin,
            spacing=spacing,
        )

        # Save
        output_path = output_dir / f"{output_prefix}protein_density.dx"
        grid.write_dx(output_path)
        self.log.info(f"Wrote {output_path}")

        # Clean up memmap if used
        if hasattr(count_grid, 'filename') and count_grid.filename:
            try:
                os.remove(count_grid.filename)
            except OSError:
                pass

        return ActionResult(
            success=True,
            output_files=[output_path],
            metadata={
                "n_frames": n_frames,
                "spacing": spacing,
                "grid_shape": grid_shape,
            },
        )

    def _run_sequential(
        self,
        trajectory: TrajectoryReader,
        mask: NDArray,
        grid_origin: NDArray,
        grid_shape: tuple[int, int, int],
        spacing: float,
        subregion: Subregion | None,
    ) -> tuple[NDArray, int]:
        """Run protein density calculation sequentially."""
        count_grid = np.zeros(grid_shape, dtype=np.uint32)

        def coord_to_index(coord: NDArray) -> tuple[int, int, int] | None:
            idx = np.floor((coord - grid_origin) / spacing).astype(int)
            if np.all(idx >= 0) and np.all(idx < grid_shape):
                return tuple(idx)
            return None

        n_frames = 0
        for frame in trajectory:
            n_frames += 1

            coords = frame.coordinates[mask]

            if subregion is not None:
                coords = subregion.filter_coords(coords)

            for coord in coords:
                idx = coord_to_index(coord)
                if idx is not None:
                    count_grid[idx] += 1

        return count_grid, n_frames

    def _run_parallel(
        self,
        trajectory: TrajectoryReader,
        mask: NDArray,
        grid_origin: NDArray,
        grid_shape: tuple[int, int, int],
        spacing: float,
        subregion: Subregion | None,
        n_workers: int,
    ) -> tuple[np.memmap, int]:
        """Run protein density calculation in parallel."""
        count_grid = _create_memmap_grid(grid_shape)
        frame_queue = Queue(maxsize=n_workers * 2)
        lock = Lock()

        workers = []
        for _ in range(n_workers):
            worker = ProteinDensityWorker(
                frame_queue=frame_queue,
                lock=lock,
                atom_mask=mask,
                count_grid=count_grid,
                grid_origin=grid_origin,
                grid_shape=grid_shape,
                grid_spacing=spacing,
                subregion=subregion,
            )
            workers.append(worker)
            worker.start()

        n_frames = 0
        for frame in trajectory:
            n_frames += 1
            frame_queue.put((frame.coordinates.copy(), n_frames))

        for _ in workers:
            frame_queue.put(None)

        for worker in workers:
            worker.join()

        return count_grid, n_frames


@register_action("density_all_ha")
class DensityAllHAAction(Action):
    """
    Calculate density grids for ALL heavy atoms in solvent.

    Unlike standard probe density which tracks specific probe atoms,
    this tracks every non-hydrogen atom in solvent residues, plus
    center of mass for each residue type.

    Useful for:
    - Detailed solvent distribution analysis
    - Comparing different atom types
    - Understanding solvent orientation

    Parameters (in run())
    ---------------------
    trajectory : TrajectoryReader
        Aligned trajectory
    reference : parmed.Structure
        Reference structure
    heavy_atom_info : dict[str, dict[str, list[int]]]
        {residue_name: {atom_name: [atom_indices]}}
    exclude_residues : list[str]
        Residue names to exclude (default: ["WAT", "HOH"])
    spacing : float
        Grid spacing
    padding : float
        Padding around reference
    output_dir : Path
        Output directory
    output_prefix : str
        Prefix for output files
    subregion : tuple | None
        Subregion for focused calculation
    n_workers : int
        Number of parallel workers

    Outputs
    -------
    - {prefix}{resname}_{atomname}_density.dx : Grid for each heavy atom type
    - {prefix}{resname}_COM_density.dx : Grid for each residue COM
    """

    name = "density_all_ha"
    description = "Calculate density grids for all heavy atoms in solvent"

    def run(
        self,
        trajectory: TrajectoryReader,
        reference=None,
        heavy_atom_info: dict[str, dict[str, list[int]]] | None = None,
        exclude_residues: list[str] | None = None,
        spacing: float = 0.5,
        padding: float = 5.0,
        output_dir: Path | None = None,
        output_prefix: str = "",
        subregion: tuple[tuple[float, float, float], tuple[float, float, float]] | None = None,
        n_workers: int = 1,
        **kwargs,
    ) -> ActionResult:
        """Execute all heavy atoms density calculation."""

        output_dir = Path(output_dir) if output_dir else Path.cwd()
        output_dir.mkdir(parents=True, exist_ok=True)

        if exclude_residues is None:
            exclude_residues = ["WAT", "HOH"]

        if heavy_atom_info is None:
            return ActionResult(
                success=False,
                error="heavy_atom_info is required",
            )

        # Filter excluded residues
        ha_info = {k: v for k, v in heavy_atom_info.items() if k not in exclude_residues}

        # Build list of all probes
        probe_names = []
        for resname, atoms in ha_info.items():
            for atomname in atoms.keys():
                probe_names.append(f"{resname}_{atomname}")
            probe_names.append(f"{resname}_COM")

        self.log.info(f"Calculating density for {len(probe_names)} atom types/COMs")

        # Parse subregion
        sub = None
        if subregion is not None:
            sub = Subregion(min_coord=subregion[0], max_coord=subregion[1])

        # Get reference coordinates
        if reference is not None:
            ref_coords = reference.coordinates
        else:
            for frame in trajectory:
                ref_coords = frame.coordinates
                break

        # Calculate grid parameters
        grid_origin, grid_shape = _calculate_grid_params(
            ref_coords, spacing, padding, sub
        )

        if n_workers > 1:
            self.log.info(f"Using parallel processing with {n_workers} workers")
            count_grids, n_frames = self._run_parallel(
                trajectory, ha_info, probe_names, grid_origin, grid_shape, spacing, sub, n_workers
            )
        else:
            self.log.info("Using sequential processing")
            count_grids, n_frames = self._run_sequential(
                trajectory, ha_info, probe_names, grid_origin, grid_shape, spacing, sub
            )

        self.log.info(f"Processed {n_frames} frames")

        # Save all grids
        output_files = []
        for probe_name in probe_names:
            count_data = count_grids[probe_name]

            grid = Grid(
                data=count_data.astype(np.float64) / n_frames,
                origin=grid_origin,
                spacing=spacing,
            )

            output_path = output_dir / f"{output_prefix}{probe_name}_density.dx"
            grid.write_dx(output_path)
            output_files.append(output_path)

            # Clean up memmap
            if hasattr(count_data, 'filename') and count_data.filename:
                try:
                    os.remove(count_data.filename)
                except OSError:
                    pass

        self.log.info(f"Wrote {len(output_files)} density grids")

        return ActionResult(
            success=True,
            output_files=output_files,
            metadata={
                "n_frames": n_frames,
                "n_probes": len(probe_names),
                "probe_names": probe_names,
            },
        )

    def _run_sequential(
        self,
        trajectory: TrajectoryReader,
        ha_info: dict[str, dict[str, list[int]]],
        probe_names: list[str],
        grid_origin: NDArray,
        grid_shape: tuple[int, int, int],
        spacing: float,
        subregion: Subregion | None,
    ) -> tuple[dict[str, NDArray], int]:
        """Run all HA density calculation sequentially."""
        count_grids = {name: np.zeros(grid_shape, dtype=np.uint32) for name in probe_names}

        def coord_to_index(coord: NDArray) -> tuple[int, int, int] | None:
            idx = np.floor((coord - grid_origin) / spacing).astype(int)
            if np.all(idx >= 0) and np.all(idx < grid_shape):
                return tuple(idx)
            return None

        n_frames = 0
        for frame in trajectory:
            n_frames += 1

            for resname, atoms in ha_info.items():
                # Process each heavy atom
                for atomname, indices in atoms.items():
                    probe_name = f"{resname}_{atomname}"
                    coords = frame.coordinates[indices]

                    if subregion is not None:
                        coords = subregion.filter_coords(coords)

                    for coord in coords:
                        idx = coord_to_index(coord)
                        if idx is not None:
                            count_grids[probe_name][idx] += 1

                # Calculate COM
                com_probe = f"{resname}_COM"
                all_indices = []
                for indices in atoms.values():
                    all_indices.extend(indices)
                if all_indices:
                    all_coords = frame.coordinates[all_indices]
                    com = all_coords.mean(axis=0)
                    if subregion is None or len(subregion.filter_coords(com.reshape(1, 3))) > 0:
                        idx = coord_to_index(com)
                        if idx is not None:
                            count_grids[com_probe][idx] += 1

        return count_grids, n_frames

    def _run_parallel(
        self,
        trajectory: TrajectoryReader,
        ha_info: dict[str, dict[str, list[int]]],
        probe_names: list[str],
        grid_origin: NDArray,
        grid_shape: tuple[int, int, int],
        spacing: float,
        subregion: Subregion | None,
        n_workers: int,
    ) -> tuple[dict[str, np.memmap], int]:
        """Run all HA density calculation in parallel."""
        count_grids = {name: _create_memmap_grid(grid_shape) for name in probe_names}
        frame_queue = Queue(maxsize=n_workers * 2)

        workers = []
        for _ in range(n_workers):
            worker = AllHADensityWorker(
                frame_queue=frame_queue,
                ha_info=ha_info,
                count_grids=count_grids,
                grid_origin=grid_origin,
                grid_shape=grid_shape,
                grid_spacing=spacing,
                subregion=subregion,
            )
            workers.append(worker)
            worker.start()

        n_frames = 0
        for frame in trajectory:
            n_frames += 1
            frame_queue.put((frame.coordinates.copy(), n_frames))

        for _ in workers:
            frame_queue.put(None)

        for worker in workers:
            worker.join()

        return count_grids, n_frames


@register_action("cpptraj_density")
class CpptrajDensityAction(Action):
    """
    Calculate density grids using cpptraj (Amber).

    This uses the external cpptraj tool for potentially faster
    density calculation, especially for large systems.

    Requires:
    - cpptraj installed and in PATH, or
    - AMBER_PTRAJ environment variable pointing to cpptraj

    Parameters (in run())
    ---------------------
    topology : Path
        Amber topology file (.prmtop)
    trajectory_pattern : str
        Glob pattern for aligned trajectory files
    probe_masks : dict[str, str]
        Probe name -> Amber mask (e.g., {"OH": ":ETA@O1"})
    grid_dimensions : tuple[int, int, int]
        Grid dimensions (nx, ny, nz)
    grid_origin : tuple[float, float, float]
        Grid origin coordinates
    grid_spacing : float
        Grid spacing in Angstroms
    include_com : bool
        Include center of mass calculations
    only_com : bool
        Only calculate COM probes
    output_dir : Path
        Output directory
    n_threads : int
        Number of cpptraj threads
    wait_completion : bool
        Wait for cpptraj to complete (default: True)

    Outputs
    -------
    - {probe}_density.dx : Density grid for each probe
    """

    name = "cpptraj_density"
    description = "Calculate density grids using cpptraj"

    def run(
        self,
        topology: Path,
        trajectory_pattern: str | list[str],
        probe_masks: dict[str, str],
        grid_dimensions: tuple[int, int, int],
        grid_origin: tuple[float, float, float],
        grid_spacing: float = 0.5,
        include_com: bool = False,
        only_com: bool = False,
        com_mask: str | None = None,
        output_dir: Path | None = None,
        output_prefix: str = "",
        n_threads: int | None = None,
        wait_completion: bool = True,
        frame_offset: int = 1,
        **kwargs,
    ) -> ActionResult:
        """Execute cpptraj-based density calculation."""
        import shutil
        import subprocess

        output_dir = Path(output_dir) if output_dir else Path.cwd()
        output_dir.mkdir(parents=True, exist_ok=True)

        # Find cpptraj
        cpptraj = os.environ.get('AMBER_PTRAJ', shutil.which('cpptraj'))
        if not cpptraj:
            return ActionResult(
                success=False,
                error="cpptraj not found. Set AMBER_PTRAJ or add cpptraj to PATH.",
            )

        self.log.info(f"Using cpptraj: {cpptraj}")

        # Build trajectory list
        if isinstance(trajectory_pattern, str):
            trajectories = sorted(Path().glob(trajectory_pattern))
        else:
            trajectories = [Path(p) for p in trajectory_pattern]

        if not trajectories:
            return ActionResult(
                success=False,
                error=f"No trajectories found matching {trajectory_pattern}",
            )

        # Calculate grid center from origin and dimensions
        grid_center = [
            grid_origin[i] + (grid_dimensions[i] * grid_spacing) / 2
            for i in range(3)
        ]

        # Determine which probes to calculate
        if only_com:
            probes_to_calc = {"COM": com_mask or probe_masks.get("COM", "")}
        else:
            probes_to_calc = dict(probe_masks)
            if include_com and com_mask:
                probes_to_calc["COM"] = com_mask

        output_files = []
        scripts = []

        for probe_name, mask in probes_to_calc.items():
            # Generate cpptraj script
            trajin_lines = "\n".join(f"trajin {traj}" for traj in trajectories)

            is_com = "COM" in probe_name.upper()
            grid_out = output_dir / f"{output_prefix}{probe_name}.dx"

            # Build grid command
            # Note: cpptraj grid command syntax
            if is_com:
                grid_cmd = (
                    f"grid {grid_out} "
                    f"{grid_dimensions[0]} {grid_spacing} "
                    f"{grid_dimensions[1]} {grid_spacing} "
                    f"{grid_dimensions[2]} {grid_spacing} "
                    f"gridcenter {grid_center[0]} {grid_center[1]} {grid_center[2]} "
                    f"{mask} center"
                )
            else:
                grid_cmd = (
                    f"grid {grid_out} "
                    f"{grid_dimensions[0]} {grid_spacing} "
                    f"{grid_dimensions[1]} {grid_spacing} "
                    f"{grid_dimensions[2]} {grid_spacing} "
                    f"gridcenter {grid_center[0]} {grid_center[1]} {grid_center[2]} "
                    f"{mask}"
                )

            script_content = f"""# cpptraj density script for {probe_name}
parm {topology}
{trajin_lines}
{grid_cmd}
run
quit
"""

            script_path = output_dir / f"{output_prefix}{probe_name}.cpptraj"
            with open(script_path, 'w') as f:
                f.write(script_content)
            scripts.append(script_path)

            self.log.info(f"Wrote cpptraj script: {script_path}")

        # Run cpptraj for each probe
        for script_path in scripts:
            probe_name = script_path.stem.replace(output_prefix, "")
            log_path = script_path.with_suffix('.log')

            cmd = [cpptraj, '-i', str(script_path)]
            if n_threads:
                # cpptraj OpenMP threads via environment
                env = os.environ.copy()
                env['OMP_NUM_THREADS'] = str(n_threads)
            else:
                env = None

            self.log.info(f"Running cpptraj for {probe_name}")

            if wait_completion:
                with open(log_path, 'w') as log_file:
                    result = subprocess.run(
                        cmd,
                        stdout=log_file,
                        stderr=subprocess.STDOUT,
                        env=env,
                        cwd=str(output_dir),
                    )
                if result.returncode != 0:
                    self.log.warning(f"cpptraj returned {result.returncode} for {probe_name}")
            else:
                # Background execution
                subprocess.Popen(
                    cmd,
                    stdout=open(log_path, 'w'),
                    stderr=subprocess.STDOUT,
                    env=env,
                    cwd=str(output_dir),
                )

            grid_out = output_dir / f"{output_prefix}{probe_name}.dx"
            if grid_out.exists():
                output_files.append(grid_out)

        return ActionResult(
            success=True,
            output_files=output_files,
            metadata={
                "cpptraj": cpptraj,
                "n_trajectories": len(trajectories),
                "n_probes": len(probes_to_calc),
                "grid_dimensions": grid_dimensions,
                "scripts": [str(s) for s in scripts],
            },
        )


def calculate_density(
    trajectory: TrajectoryReader,
    probe_indices: dict[str, NDArray],
    reference_coords: NDArray,
    spacing: float = 0.5,
    padding: float = 5.0,
    n_workers: int = 1,
    subregion: tuple[tuple[float, float, float], tuple[float, float, float]] | None = None,
) -> dict[str, Grid]:
    """
    Convenience function to calculate density grids.

    Parameters
    ----------
    trajectory : TrajectoryReader
        Aligned trajectory
    probe_indices : dict[str, NDArray]
        Probe name -> atom indices
    reference_coords : NDArray
        Reference coordinates for grid bounds (n_atoms, 3)
    spacing : float
        Grid spacing in Angstroms
    padding : float
        Padding around reference
    n_workers : int
        Number of parallel workers
    subregion : tuple | None
        Subregion for focused calculation

    Returns
    -------
    dict[str, Grid]
        Probe name -> density grid
    """
    action = DensityAction()
    result = action.run(
        trajectory=trajectory,
        reference_coords=reference_coords,
        probe_indices=probe_indices,
        spacing=spacing,
        padding=padding,
        n_workers=n_workers,
        subregion=subregion,
    )

    if not result.success:
        raise RuntimeError(result.error)

    # Read back the grids
    grids = {}
    for path in result.output_files:
        if path.suffix == '.dx':
            name = path.stem.replace('_density', '')
            grids[name] = Grid.read_dx(path)

    return grids
