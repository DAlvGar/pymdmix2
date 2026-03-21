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
... )
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import logging
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

    def __post_init__(self):
        if self.selection is None and self.atom_indices is None:
            raise ValueError("Either selection or atom_indices must be provided")


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
        **kwargs,
    ) -> ActionResult:
        """Execute density calculation."""

        output_dir = Path(output_dir) if output_dir else Path.cwd()
        output_dir.mkdir(parents=True, exist_ok=True)

        # Build probe configurations
        probes = self._build_probe_configs(
            trajectory, probe_selections, probe_indices
        )

        if not probes:
            return ActionResult(
                success=False,
                error="No probes configured. Provide probe_selections or probe_indices.",
            )

        self.log.info(f"Calculating density for {len(probes)} probes")
        for p in probes:
            n_atoms = len(p.atom_indices) if p.atom_indices is not None else "?"
            self.log.info(f"  - {p.name}: {n_atoms} atoms")

        # Get reference coordinates for grid bounds
        if reference is not None:
            ref_coords = reference.coordinates
        else:
            # Use first frame as reference
            self.log.warning("No reference provided, using first frame for grid bounds")
            for frame in trajectory:
                ref_coords = frame.coordinates
                break

        # Create grids for each probe
        grids: dict[str, Grid] = {}
        for probe in probes:
            grids[probe.name] = Grid.from_structure(
                ref_coords, spacing=spacing, padding=padding
            )

        # Process trajectory
        n_frames = 0
        for frame in trajectory:
            n_frames += 1

            if n_frames % 100 == 0:
                self.log.debug(f"Processing frame {n_frames}")

            for probe in probes:
                if probe.atom_indices is not None:
                    coords = frame.coordinates[probe.atom_indices]
                    grids[probe.name].add_counts_bulk(coords)

        self.log.info(f"Processed {n_frames} frames")

        # Convert to density and save
        output_files = []

        for probe in probes:
            grid = grids[probe.name]

            # Density
            density = grid.to_density(n_frames)
            density_path = output_dir / f"{output_prefix}{probe.name}_density.dx"
            density.write_dx(density_path)
            output_files.append(density_path)
            self.log.info(f"Wrote {density_path}")

            # Free energy (optional)
            if compute_free_energy:
                dg = density.to_free_energy(temperature)
                dg_path = output_dir / f"{output_prefix}{probe.name}_dg.dx"
                dg.write_dx(dg_path)
                output_files.append(dg_path)
                self.log.info(f"Wrote {dg_path}")

        return ActionResult(
            success=True,
            output_files=output_files,
            metadata={
                "n_frames": n_frames,
                "n_probes": len(probes),
                "probe_names": [p.name for p in probes],
                "spacing": spacing,
                "grid_shape": grids[probes[0].name].shape,
            },
        )

    def _build_probe_configs(
        self,
        trajectory: TrajectoryReader,
        probe_selections: dict[str, str] | None,
        probe_indices: dict[str, NDArray] | None,
    ) -> list[ProbeConfig]:
        """Build probe configurations from selections or indices."""
        probes = []

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

        return probes

    def validate(self, trajectory, **kwargs) -> list[str]:
        """Validate inputs."""
        errors = super().validate(trajectory, **kwargs)

        probe_selections = kwargs.get('probe_selections')
        probe_indices = kwargs.get('probe_indices')

        if not probe_selections and not probe_indices:
            errors.append("Either probe_selections or probe_indices must be provided")

        return errors


def calculate_density(
    trajectory: TrajectoryReader,
    probe_indices: dict[str, NDArray],
    reference_coords: NDArray,
    spacing: float = 0.5,
    padding: float = 5.0,
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

    Returns
    -------
    dict[str, Grid]
        Probe name -> density grid
    """
    # Create grids
    grids = {
        name: Grid.from_structure(reference_coords, spacing, padding)
        for name in probe_indices
    }

    # Count
    n_frames = 0
    for frame in trajectory:
        n_frames += 1
        for name, indices in probe_indices.items():
            coords = frame.coordinates[indices]
            grids[name].add_counts_bulk(coords)

    # Normalize
    return {
        name: grid.to_density(n_frames)
        for name, grid in grids.items()
    }
