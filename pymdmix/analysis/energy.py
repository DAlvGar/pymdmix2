"""
Energy Conversion and Averaging
===============================

Convert density grids to free energy and perform Boltzmann
averaging across multiple replicas.

The free energy is calculated using:
    ΔG = -RT ln(ρ/ρ₀)

Where ρ is the local density and ρ₀ is the bulk/reference density.

Examples
--------
>>> from pymdmix.analysis.energy import boltzmann_average, density_to_free_energy
>>>
>>> # Average multiple replica grids
>>> avg_grid = boltzmann_average(["rep1.dx", "rep2.dx", "rep3.dx"])
>>>
>>> # Convert density to free energy
>>> fe_grid = density_to_free_energy(density_grid, temperature=300)
"""

from __future__ import annotations

import logging
from collections.abc import Sequence
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

from pymdmix.analysis.base import Action, ActionResult, register_action
from pymdmix.core.grid import Grid

log = logging.getLogger(__name__)

# Physical constants
R_KCAL = 0.001987204  # kcal/(mol·K)


@dataclass
class EnergyResult:
    """Result of energy conversion."""

    grids: dict[str, Grid] = field(default_factory=dict)
    temperature: float = 300.0
    method: str = "density_to_free_energy"
    n_replicas: int = 1

    def __repr__(self) -> str:
        n = len(self.grids)
        probes = list(self.grids.keys())
        return (
            f"EnergyResult(n_probes={n}, probes={probes}, "
            f"T={self.temperature}K, replicas={self.n_replicas})"
        )


def density_to_free_energy(
    grid: Grid,
    temperature: float = 300.0,
    reference_density: float | None = None,
    min_density: float = 1e-10,
) -> Grid:
    """
    Convert density grid to free energy.

    Uses ΔG = -RT ln(ρ/ρ₀)

    Parameters
    ----------
    grid : Grid
        Density grid
    temperature : float
        Temperature in Kelvin
    reference_density : float, optional
        Reference/bulk density. If None, uses mean density.
    min_density : float
        Minimum density to avoid log(0)

    Returns
    -------
    Grid
        Free energy grid in kcal/mol
    """
    RT = R_KCAL * temperature

    # Get reference density
    if reference_density is None:
        reference_density = float(np.mean(grid.data[grid.data > min_density]))

    if reference_density <= 0:
        raise ValueError("Reference density must be positive")

    # Avoid log(0)
    density = np.maximum(grid.data, min_density)

    # Calculate free energy
    # Higher density = more negative ΔG (favorable)
    free_energy = -RT * np.log(density / reference_density)

    return Grid(
        data=free_energy,
        origin=grid.origin,
        spacing=grid.spacing,
    )


def boltzmann_average(
    grids: Sequence[Grid | str | Path],
    temperature: float = 300.0,
    weights: Sequence[float] | None = None,
) -> Grid:
    """
    Calculate Boltzmann-weighted average of multiple grids.

    For density grids, this averages the densities.
    For energy grids, this performs proper Boltzmann averaging:
        <E> = -RT ln(Σᵢ wᵢ exp(-Eᵢ/RT))

    Parameters
    ----------
    grids : sequence of Grid or paths
        Input grids to average
    temperature : float
        Temperature in Kelvin
    weights : sequence of float, optional
        Weights for each grid. If None, equal weights.

    Returns
    -------
    Grid
        Averaged grid
    """
    if len(grids) == 0:
        raise ValueError("No grids provided")

    if len(grids) == 1:
        log.warning("Only one grid provided, returning unaveraged")
        grid = grids[0]
        if isinstance(grid, (str, Path)):
            return Grid.read_dx(grid)
        return grid

    # Load grids if paths
    loaded_grids = []
    for g in grids:
        if isinstance(g, (str, Path)):
            loaded_grids.append(Grid.read_dx(g))
        else:
            loaded_grids.append(g)

    # Validate compatible shapes
    ref_shape = loaded_grids[0].shape
    ref_origin = loaded_grids[0].origin
    ref_spacing = loaded_grids[0].spacing

    for i, g in enumerate(loaded_grids[1:], 2):
        if g.shape != ref_shape:
            raise ValueError(f"Grid {i} shape {g.shape} != reference {ref_shape}")
        if abs(g.spacing - ref_spacing) > 1e-6:
            raise ValueError(f"Grid {i} spacing {g.spacing} != reference {ref_spacing}")

    # Set weights
    if weights is None:
        weights = [1.0 / len(loaded_grids)] * len(loaded_grids)
    else:
        weights = list(weights)
        if len(weights) != len(loaded_grids):
            raise ValueError("Number of weights must match number of grids")
        # Normalize weights
        total = sum(weights)
        weights = [w / total for w in weights]

    # Check if these are energy grids (contain negative values)
    is_energy = any(g.data.min() < 0 for g in loaded_grids)

    if is_energy:
        # Boltzmann averaging for energy grids
        avg_data = _boltzmann_average_energy(loaded_grids, weights, temperature)
    else:
        # Simple weighted average for density grids
        avg_data = _weighted_average_density(loaded_grids, weights)

    return Grid(
        data=avg_data,
        origin=ref_origin,
        spacing=ref_spacing,
    )


def _weighted_average_density(
    grids: list[Grid],
    weights: list[float],
) -> NDArray[np.float64]:
    """Simple weighted average of density grids."""
    result = np.zeros_like(grids[0].data)

    for grid, weight in zip(grids, weights):
        result += weight * grid.data

    return result


def _boltzmann_average_energy(
    grids: list[Grid],
    weights: list[float],
    temperature: float,
) -> NDArray[np.float64]:
    """
    Boltzmann average of energy grids.

    <E> = -RT ln(Σᵢ wᵢ exp(-Eᵢ/RT))

    This properly handles the non-linear nature of free energy averaging.
    """
    RT = R_KCAL * temperature

    # Use log-sum-exp trick for numerical stability
    # log(Σ exp(xᵢ)) = max(x) + log(Σ exp(xᵢ - max(x)))

    # Convert energies to Boltzmann factors: -E/RT
    boltz_factors = []
    for grid in grids:
        boltz_factors.append(-grid.data / RT)

    # Stack for vectorized operations
    stacked = np.stack(boltz_factors, axis=0)  # (n_grids, nx, ny, nz)
    weights_arr = np.array(weights).reshape(-1, 1, 1, 1)

    # Log-sum-exp with weights
    max_val = np.max(stacked, axis=0)
    exp_sum = np.sum(weights_arr * np.exp(stacked - max_val), axis=0)

    # Convert back to energy
    log_avg = max_val + np.log(exp_sum)
    avg_energy = -RT * log_avg

    return avg_energy


def calculate_expected_density(
    n_probe_atoms: int,
    box_volume: float,
    n_frames: int,
    voxel_volume: float,
) -> float:
    """
    Calculate expected probe density per voxel.

    Parameters
    ----------
    n_probe_atoms : int
        Number of probe atoms in system
    box_volume : float
        Total box volume in Å³
    n_frames : int
        Number of trajectory frames
    voxel_volume : float
        Volume of single grid voxel in Å³

    Returns
    -------
    float
        Expected counts per voxel
    """
    # Probability of finding a probe atom in a voxel
    p_voxel = voxel_volume / box_volume

    # Expected count per frame
    expected_per_frame = n_probe_atoms * p_voxel

    # Total expected over all frames
    return expected_per_frame * n_frames


def normalize_grid(
    grid: Grid,
    expected_density: float,
) -> Grid:
    """
    Normalize density grid by expected value.

    Result is ρ/ρ₀ where ρ₀ is the expected bulk density.
    Values > 1 indicate higher than bulk density.

    Parameters
    ----------
    grid : Grid
        Raw count grid
    expected_density : float
        Expected counts per voxel

    Returns
    -------
    Grid
        Normalized density grid
    """
    if expected_density <= 0:
        raise ValueError("Expected density must be positive")

    normalized = grid.data / expected_density

    return Grid(
        data=normalized,
        origin=grid.origin,
        spacing=grid.spacing,
    )


def replica_average(
    replicas,
    probe: str,
    temperature: float = 300.0,
    grid_suffix: str = "_density.dx",
) -> Grid:
    """
    Average density grids across replicas for a specific probe.

    Parameters
    ----------
    replicas : sequence of Replica
        Replicas to average
    probe : str
        Probe name (e.g., "WAT", "CT", "OH")
    temperature : float
        Temperature for Boltzmann averaging
    grid_suffix : str
        Suffix for density grid files

    Returns
    -------
    Grid
        Averaged grid
    """
    grid_paths = []

    for replica in replicas:
        grid_path = replica.path / f"{probe}{grid_suffix}"
        if not grid_path.exists():
            log.warning(f"Grid not found: {grid_path}")
            continue
        grid_paths.append(grid_path)

    if not grid_paths:
        raise ValueError(f"No grids found for probe {probe}")

    return boltzmann_average(grid_paths, temperature=temperature)


# =============================================================================
# EnergyAction wrapper class
# =============================================================================


@register_action("energy")
class EnergyAction(Action):
    """
    Action for density-to-free-energy conversion.

    Reads density grids from a replica directory, converts each to
    a free-energy grid using ΔG = -RT ln(ρ/ρ₀), and writes
    ``{probe}_DG.dx`` files.

    Parameters (passed as kwargs to run())
    ---------------------------------------
    replica : Replica
        Replica whose density grids should be converted.
    temperature : float
        Temperature in Kelvin (default: 300).
    grids_subdir : str
        Subdirectory inside replica.path that holds ``*_density.dx``
        files (default: "grids").
    """

    name = "energy"
    description = "Convert density grids to free energy"

    def __init__(self, temperature: float = 300.0):
        super().__init__()
        self.temperature = temperature

    def run(
        self,
        trajectory=None,
        reference=None,
        output_dir: Path | None = None,
        *,
        replica=None,
        temperature: float | None = None,
        grids_subdir: str = "grids",
        **kwargs,
    ) -> ActionResult:
        """
        Convert density grids to energy for a replica.

        Parameters
        ----------
        trajectory : ignored
            Not used; kept for interface compatibility.
        reference : ignored
            Not used; kept for interface compatibility.
        output_dir : Path | None
            Where to write energy grids.  Defaults to
            ``replica.path / grids_subdir``.
        replica : Replica
            Replica with density grids (required).
        temperature : float | None
            Override instance temperature.
        grids_subdir : str
            Sub-directory holding density DX files.

        Returns
        -------
        ActionResult
            success=True with output_files and metadata.
        """
        if replica is None:
            return ActionResult(
                success=False,
                error="replica is required for EnergyAction.run()",
            )

        temp = temperature if temperature is not None else self.temperature

        # Locate density grids
        grids_path = replica.path / grids_subdir
        if not grids_path.exists():
            return ActionResult(
                success=False,
                error=f"No grids directory found: {grids_path}",
            )

        density_files = sorted(grids_path.glob("*_density.dx"))
        if not density_files:
            return ActionResult(
                success=False,
                error=f"No density grids found in {grids_path}",
            )

        out_dir = Path(output_dir) if output_dir else grids_path
        out_dir.mkdir(parents=True, exist_ok=True)

        energy_grids: dict[str, Grid] = {}
        output_files = []

        for density_file in density_files:
            probe_name = density_file.stem.replace("_density", "")

            density_grid = Grid.read_dx(str(density_file))
            energy_grid = density_to_free_energy(density_grid, temp)

            energy_file = out_dir / f"{probe_name}_DG.dx"
            energy_grid.write_dx(str(energy_file))

            energy_grids[probe_name] = energy_grid
            output_files.append(energy_file)

        result = EnergyResult(
            grids=energy_grids,
            temperature=temp,
            method="density_to_free_energy",
            n_replicas=1,
        )

        return ActionResult(
            success=True,
            output_files=output_files,
            metadata={
                "n_probes": len(energy_grids),
                "probe_names": list(energy_grids.keys()),
                "temperature": temp,
                "energy_result": result,
            },
        )
