"""
Binding Hotspot Detection
=========================

Detect and cluster binding hotspots from density grids.

Hotspots are regions of high probe density that indicate
favorable binding sites on the protein surface.

The algorithm:
1. Load density grid(s)
2. Convert to free energy (ΔG = -RT ln(ρ/ρ0))
3. Find grid points below energy threshold
4. Cluster nearby points using hierarchical clustering
5. Rank hotspots by energy

Examples
--------
>>> from pymdmix.analysis import HotspotAction
>>> from pymdmix.core import Grid
>>>
>>> density = Grid.read_dx("probe_density.dx")
>>> action = HotspotAction()
>>> result = action.run(
...     grids={"OH": density},
...     energy_cutoff=-0.5,
...     cluster_distance=2.0,
... )
"""

from __future__ import annotations

import json
import logging
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

try:
    from scipy.cluster.hierarchy import fclusterdata
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False

from pymdmix.analysis.base import Action, ActionResult, register_action
from pymdmix.core.grid import Grid

log = logging.getLogger(__name__)


@dataclass
class Hotspot:
    """
    A binding hotspot detected from density analysis.

    Attributes
    ----------
    id : int
        Hotspot identifier
    probe : str
        Probe type (e.g., "OH", "CT")
    centroid : tuple[float, float, float]
        Center coordinates (minimum energy point)
    energy : float
        Hotspot energy in kcal/mol
    volume : float
        Hotspot volume in Å³
    n_points : int
        Number of grid points in cluster
    coords : NDArray
        All coordinates in the hotspot
    energies : NDArray
        Energy at each coordinate
    """
    id: int
    probe: str
    centroid: tuple[float, float, float]
    energy: float
    volume: float
    n_points: int
    coords: NDArray = field(repr=False)
    energies: NDArray = field(repr=False)

    @property
    def min_energy(self) -> float:
        """Minimum energy in the hotspot."""
        return float(self.energies.min())

    @property
    def mean_energy(self) -> float:
        """Mean energy in the hotspot."""
        return float(self.energies.mean())

    @property
    def extent(self) -> tuple[float, float, float]:
        """Extent in x, y, z dimensions."""
        return tuple(self.coords.max(axis=0) - self.coords.min(axis=0))

    def to_dict(self) -> dict:
        """Convert to dictionary for serialization."""
        return {
            "id": self.id,
            "probe": self.probe,
            "centroid": list(self.centroid),
            "energy": self.energy,
            "volume": self.volume,
            "n_points": self.n_points,
            "min_energy": self.min_energy,
            "mean_energy": self.mean_energy,
            "extent": list(self.extent),
        }


@register_action("hotspots")
class HotspotAction(Action):
    """
    Detect binding hotspots from density grids.

    Parameters (in run())
    ---------------------
    grids : dict[str, Grid]
        Probe name -> density grid
    energy_cutoff : float
        Energy threshold in kcal/mol (default: -0.5)
        Only points below this are considered
    cluster_distance : float
        Maximum distance for clustering in Å (default: 2.0)
    min_points : int
        Minimum points per cluster (default: 3)
    temperature : float
        Temperature for free energy (default: 300.0 K)
    output_dir : Path | None
        Output directory
    output_prefix : str
        Prefix for output files

    Outputs
    -------
    - {prefix}hotspots.json : Hotspot data
    - {prefix}hotspots.pdb : PDB with hotspot coordinates
    - {prefix}hotspots_summary.txt : Human-readable summary
    """

    name = "hotspots"
    description = "Detect binding hotspots from density grids"

    def run(
        self,
        grids: dict[str, Grid],
        energy_cutoff: float = -0.5,
        cluster_distance: float = 2.0,
        min_points: int = 3,
        temperature: float = 300.0,
        output_dir: Path | None = None,
        output_prefix: str = "",
        trajectory=None,  # Not used but required by base
        reference=None,
        **kwargs,
    ) -> ActionResult:
        """Execute hotspot detection."""

        if not HAS_SCIPY:
            return ActionResult(
                success=False,
                error="scipy is required for hotspot detection",
            )

        output_dir = Path(output_dir) if output_dir else Path.cwd()
        output_dir.mkdir(parents=True, exist_ok=True)

        self.log.info(f"Detecting hotspots from {len(grids)} grids")
        self.log.info(f"Energy cutoff: {energy_cutoff} kcal/mol")
        self.log.info(f"Cluster distance: {cluster_distance} Å")

        all_hotspots: list[Hotspot] = []
        hotspot_id = 0

        for probe_name, density_grid in grids.items():
            self.log.info(f"Processing probe: {probe_name}")

            # Convert density to free energy
            dg_grid = density_grid.to_free_energy(temperature)

            # Find points below threshold
            mask = dg_grid.data < energy_cutoff
            n_below = mask.sum()

            if n_below == 0:
                self.log.warning(f"No points below {energy_cutoff} kcal/mol for {probe_name}")
                continue

            self.log.info(f"  {n_below} points below threshold")

            # Get coordinates and energies of favorable points
            indices = np.argwhere(mask)
            coords = np.array([
                dg_grid.index_to_coord(tuple(idx))
                for idx in indices
            ])
            energies = dg_grid.data[mask]

            # Cluster points
            if len(coords) < min_points:
                self.log.warning(f"Too few points ({len(coords)}) for clustering")
                continue

            # Hierarchical clustering
            clusters = fclusterdata(
                coords,
                t=cluster_distance,
                criterion='distance',
                method='complete',
            )

            n_clusters = clusters.max()
            self.log.info(f"  Found {n_clusters} clusters")

            # Process each cluster
            for cluster_id in range(1, n_clusters + 1):
                cluster_mask = clusters == cluster_id
                cluster_coords = coords[cluster_mask]
                cluster_energies = energies[cluster_mask]

                if len(cluster_coords) < min_points:
                    continue

                # Find minimum energy point as centroid
                min_idx = cluster_energies.argmin()
                centroid = tuple(cluster_coords[min_idx])

                # Calculate volume (n_points * voxel_volume)
                voxel_volume = density_grid.spacing ** 3
                volume = len(cluster_coords) * voxel_volume

                # Volume-weighted energy
                energy = float(cluster_energies.mean())

                hotspot = Hotspot(
                    id=hotspot_id,
                    probe=probe_name,
                    centroid=centroid,
                    energy=energy,
                    volume=volume,
                    n_points=len(cluster_coords),
                    coords=cluster_coords,
                    energies=cluster_energies,
                )
                all_hotspots.append(hotspot)
                hotspot_id += 1

        # Sort by energy
        all_hotspots.sort(key=lambda h: h.energy)

        self.log.info(f"Total hotspots detected: {len(all_hotspots)}")

        # Save outputs
        output_files = []

        # JSON output
        json_path = output_dir / f"{output_prefix}hotspots.json"
        json_data = {
            "n_hotspots": len(all_hotspots),
            "energy_cutoff": energy_cutoff,
            "cluster_distance": cluster_distance,
            "temperature": temperature,
            "hotspots": [h.to_dict() for h in all_hotspots],
        }
        with open(json_path, 'w') as f:
            json.dump(json_data, f, indent=2)
        output_files.append(json_path)

        # PDB output
        pdb_path = output_dir / f"{output_prefix}hotspots.pdb"
        self._write_hotspots_pdb(all_hotspots, pdb_path)
        output_files.append(pdb_path)

        # Summary
        summary_path = output_dir / f"{output_prefix}hotspots_summary.txt"
        self._write_summary(all_hotspots, summary_path, energy_cutoff, cluster_distance)
        output_files.append(summary_path)

        return ActionResult(
            success=True,
            output_files=output_files,
            metadata={
                "n_hotspots": len(all_hotspots),
                "probes": list(grids.keys()),
                "energy_cutoff": energy_cutoff,
            },
        )

    def _write_hotspots_pdb(self, hotspots: list[Hotspot], path: Path) -> None:
        """Write hotspots as PDB file with B-factors as energy."""
        with open(path, 'w') as f:
            f.write("REMARK  Binding hotspots from pyMDMix\n")
            f.write("REMARK  B-factor column contains energy (kcal/mol)\n")

            atom_num = 1
            for hs in hotspots:
                x, y, z = hs.centroid
                energy = hs.energy
                probe = hs.probe[:3].upper()

                # Write as HETATM
                f.write(
                    f"HETATM{atom_num:5d}  C   {probe} X{hs.id:4d}    "
                    f"{x:8.3f}{y:8.3f}{z:8.3f}"
                    f"  1.00{energy:6.2f}           C\n"
                )
                atom_num += 1

            f.write("END\n")

    def _write_summary(
        self,
        hotspots: list[Hotspot],
        path: Path,
        energy_cutoff: float,
        cluster_distance: float,
    ) -> None:
        """Write human-readable summary."""
        with open(path, 'w') as f:
            f.write("Hotspot Detection Summary\n")
            f.write("=========================\n\n")
            f.write(f"Energy cutoff: {energy_cutoff} kcal/mol\n")
            f.write(f"Cluster distance: {cluster_distance} Å\n")
            f.write(f"Total hotspots: {len(hotspots)}\n\n")

            # Group by probe
            by_probe: dict[str, list[Hotspot]] = {}
            for hs in hotspots:
                by_probe.setdefault(hs.probe, []).append(hs)

            for probe, probe_hotspots in sorted(by_probe.items()):
                f.write(f"\n{probe} ({len(probe_hotspots)} hotspots)\n")
                f.write("-" * 50 + "\n")
                f.write(f"{'ID':>4} {'Energy':>8} {'Volume':>8} {'Centroid':<30}\n")

                for hs in probe_hotspots:
                    centroid_str = f"({hs.centroid[0]:.1f}, {hs.centroid[1]:.1f}, {hs.centroid[2]:.1f})"
                    f.write(f"{hs.id:>4} {hs.energy:>8.2f} {hs.volume:>8.1f} {centroid_str:<30}\n")

    def validate(self, trajectory=None, **kwargs) -> list[str]:
        """Validate inputs."""
        errors = []

        if not HAS_SCIPY:
            errors.append("scipy is required for hotspot detection")

        grids = kwargs.get('grids')
        if not grids:
            errors.append("grids dictionary is required")

        return errors


def detect_hotspots(
    grids: dict[str, Grid],
    energy_cutoff: float = -0.5,
    cluster_distance: float = 2.0,
    min_points: int = 3,
) -> list[Hotspot]:
    """
    Convenience function to detect hotspots.

    Parameters
    ----------
    grids : dict[str, Grid]
        Probe name -> density grid
    energy_cutoff : float
        Energy threshold in kcal/mol
    cluster_distance : float
        Maximum distance for clustering
    min_points : int
        Minimum points per cluster

    Returns
    -------
    list[Hotspot]
        Detected hotspots sorted by energy
    """
    action = HotspotAction()
    result = action.run(
        grids=grids,
        energy_cutoff=energy_cutoff,
        cluster_distance=cluster_distance,
        min_points=min_points,
    )

    if not result.success:
        raise RuntimeError(result.error)

    # Load hotspots from output JSON
    json_path = result.output_files[0]
    with open(json_path) as f:
        data = json.load(f)

    return data["hotspots"]
