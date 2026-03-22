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
        trajectory=None,
        reference=None,
        output_dir: Path | None = None,
        grids: dict[str, Grid] | None = None,
        energy_cutoff: float = -0.5,
        cluster_distance: float = 2.0,
        min_points: int = 3,
        temperature: float = 300.0,
        output_prefix: str = "",
        **kwargs,
    ) -> ActionResult:
        """Execute hotspot detection."""

        if grids is None:
            return ActionResult(
                success=False,
                error="grids must be provided to HotspotAction.run()",
            )

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
            coords = np.array([dg_grid.index_to_coord(tuple(idx)) for idx in indices])
            energies = dg_grid.data[mask]

            # Cluster points
            if len(coords) < min_points:
                self.log.warning(f"Too few points ({len(coords)}) for clustering")
                continue

            # Hierarchical clustering
            clusters = fclusterdata(
                coords,
                t=cluster_distance,
                criterion="distance",
                method="complete",
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
                spx, spy, spz = density_grid.spacing_tuple
                voxel_volume = spx * spy * spz
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
        with open(json_path, "w") as f:
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
        with open(path, "w") as f:
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
        with open(path, "w") as f:
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
                    centroid_str = (
                        f"({hs.centroid[0]:.1f}, {hs.centroid[1]:.1f}, {hs.centroid[2]:.1f})"
                    )
                    f.write(f"{hs.id:>4} {hs.energy:>8.2f} {hs.volume:>8.1f} {centroid_str:<30}\n")

    def validate(self, trajectory=None, **kwargs) -> list[str]:
        """Validate inputs."""
        errors = []

        if not HAS_SCIPY:
            errors.append("scipy is required for hotspot detection")

        grids = kwargs.get("grids")
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

    return data["hotspots"]  # type: ignore[no-any-return]


# =============================================================================
# HotSpotSet - Collection with clustering and filtering
# =============================================================================


@dataclass
class HotSpotSet:
    """
    Collection of hotspots with clustering and filtering capabilities.

    Replaces legacy HotSpotSet class with modern implementation.

    Attributes
    ----------
    probe : str
        Probe type for this set (e.g., "OH")
    name : str
        Set identifier
    hotspots : list[Hotspot]
        List of hotspots in the set
    info : str
        Description or metadata

    Examples
    --------
    >>> hs_set = HotSpotSet(probe="OH", name="hydroxyl_sites")
    >>> hs_set.add_hotspots(hotspots)
    >>> filtered = hs_set.prune_by_energy(-1.0)
    >>> clusters = hs_set.cluster(cutoff=2.5)
    """

    probe: str = ""
    name: str = ""
    hotspots: list[Hotspot] = field(default_factory=list)
    info: str = ""

    # Cached distance matrix and clustering
    _distance_matrix: NDArray | None = field(default=None, repr=False)
    _cluster_labels: NDArray | None = field(default=None, repr=False)

    def __post_init__(self):
        # Sort hotspots by energy (ascending)
        self.hotspots.sort(key=lambda h: h.energy)

    def __len__(self) -> int:
        return len(self.hotspots)

    def __iter__(self):
        return iter(self.hotspots)

    def __getitem__(self, idx: int) -> Hotspot:
        return self.hotspots[idx]

    @property
    def n_hotspots(self) -> int:
        """Number of hotspots in set."""
        return len(self.hotspots)

    @property
    def centroids(self) -> NDArray:
        """Array of hotspot centroid coordinates."""
        if not self.hotspots:
            return np.empty((0, 3))
        return np.array([h.centroid for h in self.hotspots])

    @property
    def energies(self) -> NDArray:
        """Array of hotspot energies."""
        return np.array([h.energy for h in self.hotspots])

    @property
    def volumes(self) -> NDArray:
        """Array of hotspot volumes."""
        return np.array([h.volume for h in self.hotspots])

    def add_hotspots(self, hotspots: Hotspot | list[Hotspot]) -> None:
        """
        Add hotspot(s) to the set.

        Parameters
        ----------
        hotspots : Hotspot or list[Hotspot]
            Hotspot(s) to add
        """
        if isinstance(hotspots, list):
            self.hotspots.extend(hotspots)
        else:
            self.hotspots.append(hotspots)

        # Re-sort and invalidate cache
        self.hotspots.sort(key=lambda h: h.energy)
        self._distance_matrix = None
        self._cluster_labels = None

    def get_by_id(self, ids: int | list[int]) -> list[Hotspot]:
        """
        Get hotspots by ID.

        Parameters
        ----------
        ids : int or list[int]
            Hotspot ID(s) to retrieve

        Returns
        -------
        list[Hotspot]
            Matching hotspots
        """
        if isinstance(ids, int):
            ids = [ids]
        return [h for h in self.hotspots if h.id in ids]

    def compute_distance_matrix(self) -> NDArray:
        """
        Compute pairwise distance matrix between hotspot centroids.

        Returns
        -------
        NDArray
            Square distance matrix (n x n)
        """
        if self._distance_matrix is not None:
            return self._distance_matrix

        if len(self.hotspots) < 2:
            return np.array([[]])

        from scipy.spatial import distance

        centroids = self.centroids
        condensed = distance.pdist(centroids)
        self._distance_matrix = distance.squareform(condensed)
        return self._distance_matrix

    def cluster(self, cutoff: float = 2.5, method: str = "average") -> NDArray:
        """
        Cluster hotspots by distance.

        Parameters
        ----------
        cutoff : float
            Distance cutoff for clustering (Å)
        method : str
            Linkage method ('average', 'single', 'complete', 'ward')

        Returns
        -------
        NDArray
            Cluster labels for each hotspot (1-indexed)
        """
        if not HAS_SCIPY:
            raise ImportError("scipy is required for clustering")

        if len(self.hotspots) < 2:
            return np.array([1] if self.hotspots else [])

        from scipy.cluster.hierarchy import fcluster, linkage
        from scipy.spatial.distance import pdist

        centroids = self.centroids
        condensed = pdist(centroids)
        Z = linkage(condensed, method=method)
        self._cluster_labels = fcluster(Z, t=cutoff, criterion="distance")

        return self._cluster_labels

    @property
    def n_clusters(self) -> int:
        """Number of clusters (after clustering)."""
        if self._cluster_labels is None:
            return 0
        return len(np.unique(self._cluster_labels))

    def get_cluster_representatives(self, cutoff: float = 2.5) -> HotSpotSet:
        """
        Get one representative hotspot per cluster (lowest energy).

        Parameters
        ----------
        cutoff : float
            Distance cutoff for clustering

        Returns
        -------
        HotSpotSet
            New set with one hotspot per cluster
        """
        labels = self.cluster(cutoff=cutoff)

        representatives = []
        for cluster_id in np.unique(labels):
            cluster_hotspots = [
                self.hotspots[i] for i, label in enumerate(labels) if label == cluster_id
            ]
            # Take lowest energy hotspot in cluster
            best = min(cluster_hotspots, key=lambda h: h.energy)
            representatives.append(best)

        return HotSpotSet(
            probe=self.probe,
            name=f"{self.name}_clustered",
            hotspots=representatives,
            info=f"Clustered with cutoff={cutoff}Å",
        )

    def prune_by_energy(self, max_energy: float) -> HotSpotSet:
        """
        Remove hotspots above energy threshold.

        Parameters
        ----------
        max_energy : float
            Maximum allowed energy (kcal/mol)

        Returns
        -------
        HotSpotSet
            New set with filtered hotspots
        """
        filtered = [h for h in self.hotspots if h.energy <= max_energy]
        return HotSpotSet(
            probe=self.probe,
            name=f"{self.name}_E<{max_energy}",
            hotspots=filtered,
            info=f"Pruned by energy <= {max_energy}",
        )

    def prune_by_volume(self, min_volume: float) -> HotSpotSet:
        """
        Remove hotspots below volume threshold.

        Parameters
        ----------
        min_volume : float
            Minimum required volume (Å³)

        Returns
        -------
        HotSpotSet
            New set with filtered hotspots
        """
        filtered = [h for h in self.hotspots if h.volume >= min_volume]
        return HotSpotSet(
            probe=self.probe,
            name=f"{self.name}_V>{min_volume}",
            hotspots=filtered,
            info=f"Pruned by volume >= {min_volume}",
        )

    def prune_by_n_points(self, min_points: int) -> HotSpotSet:
        """
        Remove hotspots with too few points.

        Parameters
        ----------
        min_points : int
            Minimum required grid points

        Returns
        -------
        HotSpotSet
            New set with filtered hotspots
        """
        filtered = [h for h in self.hotspots if h.n_points >= min_points]
        return HotSpotSet(
            probe=self.probe,
            name=f"{self.name}_N>{min_points}",
            hotspots=filtered,
            info=f"Pruned by n_points >= {min_points}",
        )

    def find_nearest(
        self,
        coord: tuple[float, float, float] | NDArray,
        max_distance: float = 5.0,
    ) -> Hotspot | None:
        """
        Find hotspot nearest to a coordinate.

        Parameters
        ----------
        coord : tuple or NDArray
            3D coordinate
        max_distance : float
            Maximum allowed distance (returns None if exceeded)

        Returns
        -------
        Hotspot | None
            Nearest hotspot or None if none within max_distance
        """
        if not self.hotspots:
            return None

        coord = np.asarray(coord)
        distances = np.linalg.norm(self.centroids - coord, axis=1)
        min_idx = distances.argmin()
        min_dist = distances[min_idx]

        if min_dist <= max_distance:
            return self.hotspots[int(min_idx)]
        return None

    def to_pdb(self, path: str | Path, only_centroids: bool = True) -> None:
        """
        Write hotspots to PDB file.

        Parameters
        ----------
        path : str or Path
            Output file path
        only_centroids : bool
            If True, write only centroids. If False, write all coordinates.
        """
        path = Path(path)

        with open(path, "w") as f:
            f.write(f"REMARK  HotSpotSet: {self.name}\n")
            f.write(f"REMARK  Probe: {self.probe}\n")
            f.write(f"REMARK  N hotspots: {len(self.hotspots)}\n")

            atom_num = 1
            for hs in self.hotspots:
                if only_centroids:
                    coords_to_write = [hs.centroid]
                else:
                    coords_to_write = hs.coords.tolist()

                for coord in coords_to_write:
                    x, y, z = coord
                    probe = hs.probe[:3].upper()
                    f.write(
                        f"HETATM{atom_num:5d}  C   {probe} X{hs.id:4d}    "
                        f"{x:8.3f}{y:8.3f}{z:8.3f}"
                        f"  1.00{hs.energy:6.2f}           C\n"
                    )
                    atom_num += 1

            f.write("END\n")

    def to_json(self, path: str | Path) -> None:
        """
        Write hotspots to JSON file.

        Parameters
        ----------
        path : str or Path
            Output file path
        """
        path = Path(path)

        data = {
            "probe": self.probe,
            "name": self.name,
            "info": self.info,
            "n_hotspots": len(self.hotspots),
            "hotspots": [h.to_dict() for h in self.hotspots],
        }

        with open(path, "w") as f:
            json.dump(data, f, indent=2)

    def summary(self) -> str:
        """Get human-readable summary."""
        lines = [
            f"HotSpotSet: {self.name}",
            f"Probe: {self.probe}",
            f"N hotspots: {len(self.hotspots)}",
            "-" * 40,
        ]

        if self.hotspots:
            lines.append(
                f"Energy range: {self.energies.min():.2f} to {self.energies.max():.2f} kcal/mol"
            )
            lines.append(f"Volume range: {self.volumes.min():.1f} to {self.volumes.max():.1f} Å³")
            lines.append("")
            lines.append(f"{'ID':>4} {'Energy':>8} {'Volume':>8} {'Points':>6}")

            for hs in self.hotspots[:10]:  # Show first 10
                lines.append(f"{hs.id:>4} {hs.energy:>8.2f} {hs.volume:>8.1f} {hs.n_points:>6}")

            if len(self.hotspots) > 10:
                lines.append(f"... and {len(self.hotspots) - 10} more")

        return "\n".join(lines)

    def __str__(self) -> str:
        return f"HotSpotSet({self.name}, {self.probe}, n={len(self.hotspots)})"

    def __repr__(self) -> str:
        return (
            f"HotSpotSet(probe={self.probe!r}, name={self.name!r}, n_hotspots={len(self.hotspots)})"
        )
