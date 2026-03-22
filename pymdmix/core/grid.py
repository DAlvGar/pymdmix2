"""
3D Grid for Density and Volumetric Data
========================================

Grid class for storing and manipulating 3D volumetric data,
particularly density grids from MD trajectory analysis.

Supports OpenDX and XPLOR formats for visualization in VMD, PyMOL, etc.

Features:
- Count-based density calculation
- Free energy conversion
- Grid operations (expand, contract, subgrid extraction)
- Radial/spherical operations
- Protein masking
- Multiple file format support (DX, XPLOR)
"""

from __future__ import annotations

import gzip
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np
from numpy.typing import NDArray

if TYPE_CHECKING:
    pass

log = logging.getLogger(__name__)


@dataclass
class Grid:
    """
    3D grid for storing volumetric data.

    Attributes
    ----------
    data : NDArray[np.float64]
        3D numpy array of grid values (nx, ny, nz)
    origin : tuple[float, float, float]
        (x, y, z) coordinates of grid origin (corner)
    spacing : float | tuple[float, float, float]
        Grid spacing in Angstroms (uniform or per-axis)
    metadata : dict
        Optional metadata (probe name, type, etc.)

    Examples
    --------
    >>> # Create grid from structure coordinates
    >>> grid = Grid.from_bounds(min_coord, max_coord, spacing=0.5)
    >>>
    >>> # Add density counts
    >>> for coord in probe_coords:
    ...     grid.add_count(coord)
    >>>
    >>> # Normalize and save
    >>> density = grid.to_density(n_frames=1000)
    >>> density.write_dx("probe_density.dx")
    """

    data: NDArray[np.float64]
    origin: tuple[float, float, float]
    spacing: float | tuple[float, float, float]
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self):
        """Validate grid data."""
        if self.data.ndim != 3:
            raise ValueError(f"Grid data must be 3D, got {self.data.ndim}D")

        # Normalize spacing to tuple
        if isinstance(self.spacing, (int, float)):
            if self.spacing <= 0:
                raise ValueError(f"Grid spacing must be positive, got {self.spacing}")
        else:
            if any(s <= 0 for s in self.spacing):
                raise ValueError(f"Grid spacing must be positive, got {self.spacing}")

    @property
    def spacing_tuple(self) -> tuple[float, float, float]:
        """Get spacing as tuple."""
        if isinstance(self.spacing, (int, float)):
            return (float(self.spacing), float(self.spacing), float(self.spacing))
        return self.spacing

    @property
    def shape(self) -> tuple[int, int, int]:
        """Grid dimensions (nx, ny, nz)."""
        return self.data.shape  # type: ignore

    @property
    def dimensions(self) -> tuple[float, float, float]:
        """Physical dimensions in Angstroms."""
        sp = self.spacing_tuple
        return (
            self.shape[0] * sp[0],
            self.shape[1] * sp[1],
            self.shape[2] * sp[2],
        )

    @property
    def center(self) -> tuple[float, float, float]:
        """Grid center coordinates."""
        return tuple(self.origin[i] + self.dimensions[i] / 2 for i in range(3))  # type: ignore

    @property
    def extent(self) -> tuple[tuple[float, float], tuple[float, float], tuple[float, float]]:
        """Grid extent as ((xmin, xmax), (ymin, ymax), (zmin, zmax))."""
        return tuple((self.origin[i], self.origin[i] + self.dimensions[i]) for i in range(3))  # type: ignore

    # =========================================================================
    # Value Access
    # =========================================================================

    def get_value(self, x: float, y: float, z: float) -> float | None:
        """
        Get grid value at a specific coordinate.

        Parameters
        ----------
        x, y, z : float
            Coordinates in Angstroms

        Returns
        -------
        float or None
            Grid value at coordinate, or None if outside grid
        """
        try:
            idx = self.coord_to_index((x, y, z))
            if all(0 <= idx[i] < self.shape[i] for i in range(3)):
                return float(self.data[idx])
        except (IndexError, ValueError):
            pass
        return None

    def set_value(self, x: float, y: float, z: float, value: float) -> bool:
        """
        Set grid value at a specific coordinate.

        Returns True if successful, False if outside grid.
        """
        try:
            idx = self.coord_to_index((x, y, z))
            if all(0 <= idx[i] < self.shape[i] for i in range(3)):
                self.data[idx] = value
                return True
        except (IndexError, ValueError):
            pass
        return False

    def __getitem__(self, key: tuple[int, int, int]) -> float:
        """Get value by index."""
        return self.data[key]

    def __setitem__(self, key: tuple[int, int, int], value: float) -> None:
        """Set value by index."""
        self.data[key] = value

    # =========================================================================
    # Factory Methods
    # =========================================================================

    @classmethod
    def from_bounds(
        cls,
        min_coord: NDArray[np.float64] | tuple[float, float, float],
        max_coord: NDArray[np.float64] | tuple[float, float, float],
        spacing: float = 0.5,
        fill_value: float = 0.0,
    ) -> Grid:
        """
        Create grid from bounding box.

        Parameters
        ----------
        min_coord : array-like
            Minimum (x, y, z) coordinates
        max_coord : array-like
            Maximum (x, y, z) coordinates
        spacing : float
            Grid spacing in Angstroms
        fill_value : float
            Initial fill value

        Returns
        -------
        Grid
            Grid covering the bounding box
        """
        min_c = np.asarray(min_coord)
        max_c = np.asarray(max_coord)

        shape = tuple(int(np.ceil((max_c[i] - min_c[i]) / spacing)) + 1 for i in range(3))

        data = np.full(shape, fill_value, dtype=np.float64)

        return cls(
            data=data,
            origin=tuple(min_c),
            spacing=spacing,
        )

    @classmethod
    def from_structure(
        cls,
        coords: NDArray[np.float64],
        spacing: float = 0.5,
        padding: float = 5.0,
        fill_value: float = 0.0,
    ) -> Grid:
        """
        Create grid encompassing structure coordinates with padding.

        Parameters
        ----------
        coords : NDArray
            Atomic coordinates (n_atoms, 3)
        spacing : float
            Grid spacing in Angstroms
        padding : float
            Padding around structure in Angstroms
        fill_value : float
            Initial fill value

        Returns
        -------
        Grid
            Grid encompassing the structure
        """
        min_coord = coords.min(axis=0) - padding
        max_coord = coords.max(axis=0) + padding
        return cls.from_bounds(min_coord, max_coord, spacing, fill_value)

    def copy(self) -> Grid:
        """Create a deep copy of the grid."""
        return Grid(
            data=self.data.copy(),
            origin=self.origin,
            spacing=self.spacing,
            metadata=self.metadata.copy(),
        )

    # =========================================================================
    # Coordinate Conversion
    # =========================================================================

    def coord_to_index(
        self,
        coord: NDArray[np.float64] | tuple[float, float, float],
    ) -> tuple[int, int, int]:
        """
        Convert coordinate to grid index.

        Parameters
        ----------
        coord : array-like
            (x, y, z) coordinate

        Returns
        -------
        tuple[int, int, int]
            Grid indices (ix, iy, iz)
        """
        c = np.asarray(coord)
        sp = self.spacing_tuple
        idx = tuple(int((c[i] - self.origin[i]) / sp[i]) for i in range(3))
        return idx  # type: ignore

    def index_to_coord(
        self,
        index: tuple[int, int, int],
    ) -> tuple[float, float, float]:
        """
        Convert grid index to coordinate.

        Parameters
        ----------
        index : tuple[int, int, int]
            Grid indices (ix, iy, iz)

        Returns
        -------
        tuple[float, float, float]
            (x, y, z) coordinate at grid point
        """
        sp = self.spacing_tuple
        return tuple(self.origin[i] + index[i] * sp[i] for i in range(3))  # type: ignore

    def get_cartesian(self, indices: NDArray[np.int64]) -> NDArray[np.float64]:
        """Convert multiple indices to cartesian coordinates."""
        sp = np.array(self.spacing_tuple)
        origin = np.array(self.origin)
        return indices * sp + origin

    def is_inside(self, coord: NDArray[np.float64] | tuple) -> bool:
        """Check if coordinate is inside grid bounds."""
        try:
            idx = self.coord_to_index(coord)
            return all(0 <= idx[i] < self.shape[i] for i in range(3))
        except (IndexError, ValueError):
            return False

    # =========================================================================
    # Count Operations
    # =========================================================================

    def add_count(self, coord: NDArray[np.float64]) -> bool:
        """
        Add count at coordinate position.

        Parameters
        ----------
        coord : NDArray
            (x, y, z) coordinate

        Returns
        -------
        bool
            True if count was added, False if outside grid
        """
        if not self.is_inside(coord):
            return False
        idx = self.coord_to_index(coord)
        self.data[idx] += 1
        return True

    def add_counts_bulk(self, coords: NDArray[np.float64]) -> int:
        """
        Add counts for multiple coordinates efficiently.

        Parameters
        ----------
        coords : NDArray
            Coordinates (n_points, 3)

        Returns
        -------
        int
            Number of counts added (points inside grid)
        """
        sp = np.array(self.spacing_tuple)
        origin = np.array(self.origin)

        # Convert all coordinates to indices
        indices = ((coords - origin) / sp).astype(int)

        # Filter to valid indices
        valid = np.all((indices >= 0) & (indices < np.array(self.shape)), axis=1)
        valid_indices = indices[valid]

        # Add counts using numpy bincount for efficiency
        if len(valid_indices) > 0:
            flat_indices = np.ravel_multi_index(valid_indices.T, self.shape)
            counts = np.bincount(flat_indices, minlength=self.data.size)
            self.data += counts.reshape(self.shape)

        return int(valid.sum())

    # =========================================================================
    # Data Transformations
    # =========================================================================

    def to_density(self, n_frames: int) -> Grid:
        """
        Convert counts to density (counts per frame).

        Parameters
        ----------
        n_frames : int
            Number of frames used for counting

        Returns
        -------
        Grid
            New grid with density values
        """
        if n_frames <= 0:
            raise ValueError(f"n_frames must be positive, got {n_frames}")

        return Grid(
            data=self.data / n_frames,
            origin=self.origin,
            spacing=self.spacing,
            metadata={**self.metadata, "type": "density"},
        )

    def to_free_energy(
        self,
        temperature: float = 300.0,
        reference: str = "mean",
        mask_value: float = 999.0,
    ) -> Grid:
        """
        Convert density to free energy (ΔG = -RT ln(ρ/ρ0)).

        Parameters
        ----------
        temperature : float
            Temperature in Kelvin
        reference : str
            Reference density: "mean", "max", or "bulk"
        mask_value : float
            Value to use for zero/negative density regions

        Returns
        -------
        Grid
            Grid with free energy values in kcal/mol
        """
        R = 0.001987204  # kcal/(mol·K)
        RT = R * temperature

        # Get reference density
        positive_data = self.data[self.data > 0]
        if len(positive_data) == 0:
            raise ValueError("No positive density values found")

        if reference == "mean":
            rho0 = positive_data.mean()
        elif reference == "max":
            rho0 = positive_data.max()
        elif reference == "bulk":
            # Use median of outer shell
            rho0 = np.median(positive_data)
        else:
            rho0 = float(reference)

        # Calculate free energy
        with np.errstate(divide="ignore", invalid="ignore"):
            dg = -RT * np.log(self.data / rho0)

        # Mask invalid values
        dg = np.where(np.isfinite(dg), dg, mask_value)

        return Grid(
            data=dg,
            origin=self.origin,
            spacing=self.spacing,
            metadata={**self.metadata, "type": "energy"},
        )

    def count_to_dg(
        self,
        expected: float,
        temperature: float = 300.0,
        volume_correction: float = 1.0,
        mask_value: float = 999.0,
    ) -> Grid:
        """
        Convert counts to free energy with volume correction.

        ΔG = -RT ln(ρ/ρ_expected * correction)

        Parameters
        ----------
        expected : float
            Expected density (bulk)
        temperature : float
            Temperature in Kelvin
        volume_correction : float
            Volume correction factor
        mask_value : float
            Value for masked regions

        Returns
        -------
        Grid
            Free energy grid
        """
        R = 0.001987204  # kcal/(mol·K)
        RT = R * temperature

        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = (self.data / expected) * volume_correction
            dg = -RT * np.log(ratio)

        dg = np.where(np.isfinite(dg), dg, mask_value)
        dg = np.where(self.data > 0, dg, mask_value)

        return Grid(
            data=dg,
            origin=self.origin,
            spacing=self.spacing,
            metadata={**self.metadata, "type": "energy"},
        )

    def average_data(self, contract: bool = False) -> float:
        """Calculate average of all grid values."""
        if contract:
            return float(self.data[self.data != 0].mean())
        return float(self.data.mean())

    def update(self, data: NDArray[np.float64]) -> None:
        """Update grid data in-place."""
        if data.shape != self.shape:
            raise ValueError(f"Shape mismatch: {data.shape} vs {self.shape}")
        self.data = data

    # =========================================================================
    # Masking Operations
    # =========================================================================

    def mask_out(
        self,
        value: float,
        mask_value: float = 0.0,
        above: bool = True,
        include: bool = True,
    ) -> Grid:
        """
        Mask grid values based on threshold.

        Parameters
        ----------
        value : float
            Threshold value
        mask_value : float
            Value to set for masked points
        above : bool
            If True, mask values above threshold; if False, mask below
        include : bool
            If True, include threshold value in mask

        Returns
        -------
        Grid
            Masked grid
        """
        if above:
            if include:
                mask = self.data >= value
            else:
                mask = self.data > value
        else:
            if include:
                mask = self.data <= value
            else:
                mask = self.data < value

        new_data = self.data.copy()
        new_data[mask] = mask_value

        return Grid(
            data=new_data,
            origin=self.origin,
            spacing=self.spacing,
            metadata=self.metadata.copy(),
        )

    def cancel_points(
        self,
        point: tuple[float, float, float],
        cutoff: float,
        value: float = 999.0,
    ) -> Grid:
        """
        Set values within cutoff of a point to specified value.

        Parameters
        ----------
        point : tuple
            Center point (x, y, z)
        cutoff : float
            Distance cutoff
        value : float
            Value to set

        Returns
        -------
        Grid
            Modified grid
        """
        indices = self.get_radial_indices(cutoff, point)
        new_data = self.data.copy()

        for idx in indices:
            new_data[idx] = value

        return Grid(
            data=new_data,
            origin=self.origin,
            spacing=self.spacing,
            metadata=self.metadata.copy(),
        )

    def merge_delete_protein(
        self,
        coords: NDArray[np.float64],
        cutoff: float = 2.0,
        value: float = 999.0,
    ) -> Grid:
        """
        Mask grid points within cutoff of protein atoms.

        Parameters
        ----------
        coords : NDArray
            Protein atom coordinates (n_atoms, 3)
        cutoff : float
            Distance cutoff
        value : float
            Value to set for masked points

        Returns
        -------
        Grid
            Masked grid
        """
        new_data = self.data.copy()

        # For each grid point, check distance to nearest atom
        for ix in range(self.shape[0]):
            for iy in range(self.shape[1]):
                for iz in range(self.shape[2]):
                    grid_coord = np.array(self.index_to_coord((ix, iy, iz)))
                    distances = np.linalg.norm(coords - grid_coord, axis=1)
                    if distances.min() < cutoff:
                        new_data[ix, iy, iz] = value

        return Grid(
            data=new_data,
            origin=self.origin,
            spacing=self.spacing,
            metadata=self.metadata.copy(),
        )

    def merge_conserve_protein(
        self,
        coords: NDArray[np.float64],
        cutoff: float = 2.0,
        value: float = 999.0,
    ) -> Grid:
        """
        Mask grid points OUTSIDE cutoff of protein atoms.

        Parameters
        ----------
        coords : NDArray
            Protein atom coordinates (n_atoms, 3)
        cutoff : float
            Distance cutoff
        value : float
            Value to set for masked points

        Returns
        -------
        Grid
            Masked grid
        """
        new_data = self.data.copy()

        for ix in range(self.shape[0]):
            for iy in range(self.shape[1]):
                for iz in range(self.shape[2]):
                    grid_coord = np.array(self.index_to_coord((ix, iy, iz)))
                    distances = np.linalg.norm(coords - grid_coord, axis=1)
                    if distances.min() >= cutoff:
                        new_data[ix, iy, iz] = value

        return Grid(
            data=new_data,
            origin=self.origin,
            spacing=self.spacing,
            metadata=self.metadata.copy(),
        )

    # =========================================================================
    # Radial/Spherical Operations
    # =========================================================================

    def get_radial_indices(
        self,
        radius: float,
        center: tuple[float, float, float] | None = None,
        min_radius: float = 0.0,
    ) -> list[tuple[int, int, int]]:
        """
        Get grid indices within radius of center point.

        Parameters
        ----------
        radius : float
            Maximum radius
        center : tuple | None
            Center point. If None, uses grid center.
        min_radius : float
            Minimum radius (for shell selection)

        Returns
        -------
        list[tuple]
            List of (ix, iy, iz) indices
        """
        if center is None:
            center = self.center

        center_arr = np.array(center)
        indices = []

        # Calculate grid points to check
        sp = self.spacing_tuple
        max_grid_dist = int(np.ceil(radius / min(sp))) + 1
        center_idx = self.coord_to_index(center)

        for ix in range(
            max(0, center_idx[0] - max_grid_dist),
            min(self.shape[0], center_idx[0] + max_grid_dist + 1),
        ):
            for iy in range(
                max(0, center_idx[1] - max_grid_dist),
                min(self.shape[1], center_idx[1] + max_grid_dist + 1),
            ):
                for iz in range(
                    max(0, center_idx[2] - max_grid_dist),
                    min(self.shape[2], center_idx[2] + max_grid_dist + 1),
                ):
                    coord = np.array(self.index_to_coord((ix, iy, iz)))
                    dist = np.linalg.norm(coord - center_arr)
                    if min_radius <= dist <= radius:
                        indices.append((ix, iy, iz))

        return indices

    def get_sphere_values(
        self,
        center: tuple[float, float, float],
        radius: float,
    ) -> NDArray[np.float64]:
        """
        Get all grid values within sphere.

        Parameters
        ----------
        center : tuple
            Sphere center
        radius : float
            Sphere radius

        Returns
        -------
        NDArray
            Array of values within sphere
        """
        indices = self.get_radial_indices(radius, center)
        return np.array([self.data[idx] for idx in indices])

    def set_sphere_values(
        self,
        center: tuple[float, float, float],
        radius: float,
        value: float,
    ) -> None:
        """
        Set all grid values within sphere.

        Parameters
        ----------
        center : tuple
            Sphere center
        radius : float
            Sphere radius
        value : float
            Value to set
        """
        indices = self.get_radial_indices(radius, center)
        for idx in indices:
            self.data[idx] = value

    def get_radial_values(
        self,
        center: tuple[float, float, float],
        r_max: float,
        r_min: float = 0.0,
    ) -> NDArray[np.float64]:
        """Get values in radial shell."""
        indices = self.get_radial_indices(r_max, center, r_min)
        return np.array([self.data[idx] for idx in indices])

    def set_radial_values(
        self,
        r_max: float,
        r_min: float,
        value: float,
        center: tuple[float, float, float] | None = None,
    ) -> None:
        """Set values in radial shell."""
        indices = self.get_radial_indices(r_max, center, r_min)
        for idx in indices:
            self.data[idx] = value

    # =========================================================================
    # Grid Resizing
    # =========================================================================

    def expand(self, buffer: int, fill_value: float = 0.0) -> Grid:
        """
        Expand grid by adding buffer points on all sides.

        Parameters
        ----------
        buffer : int
            Number of grid points to add
        fill_value : float
            Value for new points

        Returns
        -------
        Grid
            Expanded grid
        """
        new_shape = tuple(s + 2 * buffer for s in self.shape)
        new_data = np.full(new_shape, fill_value, dtype=np.float64)

        # Copy original data to center
        new_data[
            buffer : buffer + self.shape[0],
            buffer : buffer + self.shape[1],
            buffer : buffer + self.shape[2],
        ] = self.data

        # Calculate new origin
        sp = self.spacing_tuple
        new_origin = tuple(self.origin[i] - buffer * sp[i] for i in range(3))

        return Grid(
            data=new_data,
            origin=new_origin,
            spacing=self.spacing,
            metadata=self.metadata.copy(),
        )

    def contract(self, buffer: int) -> Grid:
        """
        Contract grid by removing buffer points from all sides.

        Parameters
        ----------
        buffer : int
            Number of grid points to remove

        Returns
        -------
        Grid
            Contracted grid
        """
        if any(s <= 2 * buffer for s in self.shape):
            raise ValueError(f"Cannot contract by {buffer}: grid too small")

        new_data = self.data[
            buffer:-buffer,
            buffer:-buffer,
            buffer:-buffer,
        ].copy()

        sp = self.spacing_tuple
        new_origin = tuple(self.origin[i] + buffer * sp[i] for i in range(3))

        return Grid(
            data=new_data,
            origin=new_origin,
            spacing=self.spacing,
            metadata=self.metadata.copy(),
        )

    def trim(self, threshold: float = 0.0) -> Grid:
        """
        Trim grid to smallest box containing non-zero values.

        Parameters
        ----------
        threshold : float
            Values above this are considered non-zero

        Returns
        -------
        Grid
            Trimmed grid
        """
        # Find non-zero regions
        nonzero = np.where(self.data > threshold)

        if len(nonzero[0]) == 0:
            return self.copy()

        # Get bounding box
        min_idx = [arr.min() for arr in nonzero]
        max_idx = [arr.max() for arr in nonzero]

        new_data = self.data[
            min_idx[0] : max_idx[0] + 1,
            min_idx[1] : max_idx[1] + 1,
            min_idx[2] : max_idx[2] + 1,
        ].copy()

        sp = self.spacing_tuple
        new_origin = tuple(self.origin[i] + min_idx[i] * sp[i] for i in range(3))

        return Grid(
            data=new_data,
            origin=new_origin,
            spacing=self.spacing,
            metadata=self.metadata.copy(),
        )

    def take_subgrid_box(
        self,
        min_coord: tuple[float, float, float],
        max_coord: tuple[float, float, float],
    ) -> Grid:
        """
        Extract subgrid within coordinate bounds.

        Parameters
        ----------
        min_coord : tuple
            Minimum corner (x, y, z)
        max_coord : tuple
            Maximum corner (x, y, z)

        Returns
        -------
        Grid
            Subgrid
        """
        min_idx = self.coord_to_index(min_coord)
        max_idx = self.coord_to_index(max_coord)

        # Clamp to valid range
        min_idx = tuple(max(0, min(idx, s - 1)) for idx, s in zip(min_idx, self.shape))
        max_idx = tuple(max(0, min(idx, s - 1)) for idx, s in zip(max_idx, self.shape))

        new_data = self.data[
            min_idx[0] : max_idx[0] + 1,
            min_idx[1] : max_idx[1] + 1,
            min_idx[2] : max_idx[2] + 1,
        ].copy()

        new_origin = self.index_to_coord(min_idx)

        return Grid(
            data=new_data,
            origin=new_origin,
            spacing=self.spacing,
            metadata=self.metadata.copy(),
        )

    def take_subgrid_point(
        self,
        center: tuple[float, float, float],
        distance: float,
    ) -> Grid:
        """
        Extract cubic subgrid around a point.

        Parameters
        ----------
        center : tuple
            Center point (x, y, z)
        distance : float
            Half-width of cube

        Returns
        -------
        Grid
            Subgrid
        """
        min_coord = tuple(c - distance for c in center)
        max_coord = tuple(c + distance for c in center)
        return self.take_subgrid_box(min_coord, max_coord)

    # =========================================================================
    # Statistics
    # =========================================================================

    def min_index(self) -> tuple[int, int, int]:
        """Get index of minimum value."""
        idx = np.unravel_index(self.data.argmin(), self.shape)
        return tuple(idx)  # type: ignore

    def max_index(self) -> tuple[int, int, int]:
        """Get index of maximum value."""
        idx = np.unravel_index(self.data.argmax(), self.shape)
        return tuple(idx)  # type: ignore

    def min_coord(self) -> tuple[float, float, float]:
        """Get coordinate of minimum value."""
        return self.index_to_coord(self.min_index())

    def max_coord(self) -> tuple[float, float, float]:
        """Get coordinate of maximum value."""
        return self.index_to_coord(self.max_index())

    def percentile_cutoff(
        self,
        percentile: float,
        mask_value: float | None = None,
    ) -> float:
        """
        Get value at given percentile.

        Parameters
        ----------
        percentile : float
            Percentile (0-100)
        mask_value : float | None
            If set, exclude this value from calculation

        Returns
        -------
        float
            Value at percentile
        """
        data = self.data.flatten()
        if mask_value is not None:
            data = data[data != mask_value]
        return float(np.percentile(data, percentile))

    # =========================================================================
    # File I/O - DX Format
    # =========================================================================

    def write_dx(self, path: str | Path, gzip_compress: bool = False) -> None:
        """
        Write grid in OpenDX format.

        Parameters
        ----------
        path : str or Path
            Output file path
        gzip_compress : bool
            If True, write gzipped file
        """
        path = Path(path)
        nx, ny, nz = self.shape
        n_total = nx * ny * nz
        sp = self.spacing_tuple

        lines = []

        # Header
        if self.metadata:
            for key, val in self.metadata.items():
                lines.append(f"# {key}: {val}")
        lines.append("# OpenDX density grid")
        lines.append("# Generated by pyMDMix")
        lines.append(f"object 1 class gridpositions counts {nx} {ny} {nz}")
        lines.append(f"origin {self.origin[0]:.6f} {self.origin[1]:.6f} {self.origin[2]:.6f}")
        lines.append(f"delta {sp[0]:.6f} 0.000000 0.000000")
        lines.append(f"delta 0.000000 {sp[1]:.6f} 0.000000")
        lines.append(f"delta 0.000000 0.000000 {sp[2]:.6f}")
        lines.append(f"object 2 class gridconnections counts {nx} {ny} {nz}")
        lines.append(f"object 3 class array type double rank 0 items {n_total} data follows")

        # Data in Fortran order (z fastest, then y, then x)
        data_lines = []
        count = 0
        row = []
        for ix in range(nx):
            for iy in range(ny):
                for iz in range(nz):
                    row.append(f"{self.data[ix, iy, iz]:.6e}")
                    count += 1
                    if count % 3 == 0:
                        data_lines.append(" ".join(row))
                        row = []
        if row:
            data_lines.append(" ".join(row))

        lines.extend(data_lines)
        lines.append('attribute "dep" string "positions"')
        lines.append('object "density" class field')
        lines.append('component "positions" value 1')
        lines.append('component "connections" value 2')
        lines.append('component "data" value 3')

        content = "\n".join(lines) + "\n"

        if gzip_compress:
            with gzip.open(path, "wt") as f:
                f.write(content)
        else:
            path.write_text(content)

    @classmethod
    def read_dx(cls, path: str | Path) -> Grid:
        """
        Read grid from OpenDX format.

        Parameters
        ----------
        path : str or Path
            Input file path

        Returns
        -------
        Grid
            Loaded grid
        """
        path = Path(path)

        origin = None
        spacing = [None, None, None]
        shape = None
        data_lines: list[str] = []
        reading_data = False
        metadata = {}

        opener = gzip.open if str(path).endswith(".gz") else open

        with opener(path, "rt") as f:
            for line in f:
                line = line.strip()

                if line.startswith("#"):
                    # Parse metadata from comments
                    if ":" in line:
                        key, val = line[1:].split(":", 1)
                        metadata[key.strip()] = val.strip()
                    continue

                if not line:
                    continue

                if "gridpositions counts" in line:
                    parts = line.split()
                    shape = (int(parts[-3]), int(parts[-2]), int(parts[-1]))

                elif line.startswith("origin"):
                    parts = line.split()
                    origin = (float(parts[1]), float(parts[2]), float(parts[3]))

                elif line.startswith("delta"):
                    parts = line.split()
                    deltas = [float(parts[1]), float(parts[2]), float(parts[3])]
                    for i, d in enumerate(deltas):
                        if d > 0:
                            spacing[i] = d

                elif "data follows" in line:
                    reading_data = True
                    continue

                elif reading_data:
                    if (
                        line.startswith("attribute")
                        or line.startswith("object")
                        or line.startswith("component")
                    ):
                        break
                    data_lines.append(line)

        if origin is None or shape is None:
            raise ValueError(f"Could not parse DX file: {path}")

        # Determine spacing
        sp_values = [s for s in spacing if s is not None]
        if not sp_values:
            raise ValueError(f"Could not determine spacing from DX file: {path}")
        sp = sp_values[0]  # Assume uniform if only one found

        # Parse data
        data_str = " ".join(data_lines)
        data_flat = np.fromstring(data_str, sep=" ")

        expected_size = shape[0] * shape[1] * shape[2]
        if len(data_flat) != expected_size:
            raise ValueError(f"Data size mismatch: expected {expected_size}, got {len(data_flat)}")

        # Reshape (data is in C order after reading)
        data = data_flat.reshape(shape, order="C")

        return cls(data=data, origin=origin, spacing=sp, metadata=metadata)

    # =========================================================================
    # File I/O - XPLOR Format
    # =========================================================================

    def write_xplor(
        self,
        path: str | Path,
        header: str = "",
        gzip_compress: bool = False,
    ) -> None:
        """
        Write grid in XPLOR/CNS format.

        Parameters
        ----------
        path : str or Path
            Output file path
        header : str
            Header text
        gzip_compress : bool
            If True, write gzipped file
        """
        path = Path(path)
        nx, ny, nz = self.shape
        sp = self.spacing_tuple

        lines = []

        # Header
        lines.append("")
        lines.append("       2 !NTITLE")
        lines.append(f" REMARKS {header if header else 'XPLOR grid from pyMDMix'}")
        lines.append(" REMARKS")

        # Grid info
        lines.append(f"{nx:8d}{0:8d}{nx - 1:8d}{ny:8d}{0:8d}{ny - 1:8d}{nz:8d}{0:8d}{nz - 1:8d}")

        # Cell dimensions (orthogonal)
        cell_a = nx * sp[0]
        cell_b = ny * sp[1]
        cell_c = nz * sp[2]
        lines.append(
            f"{cell_a:12.5E}{cell_b:12.5E}{cell_c:12.5E}{90.0:12.5E}{90.0:12.5E}{90.0:12.5E}"
        )
        lines.append("ZYX")

        # Data
        for iz in range(nz):
            lines.append(f"{iz:8d}")
            count = 0
            row = []
            for iy in range(ny):
                for ix in range(nx):
                    row.append(f"{self.data[ix, iy, iz]:12.5E}")
                    count += 1
                    if count % 6 == 0:
                        lines.append("".join(row))
                        row = []
            if row:
                lines.append("".join(row))

        lines.append(f"{-9999:8d}")

        # Statistics
        avg = self.data.mean()
        std = self.data.std()
        lines.append(f"{avg:12.4E}{std:12.4E}")

        content = "\n".join(lines) + "\n"

        if gzip_compress:
            with gzip.open(path, "wt") as f:
                f.write(content)
        else:
            path.write_text(content)

    @classmethod
    def read_xplor(cls, path: str | Path) -> Grid:
        """
        Read grid from XPLOR/CNS format.

        Parameters
        ----------
        path : str or Path
            Input file path

        Returns
        -------
        Grid
            Loaded grid
        """
        path = Path(path)

        opener = gzip.open if str(path).endswith(".gz") else open

        with opener(path, "rt") as f:
            lines = f.readlines()

        # Parse header
        idx = 0
        while idx < len(lines) and "NTITLE" not in lines[idx]:
            idx += 1
        idx += 1

        # Skip title lines
        while idx < len(lines) and "REMARKS" in lines[idx]:
            idx += 1

        # Parse grid dimensions
        dim_line = lines[idx].strip()
        parts = dim_line.split()
        nx = int(parts[0])
        ny = int(parts[3])
        nz = int(parts[6])
        idx += 1

        # Parse cell dimensions
        cell_line = lines[idx].strip()
        cell_parts = cell_line.split()
        cell_a = float(cell_parts[0])
        cell_b = float(cell_parts[1])
        cell_c = float(cell_parts[2])
        idx += 1

        # Skip ZYX line
        idx += 1

        # Parse data
        data = np.zeros((nx, ny, nz), dtype=np.float64)

        iz = 0
        while iz < nz:
            # Read section header
            section = int(lines[idx].strip())
            if section == -9999:
                break
            idx += 1

            # Read section data
            values = []
            while len(values) < nx * ny:
                line = lines[idx].strip()
                # Parse 6 values per line
                for i in range(0, len(line), 12):
                    try:
                        values.append(float(line[i : i + 12]))
                    except ValueError:
                        pass
                idx += 1

            # Fill data array
            vi = 0
            for iy in range(ny):
                for ix in range(nx):
                    data[ix, iy, iz] = values[vi]
                    vi += 1

            iz += 1

        # Calculate spacing
        spacing = (cell_a / nx, cell_b / ny, cell_c / nz)

        return cls(
            data=data,
            origin=(0.0, 0.0, 0.0),  # XPLOR typically uses 0 origin
            spacing=spacing[0] if spacing[0] == spacing[1] == spacing[2] else spacing,
        )

    # =========================================================================
    # Arithmetic Operations
    # =========================================================================

    def __add__(self, other: Grid | float) -> Grid:
        """Add two grids or add scalar."""
        if isinstance(other, Grid):
            if self.shape != other.shape:
                raise ValueError(f"Shape mismatch: {self.shape} vs {other.shape}")
            return Grid(
                data=self.data + other.data,
                origin=self.origin,
                spacing=self.spacing,
                metadata=self.metadata.copy(),
            )
        else:
            return Grid(
                data=self.data + other,
                origin=self.origin,
                spacing=self.spacing,
                metadata=self.metadata.copy(),
            )

    def __sub__(self, other: Grid | float) -> Grid:
        """Subtract grids or scalar."""
        if isinstance(other, Grid):
            if self.shape != other.shape:
                raise ValueError(f"Shape mismatch: {self.shape} vs {other.shape}")
            return Grid(
                data=self.data - other.data,
                origin=self.origin,
                spacing=self.spacing,
                metadata=self.metadata.copy(),
            )
        else:
            return Grid(
                data=self.data - other,
                origin=self.origin,
                spacing=self.spacing,
                metadata=self.metadata.copy(),
            )

    def __mul__(self, scalar: float) -> Grid:
        """Multiply grid by scalar."""
        return Grid(
            data=self.data * scalar,
            origin=self.origin,
            spacing=self.spacing,
            metadata=self.metadata.copy(),
        )

    def __truediv__(self, scalar: float) -> Grid:
        """Divide grid by scalar."""
        return Grid(
            data=self.data / scalar,
            origin=self.origin,
            spacing=self.spacing,
            metadata=self.metadata.copy(),
        )

    def __repr__(self) -> str:
        return (
            f"Grid(shape={self.shape}, origin={self.origin}, "
            f"spacing={self.spacing}, range=[{self.data.min():.3f}, {self.data.max():.3f}])"
        )


# =============================================================================
# Utility Functions
# =============================================================================


def average_grids(
    grids: list[Grid],
    boltzmann: bool = False,
    temperature: float = 300.0,
) -> Grid:
    """
    Average multiple grids.

    Parameters
    ----------
    grids : list[Grid]
        Grids to average (must have same shape)
    boltzmann : bool
        Use Boltzmann-weighted averaging
    temperature : float
        Temperature for Boltzmann weighting

    Returns
    -------
    Grid
        Averaged grid
    """
    if not grids:
        raise ValueError("No grids to average")

    ref = grids[0]
    for g in grids[1:]:
        if g.shape != ref.shape:
            raise ValueError("All grids must have same shape")

    if boltzmann:
        R = 0.001987204  # kcal/(mol·K)
        RT = R * temperature

        # Boltzmann weights
        weights = []
        for g in grids:
            w = np.exp(-g.data / RT)
            weights.append(w)

        total_weight = sum(weights)
        data = sum(g.data * w for g, w in zip(grids, weights)) / total_weight
    else:
        data = np.mean([g.data for g in grids], axis=0)

    return Grid(
        data=data,
        origin=ref.origin,
        spacing=ref.spacing,
        metadata={"averaged_from": len(grids)},
    )


def minimum_grids(grids: list[Grid]) -> Grid:
    """Element-wise minimum of grids."""
    if not grids:
        raise ValueError("No grids")
    ref = grids[0]
    data = np.minimum.reduce([g.data for g in grids])
    return Grid(data=data, origin=ref.origin, spacing=ref.spacing)


def maximum_grids(grids: list[Grid]) -> Grid:
    """Element-wise maximum of grids."""
    if not grids:
        raise ValueError("No grids")
    ref = grids[0]
    data = np.maximum.reduce([g.data for g in grids])
    return Grid(data=data, origin=ref.origin, spacing=ref.spacing)
