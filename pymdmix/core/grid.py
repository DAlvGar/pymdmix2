"""
3D Grid for Density and Volumetric Data
========================================

Grid class for storing and manipulating 3D volumetric data,
particularly density grids from MD trajectory analysis.

Supports OpenDX format for visualization in VMD, PyMOL, etc.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import numpy as np
from numpy.typing import NDArray


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
    spacing : float
        Grid spacing in Angstroms (uniform in all dimensions)

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
    spacing: float

    def __post_init__(self):
        """Validate grid data."""
        if self.data.ndim != 3:
            raise ValueError(f"Grid data must be 3D, got {self.data.ndim}D")
        if self.spacing <= 0:
            raise ValueError(f"Grid spacing must be positive, got {self.spacing}")

    @property
    def shape(self) -> tuple[int, int, int]:
        """Grid dimensions (nx, ny, nz)."""
        return self.data.shape  # type: ignore

    @property
    def dimensions(self) -> tuple[float, float, float]:
        """Physical dimensions in Angstroms."""
        return tuple(s * self.spacing for s in self.shape)  # type: ignore

    @property
    def center(self) -> tuple[float, float, float]:
        """Grid center coordinates."""
        return tuple(
            self.origin[i] + self.dimensions[i] / 2
            for i in range(3)
        )  # type: ignore

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
        # Convert to grid indices
        ix = int((x - self.origin[0]) / self.spacing)
        iy = int((y - self.origin[1]) / self.spacing)
        iz = int((z - self.origin[2]) / self.spacing)
        
        # Check bounds
        if (0 <= ix < self.shape[0] and
            0 <= iy < self.shape[1] and
            0 <= iz < self.shape[2]):
            return float(self.data[ix, iy, iz])
        
        return None

    @classmethod
    def from_bounds(
        cls,
        min_coord: NDArray[np.float64] | tuple[float, float, float],
        max_coord: NDArray[np.float64] | tuple[float, float, float],
        spacing: float = 0.5,
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

        Returns
        -------
        Grid
            Empty grid covering the bounding box
        """
        min_c = np.asarray(min_coord)
        max_c = np.asarray(max_coord)

        shape = tuple(
            int(np.ceil((max_c[i] - min_c[i]) / spacing)) + 1
            for i in range(3)
        )

        return cls(
            data=np.zeros(shape, dtype=np.float64),
            origin=tuple(min_c),
            spacing=spacing,
        )

    @classmethod
    def from_structure(
        cls,
        coords: NDArray[np.float64],
        spacing: float = 0.5,
        padding: float = 5.0,
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

        Returns
        -------
        Grid
            Empty grid encompassing the structure
        """
        min_coord = coords.min(axis=0) - padding
        max_coord = coords.max(axis=0) + padding
        return cls.from_bounds(min_coord, max_coord, spacing)

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
        idx = ((c - np.array(self.origin)) / self.spacing).astype(int)
        return tuple(idx)  # type: ignore

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
        return tuple(
            self.origin[i] + index[i] * self.spacing
            for i in range(3)
        )  # type: ignore

    def is_inside(self, coord: NDArray[np.float64]) -> bool:
        """Check if coordinate is inside grid bounds."""
        try:
            idx = self.coord_to_index(coord)
            return all(0 <= idx[i] < self.shape[i] for i in range(3))
        except (IndexError, ValueError):
            return False

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
        # Convert all coordinates to indices
        indices = ((coords - np.array(self.origin)) / self.spacing).astype(int)

        # Filter to valid indices
        valid = np.all(
            (indices >= 0) & (indices < np.array(self.shape)),
            axis=1
        )
        valid_indices = indices[valid]

        # Add counts using numpy bincount for efficiency
        if len(valid_indices) > 0:
            flat_indices = np.ravel_multi_index(
                valid_indices.T, self.shape
            )
            counts = np.bincount(flat_indices, minlength=self.data.size)
            self.data += counts.reshape(self.shape)

        return int(valid.sum())

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
        )

    def to_free_energy(self, temperature: float = 300.0) -> Grid:
        """
        Convert density to free energy (ΔG = -RT ln(ρ/ρ0)).

        Parameters
        ----------
        temperature : float
            Temperature in Kelvin

        Returns
        -------
        Grid
            Grid with free energy values in kcal/mol
        """
        R = 0.001987204  # kcal/(mol·K)
        RT = R * temperature

        # Avoid log(0) by setting minimum density
        rho = np.maximum(self.data, 1e-10)
        rho0 = rho.mean()  # Reference density

        dg = -RT * np.log(rho / rho0)

        return Grid(
            data=dg,
            origin=self.origin,
            spacing=self.spacing,
        )

    def write_dx(self, path: str | Path) -> None:
        """
        Write grid in OpenDX format.

        Parameters
        ----------
        path : str or Path
            Output file path
        """
        path = Path(path)
        nx, ny, nz = self.shape
        n_total = nx * ny * nz

        with open(path, 'w') as f:
            # Header
            f.write("# OpenDX density grid\n")
            f.write("# Generated by pyMDMix\n")
            f.write(f"object 1 class gridpositions counts {nx} {ny} {nz}\n")
            f.write(f"origin {self.origin[0]:.6f} {self.origin[1]:.6f} {self.origin[2]:.6f}\n")
            f.write(f"delta {self.spacing:.6f} 0.000000 0.000000\n")
            f.write(f"delta 0.000000 {self.spacing:.6f} 0.000000\n")
            f.write(f"delta 0.000000 0.000000 {self.spacing:.6f}\n")
            f.write(f"object 2 class gridconnections counts {nx} {ny} {nz}\n")
            f.write(f"object 3 class array type double rank 0 items {n_total} data follows\n")

            # Data in Fortran order (z fastest, then y, then x)
            count = 0
            for ix in range(nx):
                for iy in range(ny):
                    for iz in range(nz):
                        f.write(f"{self.data[ix, iy, iz]:.6e}")
                        count += 1
                        if count % 3 == 0:
                            f.write("\n")
                        else:
                            f.write(" ")

            if count % 3 != 0:
                f.write("\n")

            f.write('attribute "dep" string "positions"\n')
            f.write('object "density" class field\n')
            f.write('component "positions" value 1\n')
            f.write('component "connections" value 2\n')
            f.write('component "data" value 3\n')

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
        spacing = None
        shape = None
        data_lines: list[str] = []
        reading_data = False

        with open(path, 'r') as f:
            for line in f:
                line = line.strip()

                if line.startswith('#') or not line:
                    continue

                if 'gridpositions counts' in line:
                    parts = line.split()
                    shape = (int(parts[-3]), int(parts[-2]), int(parts[-1]))

                elif line.startswith('origin'):
                    parts = line.split()
                    origin = (float(parts[1]), float(parts[2]), float(parts[3]))

                elif line.startswith('delta'):
                    parts = line.split()
                    # Take the non-zero delta value
                    deltas = [float(parts[1]), float(parts[2]), float(parts[3])]
                    for d in deltas:
                        if d > 0:
                            spacing = d
                            break

                elif 'data follows' in line:
                    reading_data = True
                    continue

                elif reading_data:
                    if line.startswith('attribute') or line.startswith('object') or line.startswith('component'):
                        break
                    data_lines.append(line)

        if origin is None or spacing is None or shape is None:
            raise ValueError(f"Could not parse DX file: {path}")

        # Parse data
        data_str = ' '.join(data_lines)
        data_flat = np.fromstring(data_str, sep=' ')

        expected_size = shape[0] * shape[1] * shape[2]
        if len(data_flat) != expected_size:
            raise ValueError(
                f"Data size mismatch: expected {expected_size}, got {len(data_flat)}"
            )

        # Reshape (data is in Fortran order: z fastest)
        data = data_flat.reshape(shape, order='C')

        return cls(data=data, origin=origin, spacing=spacing)

    def __add__(self, other: Grid) -> Grid:
        """Add two grids element-wise."""
        if self.shape != other.shape:
            raise ValueError(f"Shape mismatch: {self.shape} vs {other.shape}")
        if self.origin != other.origin:
            raise ValueError(f"Origin mismatch: {self.origin} vs {other.origin}")
        if self.spacing != other.spacing:
            raise ValueError(f"Spacing mismatch: {self.spacing} vs {other.spacing}")

        return Grid(
            data=self.data + other.data,
            origin=self.origin,
            spacing=self.spacing,
        )

    def __mul__(self, scalar: float) -> Grid:
        """Multiply grid by scalar."""
        return Grid(
            data=self.data * scalar,
            origin=self.origin,
            spacing=self.spacing,
        )

    def __repr__(self) -> str:
        return (
            f"Grid(shape={self.shape}, origin={self.origin}, "
            f"spacing={self.spacing}, range=[{self.data.min():.3f}, {self.data.max():.3f}])"
        )
