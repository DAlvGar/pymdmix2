"""
Trajectory Reading
==================

Abstract trajectory interface with multiple backends:
- MDAnalysis (universal: Amber, Gromacs, NAMD, etc.)
- Native Amber NetCDF (minimal dependencies)

Examples
--------
>>> # Auto-detect best backend
>>> traj = open_trajectory("system.prmtop", "md.nc")
>>> for frame in traj:
...     print(frame.coordinates.shape)

>>> # Use MDAnalysis for Gromacs trajectory
>>> traj = open_trajectory("system.gro", "md.xtc", backend="mdanalysis")
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, Protocol, runtime_checkable
import logging
import numpy as np
from numpy.typing import NDArray

log = logging.getLogger(__name__)

# Check for optional dependencies
try:
    import MDAnalysis as mda
    HAS_MDANALYSIS = True
except ImportError:
    HAS_MDANALYSIS = False
    mda = None

try:
    import netCDF4
    HAS_NETCDF4 = True
except ImportError:
    HAS_NETCDF4 = False
    netCDF4 = None


@dataclass
class Frame:
    """
    Single trajectory frame.

    Attributes
    ----------
    coordinates : NDArray[np.float64]
        Atomic coordinates (n_atoms, 3) in Angstroms
    time : float | None
        Simulation time in picoseconds
    box : NDArray[np.float64] | None
        Box dimensions [a, b, c, alpha, beta, gamma]
    """
    coordinates: NDArray[np.float64]
    time: float | None = None
    box: NDArray[np.float64] | None = None

    @property
    def n_atoms(self) -> int:
        """Number of atoms in frame."""
        return self.coordinates.shape[0]


@runtime_checkable
class TrajectoryReader(Protocol):
    """Protocol defining trajectory reader interface."""

    @property
    def n_frames(self) -> int:
        """Total number of frames."""
        ...

    @property
    def n_atoms(self) -> int:
        """Number of atoms per frame."""
        ...

    def __iter__(self) -> Iterator[Frame]:
        """Iterate over frames."""
        ...

    def __len__(self) -> int:
        """Number of frames."""
        ...


class BaseTrajectoryReader(ABC):
    """Abstract base class for trajectory readers."""

    def __init__(self, topology: str | Path, trajectory: str | Path):
        self.topology_path = Path(topology)
        self.trajectory_path = Path(trajectory)

        if not self.topology_path.exists():
            raise FileNotFoundError(f"Topology not found: {self.topology_path}")
        if not self.trajectory_path.exists():
            raise FileNotFoundError(f"Trajectory not found: {self.trajectory_path}")

    @property
    @abstractmethod
    def n_frames(self) -> int:
        """Total number of frames."""
        ...

    @property
    @abstractmethod
    def n_atoms(self) -> int:
        """Number of atoms per frame."""
        ...

    @abstractmethod
    def __iter__(self) -> Iterator[Frame]:
        """Iterate over frames."""
        ...

    def __len__(self) -> int:
        """Number of frames."""
        return self.n_frames

    def __repr__(self) -> str:
        return (
            f"{self.__class__.__name__}("
            f"trajectory='{self.trajectory_path.name}', "
            f"n_frames={self.n_frames}, n_atoms={self.n_atoms})"
        )


class MDAnalysisReader(BaseTrajectoryReader):
    """
    Universal trajectory reader using MDAnalysis.

    Supports: NetCDF (.nc), DCD, XTC, TRR, and many more.

    Parameters
    ----------
    topology : str or Path
        Topology file (prmtop, pdb, gro, psf, etc.)
    trajectory : str or Path
        Trajectory file (nc, dcd, xtc, trr, etc.)

    Examples
    --------
    >>> reader = MDAnalysisReader("system.prmtop", "md.nc")
    >>> for frame in reader:
    ...     print(frame.coordinates.mean(axis=0))
    """

    def __init__(self, topology: str | Path, trajectory: str | Path):
        if not HAS_MDANALYSIS:
            raise ImportError(
                "MDAnalysis is required for this reader. "
                "Install with: pip install MDAnalysis"
            )

        super().__init__(topology, trajectory)
        self._universe = mda.Universe(
            str(self.topology_path),
            str(self.trajectory_path)
        )
        log.debug(f"Opened trajectory with MDAnalysis: {self.trajectory_path}")

    @property
    def n_frames(self) -> int:
        return len(self._universe.trajectory)

    @property
    def n_atoms(self) -> int:
        return len(self._universe.atoms)

    def __iter__(self) -> Iterator[Frame]:
        for ts in self._universe.trajectory:
            yield Frame(
                coordinates=ts.positions.copy(),
                time=ts.time if hasattr(ts, 'time') else None,
                box=ts.dimensions.copy() if ts.dimensions is not None else None,
            )

    def select_atoms(self, selection: str) -> NDArray[np.int64]:
        """
        Get atom indices matching MDAnalysis selection.

        Parameters
        ----------
        selection : str
            MDAnalysis selection string

        Returns
        -------
        NDArray[np.int64]
            Atom indices matching selection

        Examples
        --------
        >>> indices = reader.select_atoms("resname WAT and name O")
        """
        return self._universe.select_atoms(selection).indices

    @property
    def universe(self):
        """Access underlying MDAnalysis Universe."""
        return self._universe


class AmberNetCDFReader(BaseTrajectoryReader):
    """
    Direct Amber NetCDF trajectory reader.

    Minimal dependency fallback when MDAnalysis is not available.
    Only supports Amber NetCDF format (.nc).

    Parameters
    ----------
    topology : str or Path
        Amber topology file (.prmtop)
    trajectory : str or Path
        Amber NetCDF trajectory (.nc)
    """

    def __init__(self, topology: str | Path, trajectory: str | Path):
        if not HAS_NETCDF4:
            raise ImportError(
                "netCDF4 is required for this reader. "
                "Install with: pip install netCDF4"
            )

        super().__init__(topology, trajectory)

        self._nc = netCDF4.Dataset(str(self.trajectory_path), 'r')
        self._n_frames = self._nc.dimensions['frame'].size
        self._n_atoms = self._nc.dimensions['atom'].size

        log.debug(f"Opened Amber NetCDF trajectory: {self.trajectory_path}")

    @property
    def n_frames(self) -> int:
        return self._n_frames

    @property
    def n_atoms(self) -> int:
        return self._n_atoms

    def __iter__(self) -> Iterator[Frame]:
        coords_var = self._nc.variables['coordinates']

        # Check for optional variables
        time_var = self._nc.variables.get('time')
        box_var = self._nc.variables.get('cell_lengths')
        angles_var = self._nc.variables.get('cell_angles')

        for i in range(self._n_frames):
            coords = coords_var[i, :, :].astype(np.float64)

            time = float(time_var[i]) if time_var is not None else None

            box = None
            if box_var is not None and angles_var is not None:
                lengths = box_var[i, :]
                angles = angles_var[i, :]
                box = np.concatenate([lengths, angles]).astype(np.float64)

            yield Frame(coordinates=coords, time=time, box=box)

    def __del__(self):
        """Close NetCDF file on cleanup."""
        if hasattr(self, '_nc') and self._nc is not None:
            try:
                self._nc.close()
            except Exception:
                pass


class FrameSliceReader:
    """
    Wrapper for reading a slice of frames from a trajectory.

    Parameters
    ----------
    reader : TrajectoryReader
        Base trajectory reader
    start : int
        Starting frame index
    stop : int | None
        Stopping frame index (exclusive)
    step : int
        Frame step
    """

    def __init__(
        self,
        reader: TrajectoryReader,
        start: int = 0,
        stop: int | None = None,
        step: int = 1,
    ):
        self._reader = reader
        self._start = start
        self._stop = stop if stop is not None else reader.n_frames
        self._step = step

    @property
    def n_frames(self) -> int:
        return len(range(self._start, self._stop, self._step))

    @property
    def n_atoms(self) -> int:
        return self._reader.n_atoms

    def __len__(self) -> int:
        return self.n_frames

    def __iter__(self) -> Iterator[Frame]:
        for i, frame in enumerate(self._reader):
            if i < self._start:
                continue
            if i >= self._stop:
                break
            if (i - self._start) % self._step == 0:
                yield frame


def open_trajectory(
    topology: str | Path,
    trajectory: str | Path,
    backend: str = "auto",
) -> TrajectoryReader:
    """
    Open trajectory with appropriate reader.

    Parameters
    ----------
    topology : str or Path
        Topology file (prmtop, pdb, gro, psf, etc.)
    trajectory : str or Path
        Trajectory file (nc, dcd, xtc, trr, etc.)
    backend : str
        Backend to use: 'auto', 'mdanalysis', or 'amber'

    Returns
    -------
    TrajectoryReader
        Trajectory reader instance

    Examples
    --------
    >>> # Auto-detect backend
    >>> traj = open_trajectory("system.prmtop", "md.nc")

    >>> # Force MDAnalysis for Gromacs
    >>> traj = open_trajectory("system.gro", "md.xtc", backend="mdanalysis")

    >>> # Use minimal Amber reader
    >>> traj = open_trajectory("system.prmtop", "md.nc", backend="amber")
    """
    traj_path = Path(trajectory)
    traj_suffix = traj_path.suffix.lower()

    if backend == "auto":
        # Prefer MDAnalysis if available
        if HAS_MDANALYSIS:
            backend = "mdanalysis"
        elif traj_suffix in ('.nc', '.ncdf') and HAS_NETCDF4:
            backend = "amber"
        else:
            raise ImportError(
                f"No suitable backend available for {traj_suffix} format. "
                "Install MDAnalysis: pip install MDAnalysis"
            )

    if backend == "mdanalysis":
        return MDAnalysisReader(topology, trajectory)
    elif backend == "amber":
        if traj_suffix not in ('.nc', '.ncdf'):
            raise ValueError(
                f"Amber backend only supports NetCDF (.nc), got {traj_suffix}"
            )
        return AmberNetCDFReader(topology, trajectory)
    else:
        raise ValueError(f"Unknown backend: {backend}")


def get_available_backends() -> list[str]:
    """Get list of available trajectory backends."""
    backends = []
    if HAS_MDANALYSIS:
        backends.append("mdanalysis")
    if HAS_NETCDF4:
        backends.append("amber")
    return backends
