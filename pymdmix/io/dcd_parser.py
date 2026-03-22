"""
NAMD DCD Trajectory Parser
==========================

Binary DCD file parser for NAMD trajectories.
Handles 32-bit big-endian DCD format with box information.

Examples
--------
>>> from pymdmix.io.dcd_parser import DCDReader
>>> with DCDReader("trajectory.dcd") as reader:
...     print(f"Frames: {reader.n_frames}, Atoms: {reader.n_atoms}")
...     for frame in reader:
...         print(frame.coordinates.shape)
"""

from __future__ import annotations

import logging
import struct
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path
from typing import BinaryIO

import numpy as np
from numpy.typing import NDArray

log = logging.getLogger(__name__)


class DCDError(Exception):
    """Error reading DCD file."""

    pass


class DCDFormatError(DCDError):
    """Invalid DCD format."""

    pass


@dataclass
class DCDHeader:
    """
    DCD file header information.

    Attributes
    ----------
    n_frames : int
        Number of coordinate sets (frames)
    istart : int
        Starting timestep
    nsavc : int
        Steps between DCD saves
    n_total_steps : int
        Total number of simulation steps
    delta : float
        Timestep of simulation
    has_extra_block : bool
        Whether file has CHARMM extra block (unit cell info)
    has_4dims : bool
        Whether file has 4th dimension
    n_atoms : int
        Number of atoms
    title : list[str]
        Title lines from DCD header
    """

    n_frames: int
    istart: int
    nsavc: int
    n_total_steps: int
    delta: float
    has_extra_block: bool
    has_4dims: bool
    n_atoms: int
    title: list[str]

    @property
    def header_size(self) -> int:
        """Size of header in bytes."""
        return 116 + (80 * len(self.title))

    @property
    def frame_size(self) -> int:
        """Size of each frame in bytes."""
        # 4 bytes per float * 3 axes * n_atoms + enclosing integers (6 * 4)
        f_size = (3 * 4 * self.n_atoms) + 24
        if self.has_extra_block:
            f_size += 56
        return f_size


@dataclass
class DCDFrame:
    """
    Single DCD trajectory frame.

    Attributes
    ----------
    coordinates : NDArray[np.float32]
        Atomic coordinates (n_atoms, 3) in Angstroms
    unitcell : NDArray[np.float64] | None
        Unit cell parameters if present
    """

    coordinates: NDArray[np.float32]
    unitcell: NDArray[np.float64] | None = None

    @property
    def n_atoms(self) -> int:
        """Number of atoms."""
        return self.coordinates.shape[0]


class DCDReader:
    """
    NAMD DCD trajectory file reader.

    Reads 32-bit big-endian DCD format with optional box information.
    Supports random frame access and iteration.

    Parameters
    ----------
    path : str or Path
        Path to DCD file

    Examples
    --------
    >>> # Read all frames
    >>> with DCDReader("trajectory.dcd") as reader:
    ...     coords = reader.read_all()
    ...     print(coords.shape)  # (n_frames, n_atoms, 3)

    >>> # Iterate through frames
    >>> with DCDReader("trajectory.dcd") as reader:
    ...     for frame in reader:
    ...         center_of_mass = frame.coordinates.mean(axis=0)

    >>> # Random access
    >>> with DCDReader("trajectory.dcd") as reader:
    ...     frame_100 = reader[100]
    """

    def __init__(self, path: str | Path):
        self.path = Path(path)
        if not self.path.exists():
            raise FileNotFoundError(f"DCD file not found: {self.path}")

        self._file: BinaryIO | None = None
        self._header: DCDHeader | None = None
        self._open()
        self._read_header()

        log.debug(f"Opened DCD: {self.path} ({self.n_frames} frames, {self.n_atoms} atoms)")

    def _open(self) -> None:
        """Open the DCD file."""
        self._file = open(self.path, "rb", buffering=0)

    def _read_header(self) -> None:
        """Read and parse DCD header."""
        f = self._file
        if f.tell() != 0:
            f.seek(0)

        # First header block: 92 bytes
        header_bytes = f.read(92)
        if len(header_bytes) < 92:
            raise DCDFormatError(
                f"DCD file too small: expected at least 92 header bytes, got {len(header_bytes)}"
            )
        header_data = struct.unpack(">I 4s 9I f 11I", header_bytes)

        # Validate format
        # header_data[0] = opening block size (84)
        # header_data[1] = magic (b'CORD')
        # header_data[13] = CHARMM version/4dims flag (non-zero for CHARMM format)
        # header_data[-1] = closing block size (84)
        if not (
            header_data[0] == 84
            and header_data[1] == b"CORD"
            and header_data[-1] == 84
            and header_data[13] != 0
        ):
            raise DCDFormatError("Invalid DCD format. Expected 32-bit big-endian CORD format.")

        n_frames = header_data[2]
        istart = header_data[3]
        nsavc = header_data[4]
        n_total_steps = header_data[5]
        delta = header_data[11]
        has_extra_block = bool(header_data[12])
        has_4dims = bool(header_data[13])

        # Read title block
        title_header = struct.unpack(">I", f.read(4))[0]
        if (title_header - 4) % 80 != 0:
            raise DCDFormatError("Invalid title block in DCD header")

        n_title = struct.unpack(">I", f.read(4))[0]
        title = []
        for _ in range(n_title):
            title_line = struct.unpack(">80s", f.read(80))[0]
            title.append(title_line.decode("ascii", errors="ignore").strip())
        f.read(4)  # Skip closing block

        # Read atom count block
        atom_block = struct.unpack(">3I", f.read(12))
        if atom_block[-1] != 4:
            raise DCDFormatError("Invalid atom count block in DCD header")

        n_atoms = atom_block[1]

        self._header = DCDHeader(
            n_frames=n_frames,
            istart=istart,
            nsavc=nsavc,
            n_total_steps=n_total_steps,
            delta=delta,
            has_extra_block=has_extra_block,
            has_4dims=has_4dims,
            n_atoms=n_atoms,
            title=title,
        )

    @property
    def header(self) -> DCDHeader:
        """DCD header information."""
        return self._header

    @property
    def n_frames(self) -> int:
        """Number of frames."""
        return self._header.n_frames

    @property
    def n_atoms(self) -> int:
        """Number of atoms."""
        return self._header.n_atoms

    def _read_extra_block(self) -> NDArray[np.float64]:
        """Read CHARMM extra block with unit cell information."""
        f = self._file
        block_size = struct.unpack(">I", f.read(4))[0]
        if block_size != 48:
            raise DCDFormatError("Invalid extra block size in DCD frame")

        unitcell = np.frombuffer(f.read(48), dtype=">f8")
        f.read(4)  # Skip closing block
        return unitcell

    def _read_frame(self) -> DCDFrame:
        """Read a single frame at current file position."""
        f = self._file

        # Read unit cell if present
        unitcell = None
        if self._header.has_extra_block:
            unitcell = self._read_extra_block()

        # Read coordinates (x, y, z stored separately)
        n_atoms = self._header.n_atoms
        float_size = struct.calcsize("f")

        coords = np.zeros((n_atoms, 3), dtype=">f4")

        f.read(4)  # Opening block
        coords[:, 0] = np.frombuffer(f.read(float_size * n_atoms), dtype=">f4")
        f.read(8)  # Closing + opening blocks
        coords[:, 1] = np.frombuffer(f.read(float_size * n_atoms), dtype=">f4")
        f.read(8)  # Closing + opening blocks
        coords[:, 2] = np.frombuffer(f.read(float_size * n_atoms), dtype=">f4")
        f.read(4)  # Closing block

        return DCDFrame(coordinates=coords, unitcell=unitcell)

    def read_all(self) -> NDArray[np.float32]:
        """
        Read all frames into a single array.

        Returns
        -------
        NDArray[np.float32]
            Coordinates array with shape (n_frames, n_atoms, 3)
        """
        f = self._file
        f.seek(self._header.header_size)

        all_frames = np.zeros((self._header.n_frames, self._header.n_atoms, 3), dtype=">f4")

        for i in range(self._header.n_frames):
            frame = self._read_frame()
            all_frames[i, :] = frame.coordinates

        return all_frames

    def get_frame(self, index: int) -> DCDFrame:
        """
        Read a specific frame by index.

        Parameters
        ----------
        index : int
            Frame index (0-based)

        Returns
        -------
        DCDFrame
            Frame data
        """
        if index < 0 or index >= self._header.n_frames:
            raise IndexError(f"Frame index {index} out of range [0, {self._header.n_frames})")

        # Seek to frame position
        position = self._header.header_size + (self._header.frame_size * index)
        self._file.seek(position)

        return self._read_frame()

    def __getitem__(self, index: int) -> DCDFrame:
        """Random access to frames."""
        return self.get_frame(index)

    def __iter__(self) -> Iterator[DCDFrame]:
        """Iterate over all frames."""
        self._file.seek(self._header.header_size)
        for _ in range(self._header.n_frames):
            yield self._read_frame()

    def __len__(self) -> int:
        """Number of frames."""
        return self._header.n_frames

    def __enter__(self) -> DCDReader:
        """Context manager entry."""
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        """Context manager exit."""
        self.close()

    def close(self) -> None:
        """Close the DCD file."""
        if self._file is not None:
            self._file.close()
            self._file = None

    def __repr__(self) -> str:
        return f"DCDReader('{self.path.name}', n_frames={self.n_frames}, n_atoms={self.n_atoms})"


def read_dcd(path: str | Path) -> NDArray[np.float32]:
    """
    Convenience function to read all DCD coordinates.

    Parameters
    ----------
    path : str or Path
        Path to DCD file

    Returns
    -------
    NDArray[np.float32]
        Coordinates array with shape (n_frames, n_atoms, 3)
    """
    with DCDReader(path) as reader:
        return reader.read_all()
