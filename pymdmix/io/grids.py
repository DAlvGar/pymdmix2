"""
Grid file format readers and writers.

Supports:
- DX (OpenDX) - standard format for density grids
- MRC/CCP4 - electron density format (common in cryo-EM/crystallography)

The DX format is the primary format used by pyMDMix.
MRC support allows interoperability with visualization tools like ChimeraX.
"""

from __future__ import annotations

import struct
from dataclasses import dataclass
from enum import Enum
from pathlib import Path

import numpy as np

from pymdmix.core.grid import Grid


class GridFormat(Enum):
    """Supported grid file formats."""

    DX = "dx"
    MRC = "mrc"
    CCP4 = "ccp4"

    @classmethod
    def from_path(cls, path: str | Path) -> GridFormat:
        """Detect format from file extension."""
        suffix = Path(path).suffix.lower().lstrip(".")
        if suffix in ("dx", "opendx"):
            return cls.DX
        elif suffix in ("mrc", "map", "ccp4"):
            return cls.MRC
        else:
            raise ValueError(f"Unknown grid format: {suffix}")


# =============================================================================
# DX Format (OpenDX)
# =============================================================================


def read_dx(path: str | Path) -> Grid:
    """Read a DX format grid file.

    Args:
        path: Path to DX file

    Returns:
        Grid object

    Note:
        This wraps Grid.read_dx() for convenience.
    """
    return Grid.read_dx(Path(path))


def write_dx(grid: Grid, path: str | Path, precision: int = 6) -> None:
    """Write a grid to DX format.

    Args:
        grid: Grid to write
        path: Output path
        precision: Decimal precision for values (ignored, uses Grid method)

    Note:
        This wraps Grid.write_dx() for convenience.
    """
    grid.write_dx(Path(path))


# =============================================================================
# MRC/CCP4 Format
# =============================================================================


@dataclass
class MRCHeader:
    """MRC file header information."""

    nx: int
    ny: int
    nz: int
    mode: int  # Data type (2 = float32)
    nxstart: int
    nystart: int
    nzstart: int
    mx: int
    my: int
    mz: int
    xlen: float
    ylen: float
    zlen: float
    alpha: float
    beta: float
    gamma: float
    mapc: int  # Column axis (1=X, 2=Y, 3=Z)
    mapr: int  # Row axis
    maps: int  # Section axis
    dmin: float
    dmax: float
    dmean: float
    ispg: int  # Space group
    nsymbt: int  # Extended header size
    origin: tuple[float, float, float]


def read_mrc(path: str | Path) -> Grid:
    """Read an MRC/CCP4 format grid file.

    Args:
        path: Path to MRC/CCP4 file

    Returns:
        Grid object

    Raises:
        ValueError: If file format is invalid

    Note:
        Supports MRC2014 format with machine stamp detection.
        pyMDMix grids use uniform spacing, so the MRC spacing is averaged
        if anisotropic.
    """
    path = Path(path)

    with open(path, "rb") as f:
        # Read header (1024 bytes)
        header_data = f.read(1024)

        if len(header_data) < 1024:
            raise ValueError("Invalid MRC file: header too short")

        # Detect endianness from machine stamp (bytes 212-215)
        machine_stamp = header_data[212:216]

        if machine_stamp[0:2] == b"\x44\x44":
            endian = "<"  # Little endian
        elif machine_stamp[0:2] == b"\x11\x11":
            endian = ">"  # Big endian
        else:
            # Try to detect from NX value
            nx_le = struct.unpack("<i", header_data[0:4])[0]
            endian = "<" if 0 < nx_le < 10000 else ">"

        # Parse header
        nx, ny, nz = struct.unpack(endian + "iii", header_data[0:12])
        mode = struct.unpack(endian + "i", header_data[12:16])[0]

        if mode != 2:
            raise ValueError(f"Only mode 2 (float32) supported, got mode {mode}")

        # Grid start indices
        nxstart, nystart, nzstart = struct.unpack(endian + "iii", header_data[16:28])

        # Unit cell dimensions
        mx, my, mz = struct.unpack(endian + "iii", header_data[28:40])
        xlen, ylen, zlen = struct.unpack(endian + "fff", header_data[40:52])

        # Extended header size
        nsymbt = struct.unpack(endian + "i", header_data[92:96])[0]

        # Origin (MRC2014 style, bytes 196-207)
        origin = struct.unpack(endian + "fff", header_data[196:208])

        # Skip extended header
        if nsymbt > 0:
            f.read(nsymbt)

        # Read data
        data_size = nx * ny * nz * 4  # float32 = 4 bytes
        data_bytes = f.read(data_size)

        if len(data_bytes) < data_size:
            raise ValueError(f"Incomplete data: expected {data_size}, got {len(data_bytes)}")

        data = np.frombuffer(data_bytes, dtype=np.float32)
        if endian == ">":
            data = data.byteswap()

        # Reshape (MRC is column-major: X fastest, Z slowest)
        data = data.reshape((nz, ny, nx))

        # Transpose to (x, y, z) order
        data = data.transpose((2, 1, 0)).copy().astype(np.float64)

        # Calculate grid parameters
        # Spacing - use average if anisotropic
        dx = xlen / nx if nx > 0 else 1.0
        dy = ylen / ny if ny > 0 else 1.0
        dz = zlen / nz if nz > 0 else 1.0
        spacing = (dx + dy + dz) / 3.0

        # Origin
        ox = origin[0] + nxstart * dx
        oy = origin[1] + nystart * dy
        oz = origin[2] + nzstart * dz

        # Create Grid using dataclass directly
        grid = Grid(
            data=data,
            origin=(ox, oy, oz),
            spacing=spacing,
        )

        return grid


def write_mrc(grid: Grid, path: str | Path, label: str = "pyMDMix") -> None:
    """Write a grid to MRC/CCP4 format.

    Args:
        grid: Grid to write
        path: Output path
        label: Label string (max 80 chars)

    Note:
        Writes MRC2014 format with float32 data.
        Uses uniform spacing from Grid.
    """
    path = Path(path)

    nx, ny, nz = grid.shape
    ox, oy, oz = grid.origin
    spacing = grid.spacing

    # Transpose data from (x, y, z) to (z, y, x) for MRC format
    data = grid.data.transpose((2, 1, 0)).astype(np.float32)

    with open(path, "wb") as f:
        # Write header (1024 bytes)
        header = bytearray(1024)

        # NX, NY, NZ (columns, rows, sections)
        struct.pack_into("<iii", header, 0, nx, ny, nz)

        # Mode (2 = float32)
        struct.pack_into("<i", header, 12, 2)

        # Start indices (0)
        struct.pack_into("<iii", header, 16, 0, 0, 0)

        # Grid sampling (same as dimensions)
        struct.pack_into("<iii", header, 28, nx, ny, nz)

        # Cell dimensions in Angstroms (uniform spacing)
        xlen = nx * spacing
        ylen = ny * spacing
        zlen = nz * spacing
        struct.pack_into("<fff", header, 40, xlen, ylen, zlen)

        # Cell angles (90, 90, 90 for orthogonal)
        struct.pack_into("<fff", header, 52, 90.0, 90.0, 90.0)

        # Axis mapping (1=X, 2=Y, 3=Z)
        struct.pack_into("<iii", header, 64, 1, 2, 3)

        # Statistics
        dmin = float(data.min())
        dmax = float(data.max())
        dmean = float(data.mean())
        struct.pack_into("<fff", header, 76, dmin, dmax, dmean)

        # Space group (1 for P1)
        struct.pack_into("<i", header, 88, 1)

        # Extended header size (0)
        struct.pack_into("<i", header, 92, 0)

        # MRC2014 marker
        header[104:108] = b"MAP "

        # Machine stamp (little-endian)
        header[212:216] = b"\x44\x44\x00\x00"

        # RMS deviation
        rms = float(data.std())
        struct.pack_into("<f", header, 216, rms)

        # Number of labels
        struct.pack_into("<i", header, 220, 1)

        # Label (80 chars)
        label_bytes = label[:80].ljust(80).encode("ascii", errors="replace")
        header[224:304] = label_bytes

        # Origin (MRC2014 style)
        struct.pack_into("<fff", header, 196, ox, oy, oz)

        f.write(header)

        # Write data
        f.write(data.tobytes())


# =============================================================================
# Format Conversion
# =============================================================================


def convert_grid(
    input_path: str | Path,
    output_path: str | Path,
    output_format: GridFormat | None = None,
) -> None:
    """Convert a grid between formats.

    Args:
        input_path: Input grid file
        output_path: Output grid file
        output_format: Target format (auto-detected from extension if None)

    Examples:
        >>> convert_grid("density.dx", "density.mrc")
        >>> convert_grid("density.mrc", "density.dx", GridFormat.DX)
    """
    input_path = Path(input_path)
    output_path = Path(output_path)

    # Detect formats
    input_format = GridFormat.from_path(input_path)
    if output_format is None:
        output_format = GridFormat.from_path(output_path)

    # Read input
    if input_format == GridFormat.DX:
        grid = read_dx(input_path)
    elif input_format in (GridFormat.MRC, GridFormat.CCP4):
        grid = read_mrc(input_path)
    else:
        raise ValueError(f"Unsupported input format: {input_format}")

    # Write output
    if output_format == GridFormat.DX:
        write_dx(grid, output_path)
    elif output_format in (GridFormat.MRC, GridFormat.CCP4):
        write_mrc(grid, output_path)
    else:
        raise ValueError(f"Unsupported output format: {output_format}")


# =============================================================================
# Utilities
# =============================================================================


def grid_info(path: str | Path) -> dict:
    """Get information about a grid file without loading all data.

    Args:
        path: Grid file path

    Returns:
        Dictionary with grid metadata
    """
    path = Path(path)
    fmt = GridFormat.from_path(path)

    if fmt == GridFormat.DX:
        # Parse DX header
        with open(path) as f:
            info = {"format": "DX"}
            for line in f:
                if line.startswith("#"):
                    continue
                if "object 1" in line and "gridpositions" in line:
                    parts = line.split()
                    info["shape"] = (int(parts[-3]), int(parts[-2]), int(parts[-1]))
                elif "origin" in line:
                    parts = line.split()
                    info["origin"] = (float(parts[1]), float(parts[2]), float(parts[3]))
                elif "delta" in line:
                    parts = line.split()
                    # First non-zero is the spacing for this dimension
                    deltas = [float(parts[i]) for i in range(1, 4)]
                    if "spacing" not in info:
                        info["spacing"] = []
                    spacing = max(abs(d) for d in deltas)
                    if spacing > 0:
                        info["spacing"].append(spacing)
                elif "object 3" in line:
                    break
            if "spacing" in info and len(info["spacing"]) > 0:
                info["spacing"] = info["spacing"][0]  # Uniform spacing
            return info

    elif fmt in (GridFormat.MRC, GridFormat.CCP4):
        with open(path, "rb") as f:
            header = f.read(1024)
            nx, ny, nz = struct.unpack("<iii", header[0:12])
            xlen, ylen, zlen = struct.unpack("<fff", header[40:52])
            origin = struct.unpack("<fff", header[196:208])
            dmin, dmax, dmean = struct.unpack("<fff", header[76:88])

            return {
                "format": "MRC",
                "shape": (nx, ny, nz),
                "origin": origin,
                "cell_dimensions": (xlen, ylen, zlen),
                "spacing": (xlen / nx + ylen / ny + zlen / nz) / 3.0,
                "min": dmin,
                "max": dmax,
                "mean": dmean,
            }

    else:
        raise ValueError(f"Unsupported format: {fmt}")
