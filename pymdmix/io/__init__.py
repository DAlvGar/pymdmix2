"""
pyMDMix I/O - File format readers and writers.

Supports:
- DX (OpenDX) format for density grids
- MRC/CCP4 format for electron density
- PDB format for structures
- JSON/YAML for configuration
- System and MD settings configuration files
- Amber OFF/lib files
"""

from __future__ import annotations

from pymdmix.io.dcd_parser import (
    DCDError,
    DCDFormatError,
    DCDFrame,
    DCDHeader,
    DCDReader,
    read_dcd,
)
from pymdmix.io.grids import (
    GridFormat,
    convert_grid,
    read_dx,
    read_mrc,
    write_dx,
    write_mrc,
)
from pymdmix.io.off_manager import (
    Atom,
    OFFManager,
    OFFManagerError,
    OFFSectionError,
    OFFUnitNotFoundError,
    Residue,
)
from pymdmix.io.parsers import (
    BadFile,
    BadSolvent,
    MDSettingsConfigFileParser,
    MDSettingsParserError,
    ParserError,
    SystemConfigFileParser,
    SystemParserError,
)
from pymdmix.io.plotting import (
    plot_convergence,
    plot_energy_grid,
    plot_rmsd,
)

__all__ = [
    # Grid I/O
    "read_dx",
    "write_dx",
    "read_mrc",
    "write_mrc",
    "convert_grid",
    "GridFormat",
    # Plotting
    "plot_energy_grid",
    "plot_convergence",
    "plot_rmsd",
    # Config parsers
    "SystemConfigFileParser",
    "MDSettingsConfigFileParser",
    "ParserError",
    "SystemParserError",
    "BadFile",
    "MDSettingsParserError",
    "BadSolvent",
    # OFF file handling
    "OFFManager",
    "OFFManagerError",
    "OFFSectionError",
    "OFFUnitNotFoundError",
    "Atom",
    "Residue",
    # DCD parsing
    "DCDReader",
    "DCDFrame",
    "DCDHeader",
    "DCDError",
    "DCDFormatError",
    "read_dcd",
]
