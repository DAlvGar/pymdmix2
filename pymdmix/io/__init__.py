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

from pymdmix.io.grids import (
    read_dx,
    write_dx,
    read_mrc,
    write_mrc,
    convert_grid,
    GridFormat,
)
from pymdmix.io.plotting import (
    plot_energy_grid,
    plot_convergence,
    plot_rmsd,
)
from pymdmix.io.parsers import (
    SystemConfigFileParser,
    MDSettingsConfigFileParser,
    ParserError,
    SystemParserError,
    BadFile,
    MDSettingsParserError,
    BadSolvent,
)
from pymdmix.io.off_manager import (
    OFFManager,
    OFFManagerError,
    OFFSectionError,
    OFFUnitNotFoundError,
    Atom,
    Residue,
)
from pymdmix.io.dcd_parser import (
    DCDReader,
    DCDFrame,
    DCDHeader,
    DCDError,
    DCDFormatError,
    read_dcd,
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
