"""
Structure Setup
===============

Tools for preparing molecular structures for MDMix simulations.

- prepare: Structure cleaning, capping, protonation
- solvate: Solvation with organic solvent mixtures
- topology: LEaP integration for topology generation

Classes
-------
- PQRParseFile: Parse PQR files from PDB2PQR
- PDB2PQRInterface: Web interface to PDB2PQR server
- AmberPDBCleaner: Clean PDB for Amber simulations
- AutoPrepare: Full preparation workflow
"""

from pymdmix.setup.prepare import (
    # PDB cleaning
    AmberPDBCleaner,
    # Main workflow
    AutoPrepare,
    AutoPrepareError,
    ConnectionError,
    FormChangeError,
    JobError,
    PDB2PQRError,
    # PDB2PQR interface
    PDB2PQRInterface,
    PQRAtom,
    # PQR parsing
    PQRParseFile,
    PrepareResult,
    add_caps,
    cap_with_leap,
    center_structure,
    find_and_fix_disulfides,
    # Structure preparation
    prepare_structure,
    remove_clashing_waters,
    renumber_residues,
    standardize_atom_names,
)

__all__ = [
    # Main workflow
    "AutoPrepare",
    "AutoPrepareError",
    # Structure preparation
    "prepare_structure",
    "PrepareResult",
    "add_caps",
    "cap_with_leap",
    "find_and_fix_disulfides",
    "remove_clashing_waters",
    "renumber_residues",
    "standardize_atom_names",
    "center_structure",
    # PDB cleaning
    "AmberPDBCleaner",
    # PQR parsing
    "PQRParseFile",
    "PQRAtom",
    # PDB2PQR interface
    "PDB2PQRInterface",
    "PDB2PQRError",
    "ConnectionError",
    "FormChangeError",
    "JobError",
]
