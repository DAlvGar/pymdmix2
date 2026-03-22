"""Core data structures and utilities."""

from pymdmix.core.grid import Grid
from pymdmix.core.trajectory import open_trajectory, Frame, TrajectoryReader
from pymdmix.core.structure import load_structure, get_protein_mask, get_water_mask
from pymdmix.core.solvent import Solvent, Probe, SolventLibrary
from pymdmix.core.containers import Atom, Residue
from pymdmix.core.system import (
    System,
    SolvatedSystem,
    SystemError,
    BadFileError,
    FileLock,
    load_system,
)

__all__ = [
    "Grid",
    "open_trajectory",
    "Frame",
    "TrajectoryReader",
    "load_structure",
    "get_protein_mask",
    "get_water_mask",
    "Solvent",
    "Probe",
    "SolventLibrary",
    "Atom",
    "Residue",
    # System management
    "System",
    "SolvatedSystem",
    "SystemError",
    "BadFileError",
    "FileLock",
    "load_system",
]
