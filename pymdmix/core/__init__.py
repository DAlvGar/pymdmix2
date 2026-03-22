"""Core data structures and utilities."""

from pymdmix.core.grid import Grid
from pymdmix.core.trajectory import open_trajectory, Frame, TrajectoryReader
from pymdmix.core.structure import (
    # I/O
    load_structure,
    save_structure,
    # Masks
    get_protein_mask,
    get_water_mask,
    get_nucleic_mask,
    get_ion_mask,
    get_solute_mask,
    get_solvent_mask,
    get_heavy_atom_mask,
    get_hydrogen_mask,
    get_backbone_mask,
    get_residue_mask,
    get_atom_mask,
    # Selection
    select_atoms,
    select_chains,
    # Chain operations
    get_chain_ids,
    get_chain_start_indices,
    # Residue info
    count_residues,
    get_residue_names,
    # Disulfides
    find_disulfides,
    rename_cys_to_cyx,
    # Centering
    center_structure,
    # Profile conversions
    atom_to_residue_values,
    residue_to_atom_values,
    # Probe/density helpers
    get_probe_coords,
    iter_residue_coords,
    get_residue_com_coords,
    detect_solvent_type,
    # Alignment
    align_structures,
    align_to_reference,
    compute_rmsd,
)
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
    # Grid
    "Grid",
    # Trajectory
    "open_trajectory",
    "Frame",
    "TrajectoryReader",
    # Structure - I/O
    "load_structure",
    "save_structure",
    # Structure - Masks
    "get_protein_mask",
    "get_water_mask",
    "get_nucleic_mask",
    "get_ion_mask",
    "get_solute_mask",
    "get_solvent_mask",
    "get_heavy_atom_mask",
    "get_hydrogen_mask",
    "get_backbone_mask",
    "get_residue_mask",
    "get_atom_mask",
    # Structure - Selection
    "select_atoms",
    "select_chains",
    # Structure - Chains
    "get_chain_ids",
    "get_chain_start_indices",
    # Structure - Residue info
    "count_residues",
    "get_residue_names",
    # Structure - Disulfides
    "find_disulfides",
    "rename_cys_to_cyx",
    # Structure - Centering
    "center_structure",
    # Structure - Profile conversions
    "atom_to_residue_values",
    "residue_to_atom_values",
    # Structure - Probe/density
    "get_probe_coords",
    "iter_residue_coords",
    "get_residue_com_coords",
    "detect_solvent_type",
    # Structure - Alignment
    "align_structures",
    "align_to_reference",
    "compute_rmsd",
    # Solvent
    "Solvent",
    "Probe",
    "SolventLibrary",
    # Containers
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
