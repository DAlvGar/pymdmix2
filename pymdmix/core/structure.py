"""
Structure Handling
==================

Utilities for working with molecular structures using parmed.
Provides convenience functions while keeping parmed.Structure as the core type.

Examples
--------
>>> struct = load_structure("protein.pdb")
>>> protein_mask = get_protein_mask(struct)
>>> protein_coords = struct.coordinates[protein_mask]
"""

from __future__ import annotations

from pathlib import Path
from typing import Sequence
import logging
import numpy as np
from numpy.typing import NDArray

import parmed

log = logging.getLogger(__name__)

# Standard residue sets for masking
PROTEIN_RESIDUES = frozenset({
    # Standard amino acids
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY',
    'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER',
    'THR', 'TRP', 'TYR', 'VAL',
    # Amber protonation states
    'CYX', 'CYM',  # Cysteine variants
    'HID', 'HIE', 'HIP',  # Histidine variants
    'ASH', 'GLH',  # Protonated Asp/Glu
    'LYN',  # Neutral lysine
    # Caps
    'ACE', 'NME', 'NHE',
})

WATER_RESIDUES = frozenset({
    'WAT', 'HOH', 'TIP3', 'TIP4', 'TIP5', 'SPC', 'SPCE',
    'T3P', 'T4P', 'T4E', 'T5P', 'OPC', 'OPC3',
})

DNA_RESIDUES = frozenset({
    'DA', 'DT', 'DG', 'DC', 'DU',
    'DA5', 'DA3', 'DAN',
    'DT5', 'DT3', 'DTN',
    'DG5', 'DG3', 'DGN',
    'DC5', 'DC3', 'DCN',
})

RNA_RESIDUES = frozenset({
    'A', 'U', 'G', 'C',
    'RA', 'RU', 'RG', 'RC',
    'A5', 'A3', 'AN',
    'U5', 'U3', 'UN',
    'G5', 'G3', 'GN',
    'C5', 'C3', 'CN',
})

ION_RESIDUES = frozenset({
    'Na+', 'NA', 'SOD',
    'K+', 'K', 'POT',
    'Cl-', 'CL', 'CLA',
    'Mg2+', 'MG', 'MG2',
    'Ca2+', 'CA', 'CAL',
    'Zn2+', 'ZN',
    'Fe2+', 'FE', 'FE2',
})


def load_structure(path: str | Path) -> parmed.Structure:
    """
    Load structure from file.

    Parameters
    ----------
    path : str or Path
        Path to structure file (PDB, mol2, prmtop, etc.)

    Returns
    -------
    parmed.Structure
        Loaded structure

    Examples
    --------
    >>> struct = load_structure("protein.pdb")
    >>> print(f"Loaded {len(struct.atoms)} atoms")
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Structure file not found: {path}")

    struct = parmed.load_file(str(path))
    log.debug(f"Loaded structure from {path}: {len(struct.atoms)} atoms")

    return struct


def save_structure(
    struct: parmed.Structure,
    path: str | Path,
    overwrite: bool = False,
) -> None:
    """
    Save structure to file.

    Parameters
    ----------
    struct : parmed.Structure
        Structure to save
    path : str or Path
        Output file path
    overwrite : bool
        Overwrite existing file
    """
    path = Path(path)
    if path.exists() and not overwrite:
        raise FileExistsError(f"File exists: {path}. Use overwrite=True.")

    struct.save(str(path), overwrite=overwrite)
    log.debug(f"Saved structure to {path}")


def get_protein_mask(struct: parmed.Structure) -> NDArray[np.bool_]:
    """
    Get boolean mask for protein atoms.

    Parameters
    ----------
    struct : parmed.Structure
        Input structure

    Returns
    -------
    NDArray[np.bool_]
        Boolean mask (True for protein atoms)
    """
    return np.array([
        a.residue.name in PROTEIN_RESIDUES
        for a in struct.atoms
    ])


def get_water_mask(struct: parmed.Structure) -> NDArray[np.bool_]:
    """
    Get boolean mask for water atoms.

    Parameters
    ----------
    struct : parmed.Structure
        Input structure

    Returns
    -------
    NDArray[np.bool_]
        Boolean mask (True for water atoms)
    """
    return np.array([
        a.residue.name in WATER_RESIDUES
        for a in struct.atoms
    ])


def get_nucleic_mask(struct: parmed.Structure) -> NDArray[np.bool_]:
    """
    Get boolean mask for nucleic acid atoms (DNA + RNA).

    Parameters
    ----------
    struct : parmed.Structure
        Input structure

    Returns
    -------
    NDArray[np.bool_]
        Boolean mask (True for nucleic acid atoms)
    """
    nucleic = DNA_RESIDUES | RNA_RESIDUES
    return np.array([
        a.residue.name in nucleic
        for a in struct.atoms
    ])


def get_ion_mask(struct: parmed.Structure) -> NDArray[np.bool_]:
    """
    Get boolean mask for ion atoms.

    Parameters
    ----------
    struct : parmed.Structure
        Input structure

    Returns
    -------
    NDArray[np.bool_]
        Boolean mask (True for ion atoms)
    """
    return np.array([
        a.residue.name in ION_RESIDUES
        for a in struct.atoms
    ])


def get_solute_mask(struct: parmed.Structure) -> NDArray[np.bool_]:
    """
    Get boolean mask for solute atoms (protein + nucleic acids).

    Parameters
    ----------
    struct : parmed.Structure
        Input structure

    Returns
    -------
    NDArray[np.bool_]
        Boolean mask (True for solute atoms)
    """
    return get_protein_mask(struct) | get_nucleic_mask(struct)


def get_solvent_mask(struct: parmed.Structure) -> NDArray[np.bool_]:
    """
    Get boolean mask for solvent atoms (water + ions).

    Parameters
    ----------
    struct : parmed.Structure
        Input structure

    Returns
    -------
    NDArray[np.bool_]
        Boolean mask (True for solvent atoms)
    """
    return get_water_mask(struct) | get_ion_mask(struct)


def get_residue_mask(
    struct: parmed.Structure,
    resnames: str | Sequence[str],
) -> NDArray[np.bool_]:
    """
    Get boolean mask for specific residue names.

    Parameters
    ----------
    struct : parmed.Structure
        Input structure
    resnames : str or sequence of str
        Residue name(s) to match

    Returns
    -------
    NDArray[np.bool_]
        Boolean mask (True for matching atoms)
    """
    if isinstance(resnames, str):
        resnames = {resnames}
    else:
        resnames = set(resnames)

    return np.array([
        a.residue.name in resnames
        for a in struct.atoms
    ])


def get_atom_mask(
    struct: parmed.Structure,
    atom_names: str | Sequence[str],
) -> NDArray[np.bool_]:
    """
    Get boolean mask for specific atom names.

    Parameters
    ----------
    struct : parmed.Structure
        Input structure
    atom_names : str or sequence of str
        Atom name(s) to match

    Returns
    -------
    NDArray[np.bool_]
        Boolean mask (True for matching atoms)
    """
    if isinstance(atom_names, str):
        atom_names = {atom_names}
    else:
        atom_names = set(atom_names)

    return np.array([
        a.name in atom_names
        for a in struct.atoms
    ])


def get_heavy_atom_mask(struct: parmed.Structure) -> NDArray[np.bool_]:
    """
    Get boolean mask for heavy (non-hydrogen) atoms.

    Parameters
    ----------
    struct : parmed.Structure
        Input structure

    Returns
    -------
    NDArray[np.bool_]
        Boolean mask (True for heavy atoms)
    """
    return np.array([
        a.atomic_number > 1
        for a in struct.atoms
    ])


def select_atoms(
    struct: parmed.Structure,
    mask: NDArray[np.bool_],
) -> parmed.Structure:
    """
    Create new structure with selected atoms.

    Parameters
    ----------
    struct : parmed.Structure
        Input structure
    mask : NDArray[np.bool_]
        Boolean mask for atom selection

    Returns
    -------
    parmed.Structure
        New structure with selected atoms
    """
    indices = np.where(mask)[0].tolist()
    return struct[indices]


def center_structure(struct: parmed.Structure) -> None:
    """
    Center structure at origin (in-place).

    Parameters
    ----------
    struct : parmed.Structure
        Structure to center (modified in place)
    """
    coords = struct.coordinates
    struct.coordinates = coords - coords.mean(axis=0)


def find_disulfides(
    struct: parmed.Structure,
    cutoff: float = 2.5,
) -> list[tuple[int, int]]:
    """
    Find disulfide bonds based on SG-SG distance.

    Parameters
    ----------
    struct : parmed.Structure
        Input structure
    cutoff : float
        Maximum SG-SG distance for disulfide bond (Angstroms)

    Returns
    -------
    list[tuple[int, int]]
        List of (residue_idx1, residue_idx2) pairs
    """
    # Find CYS/CYX SG atoms
    sg_atoms = [
        (i, a) for i, a in enumerate(struct.atoms)
        if a.residue.name in ('CYS', 'CYX') and a.name == 'SG'
    ]

    if len(sg_atoms) < 2:
        return []

    disulfides = []
    coords = struct.coordinates

    for i, (idx1, a1) in enumerate(sg_atoms):
        for idx2, a2 in sg_atoms[i+1:]:
            dist = np.linalg.norm(coords[idx1] - coords[idx2])
            if dist < cutoff:
                disulfides.append((a1.residue.idx, a2.residue.idx))

    log.debug(f"Found {len(disulfides)} disulfide bonds")
    return disulfides


def rename_cys_to_cyx(
    struct: parmed.Structure,
    disulfides: list[tuple[int, int]] | None = None,
) -> None:
    """
    Rename CYS residues involved in disulfides to CYX (in-place).

    Parameters
    ----------
    struct : parmed.Structure
        Structure to modify
    disulfides : list[tuple[int, int]] | None
        Disulfide pairs. If None, auto-detect.
    """
    if disulfides is None:
        disulfides = find_disulfides(struct)

    # Collect residue indices involved in disulfides
    ss_residues = set()
    for res1, res2 in disulfides:
        ss_residues.add(res1)
        ss_residues.add(res2)

    # Rename
    for residue in struct.residues:
        if residue.idx in ss_residues and residue.name == 'CYS':
            residue.name = 'CYX'
            log.debug(f"Renamed residue {residue.idx} from CYS to CYX")


def get_residue_names(struct: parmed.Structure) -> list[str]:
    """Get list of unique residue names in structure."""
    return list(set(r.name for r in struct.residues))


def get_chain_ids(struct: parmed.Structure) -> list[str]:
    """Get list of unique chain IDs in structure."""
    return list(set(r.chain for r in struct.residues if r.chain))


def count_residues(struct: parmed.Structure) -> dict[str, int]:
    """
    Count residues by name.

    Parameters
    ----------
    struct : parmed.Structure
        Input structure

    Returns
    -------
    dict[str, int]
        Residue name -> count mapping
    """
    from collections import Counter
    return dict(Counter(r.name for r in struct.residues))
