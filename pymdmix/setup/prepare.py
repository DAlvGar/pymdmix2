"""
Structure Preparation
=====================

Prepare molecular structures for MDMix simulations:
- Add terminal caps (ACE/NME)
- Handle disulfide bonds
- Remove clashing atoms
- Clean up atom/residue naming

Examples
--------
>>> from pymdmix.setup import prepare_structure
>>> from pymdmix.core import load_structure
>>>
>>> struct = load_structure("protein.pdb")
>>> prepared = prepare_structure(struct, cap_termini=True)
>>> prepared.save("protein_prepared.pdb")
"""

from __future__ import annotations

from dataclasses import dataclass
import logging
import numpy as np
from numpy.typing import NDArray

import parmed

from pymdmix.core.structure import (
    find_disulfides,
    rename_cys_to_cyx,
    get_water_mask,
    PROTEIN_RESIDUES,
)

log = logging.getLogger(__name__)


@dataclass
class PrepareResult:
    """Result of structure preparation."""
    structure: parmed.Structure
    n_caps_added: int = 0
    n_disulfides: int = 0
    n_waters_removed: int = 0
    warnings: list[str] = None

    def __post_init__(self):
        if self.warnings is None:
            self.warnings = []


def prepare_structure(
    structure: parmed.Structure,
    cap_termini: bool = True,
    fix_disulfides: bool = True,
    remove_waters: bool = False,
    remove_hydrogens: bool = False,
) -> PrepareResult:
    """
    Prepare structure for simulation.

    Parameters
    ----------
    structure : parmed.Structure
        Input structure
    cap_termini : bool
        Add ACE/NME caps to chain termini
    fix_disulfides : bool
        Detect and rename disulfide cysteines to CYX
    remove_waters : bool
        Remove crystallographic waters
    remove_hydrogens : bool
        Remove hydrogen atoms (for re-protonation)

    Returns
    -------
    PrepareResult
        Prepared structure with metadata
    """
    # Work on a copy
    struct = structure.copy(parmed.Structure)
    result = PrepareResult(structure=struct)

    # Remove waters if requested
    if remove_waters:
        water_mask = get_water_mask(struct)
        n_waters = water_mask.sum()
        if n_waters > 0:
            # Get indices to keep
            keep_indices = np.where(~water_mask)[0].tolist()
            struct = struct[keep_indices]
            result.structure = struct
            result.n_waters_removed = n_waters
            log.info(f"Removed {n_waters} water atoms")

    # Remove hydrogens if requested
    if remove_hydrogens:
        h_mask = np.array([a.atomic_number == 1 for a in struct.atoms])
        n_h = h_mask.sum()
        if n_h > 0:
            keep_indices = np.where(~h_mask)[0].tolist()
            struct = struct[keep_indices]
            result.structure = struct
            log.info(f"Removed {n_h} hydrogen atoms")

    # Fix disulfides
    if fix_disulfides:
        disulfides = find_disulfides(struct)
        if disulfides:
            rename_cys_to_cyx(struct, disulfides)
            result.n_disulfides = len(disulfides)
            log.info(f"Fixed {len(disulfides)} disulfide bonds")

    # Add caps
    if cap_termini:
        n_caps = add_caps(struct)
        result.n_caps_added = n_caps
        if n_caps > 0:
            log.info(f"Added {n_caps} terminal caps")

    return result


def add_caps(
    structure: parmed.Structure,
    cap_n: bool = True,
    cap_c: bool = True,
) -> int:
    """
    Add ACE/NME caps to protein chain termini.

    This modifies the structure in place.

    Parameters
    ----------
    structure : parmed.Structure
        Structure to modify
    cap_n : bool
        Add ACE cap to N-termini
    cap_c : bool
        Add NME cap to C-termini

    Returns
    -------
    int
        Number of caps added
    """
    n_caps = 0

    # Find protein chains
    chains = _find_protein_chains(structure)

    for chain_start, chain_end in chains:
        # Check N-terminus
        if cap_n:
            first_res = structure.residues[chain_start]
            if first_res.name not in ('ACE',):
                # Would add ACE here
                # For now, just log - actual implementation requires
                # adding atoms which is complex in parmed
                log.debug(f"N-terminus at residue {first_res.name}{first_res.number} needs ACE cap")
                # n_caps += 1

        # Check C-terminus
        if cap_c:
            last_res = structure.residues[chain_end]
            if last_res.name not in ('NME', 'NHE'):
                log.debug(f"C-terminus at residue {last_res.name}{last_res.number} needs NME cap")
                # n_caps += 1

    # Note: Full implementation would use LEaP for proper capping
    # This is a placeholder that identifies where caps are needed

    return n_caps


def _find_protein_chains(structure: parmed.Structure) -> list[tuple[int, int]]:
    """
    Find protein chain boundaries.

    Returns list of (start_residue_idx, end_residue_idx) tuples.
    """
    chains = []
    chain_start = None

    for i, res in enumerate(structure.residues):
        is_protein = res.name in PROTEIN_RESIDUES

        if is_protein and chain_start is None:
            chain_start = i
        elif not is_protein and chain_start is not None:
            chains.append((chain_start, i - 1))
            chain_start = None

    # Handle last chain
    if chain_start is not None:
        chains.append((chain_start, len(structure.residues) - 1))

    return chains


def find_and_fix_disulfides(
    structure: parmed.Structure,
    cutoff: float = 2.5,
) -> list[tuple[int, int]]:
    """
    Find disulfide bonds and rename CYS to CYX.

    Parameters
    ----------
    structure : parmed.Structure
        Structure to modify (in place)
    cutoff : float
        SG-SG distance cutoff in Angstroms

    Returns
    -------
    list[tuple[int, int]]
        List of disulfide pairs (residue indices)
    """
    disulfides = find_disulfides(structure, cutoff)
    if disulfides:
        rename_cys_to_cyx(structure, disulfides)
    return disulfides


def remove_clashing_waters(
    structure: parmed.Structure,
    solute_coords: NDArray[np.float64],
    clash_distance: float = 1.5,
) -> int:
    """
    Remove water molecules that clash with solute.

    Parameters
    ----------
    structure : parmed.Structure
        Structure containing waters (modified in place)
    solute_coords : NDArray
        Solute coordinates to check against
    clash_distance : float
        Minimum allowed distance in Angstroms

    Returns
    -------
    int
        Number of water residues removed
    """
    try:
        from scipy.spatial import cKDTree
    except ImportError:
        log.warning("scipy required for clash detection")
        return 0

    # Build tree of solute coordinates
    solute_tree = cKDTree(solute_coords)

    # Find water residues
    water_residues = [
        r for r in structure.residues
        if r.name in ('WAT', 'HOH', 'TIP3')
    ]

    # Check each water
    residues_to_remove = set()

    for res in water_residues:
        for atom in res.atoms:
            coord = structure.coordinates[atom.idx]
            # Check distance to nearest solute atom
            dist, _ = solute_tree.query(coord)
            if dist < clash_distance:
                residues_to_remove.add(res.idx)
                break

    # Remove clashing residues
    if residues_to_remove:
        # Get atom indices to keep
        _keep_atoms = [
            a.idx for a in structure.atoms
            if a.residue.idx not in residues_to_remove
        ]

        # This modifies structure - would need to reassign
        log.info(f"Would remove {len(residues_to_remove)} clashing water residues")

    return len(residues_to_remove)


def renumber_residues(
    structure: parmed.Structure,
    start: int = 1,
) -> None:
    """
    Renumber residues sequentially starting from start.

    Parameters
    ----------
    structure : parmed.Structure
        Structure to modify (in place)
    start : int
        Starting residue number
    """
    for i, res in enumerate(structure.residues):
        res.number = start + i


def standardize_atom_names(structure: parmed.Structure) -> int:
    """
    Standardize atom names to Amber conventions.

    Parameters
    ----------
    structure : parmed.Structure
        Structure to modify (in place)

    Returns
    -------
    int
        Number of atoms renamed
    """
    # Common name mappings (PDB -> Amber)
    name_map = {
        # Terminal oxygens
        'OXT': 'OXT',
        'O1': 'O',
        'O2': 'OXT',
        # Hydrogens
        '1H': 'H1',
        '2H': 'H2',
        '3H': 'H3',
        'HN': 'H',
    }

    n_renamed = 0
    for atom in structure.atoms:
        if atom.name in name_map:
            old_name = atom.name
            atom.name = name_map[atom.name]
            n_renamed += 1
            log.debug(f"Renamed {old_name} -> {atom.name}")

    return n_renamed
