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


def get_hydrogen_mask(struct: parmed.Structure) -> NDArray[np.bool_]:
    """
    Get boolean mask for hydrogen atoms.

    Equivalent to Biskit's maskH().

    Parameters
    ----------
    struct : parmed.Structure
        Input structure

    Returns
    -------
    NDArray[np.bool_]
        Boolean mask (True for hydrogen atoms)
    """
    return np.array([
        a.atomic_number == 1
        for a in struct.atoms
    ])


# Backbone atom names for proteins
BACKBONE_ATOMS = frozenset({'N', 'CA', 'C', 'O'})


def get_backbone_mask(struct: parmed.Structure) -> NDArray[np.bool_]:
    """
    Get boolean mask for protein backbone atoms (N, CA, C, O).

    Equivalent to Biskit's maskBB().

    Parameters
    ----------
    struct : parmed.Structure
        Input structure

    Returns
    -------
    NDArray[np.bool_]
        Boolean mask (True for backbone atoms)
    """
    return np.array([
        a.residue.name in PROTEIN_RESIDUES and a.name in BACKBONE_ATOMS
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


def get_chain_start_indices(struct: parmed.Structure) -> list[int]:
    """
    Get atom indices where each chain starts.

    Equivalent to Biskit's chainIndex().

    Parameters
    ----------
    struct : parmed.Structure
        Input structure

    Returns
    -------
    list[int]
        Atom indices marking start of each chain
    """
    if not struct.atoms:
        return []

    indices = [0]
    current_chain = struct.atoms[0].residue.chain

    for i, atom in enumerate(struct.atoms):
        if atom.residue.chain != current_chain:
            indices.append(i)
            current_chain = atom.residue.chain

    return indices


def select_chains(
    struct: parmed.Structure,
    chain_ids: str | Sequence[str],
) -> parmed.Structure:
    """
    Select atoms belonging to specified chains.

    Equivalent to Biskit's takeChains().

    Parameters
    ----------
    struct : parmed.Structure
        Input structure
    chain_ids : str or sequence of str
        Chain ID(s) to select

    Returns
    -------
    parmed.Structure
        New structure with only specified chains
    """
    if isinstance(chain_ids, str):
        chain_set = {chain_ids}
    else:
        chain_set = set(chain_ids)

    mask = np.array([
        a.residue.chain in chain_set
        for a in struct.atoms
    ])
    return select_atoms(struct, mask)


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


# =============================================================================
# Profile Conversions (Biskit atom2resProfile / res2atomProfile equivalents)
# =============================================================================

def atom_to_residue_values(
    struct: parmed.Structure,
    atom_values: Sequence,
) -> list:
    """
    Convert per-atom values to per-residue (takes first atom of each residue).

    Equivalent to Biskit's atom2resProfile.

    Parameters
    ----------
    struct : parmed.Structure
        Input structure
    atom_values : Sequence
        Per-atom values (length must match number of atoms)

    Returns
    -------
    list
        Per-residue values (one per residue)

    Examples
    --------
    >>> bfactors = [a.bfactor for a in struct.atoms]
    >>> res_bfactors = atom_to_residue_values(struct, bfactors)
    """
    if len(atom_values) != len(struct.atoms):
        raise ValueError(
            f"Length mismatch: {len(atom_values)} values for {len(struct.atoms)} atoms"
        )

    result = []
    for res in struct.residues:
        first_atom_idx = res.atoms[0].idx
        result.append(atom_values[first_atom_idx])
    return result


def residue_to_atom_values(
    struct: parmed.Structure,
    res_values: Sequence,
) -> NDArray:
    """
    Expand per-residue values to per-atom array.

    Equivalent to Biskit's res2atomProfile.

    Parameters
    ----------
    struct : parmed.Structure
        Input structure
    res_values : Sequence
        Per-residue values (length must match number of residues)

    Returns
    -------
    NDArray
        Per-atom values

    Examples
    --------
    >>> res_numbers = list(range(len(struct.residues)))
    >>> atom_numbers = residue_to_atom_values(struct, res_numbers)
    """
    if len(res_values) != len(struct.residues):
        raise ValueError(
            f"Length mismatch: {len(res_values)} values for {len(struct.residues)} residues"
        )

    result = np.zeros(len(struct.atoms))
    for i, res in enumerate(struct.residues):
        for atom in res.atoms:
            result[atom.idx] = res_values[i]
    return result


# =============================================================================
# Probe/Density Helpers (for trajectory analysis)
# =============================================================================

def get_probe_coords(
    struct: parmed.Structure,
    residue_name: str,
    atom_names: str | Sequence[str],
) -> NDArray[np.float64]:
    """
    Extract coordinates of probe atoms for density analysis.

    Parameters
    ----------
    struct : parmed.Structure
        Input structure (typically a trajectory frame)
    residue_name : str
        Residue name to select (e.g., 'ETA' for ethanol)
    atom_names : str or sequence of str
        Atom name(s) within the residue (e.g., 'O' for hydroxyl oxygen)

    Returns
    -------
    NDArray[np.float64]
        Coordinates of matching atoms, shape (N, 3)

    Examples
    --------
    >>> # Get ethanol oxygen positions
    >>> coords = get_probe_coords(frame, 'ETA', 'O')
    """
    mask = get_residue_mask(struct, residue_name) & get_atom_mask(struct, atom_names)
    return struct.coordinates[mask]


def iter_residue_coords(
    struct: parmed.Structure,
    residue_name: str,
):
    """
    Iterate over coordinates of each residue with given name.

    Parameters
    ----------
    struct : parmed.Structure
        Input structure
    residue_name : str
        Residue name to iterate over

    Yields
    ------
    NDArray[np.float64]
        Coordinates of atoms in each residue, shape (n_atoms_in_res, 3)

    Examples
    --------
    >>> for res_coords in iter_residue_coords(frame, 'ETA'):
    ...     com = res_coords.mean(axis=0)
    """
    coords = struct.coordinates
    for res in struct.residues:
        if res.name == residue_name:
            indices = [a.idx for a in res.atoms]
            yield coords[indices]


def get_residue_com_coords(
    struct: parmed.Structure,
    residue_name: str,
    mass_weighted: bool = False,
) -> NDArray[np.float64]:
    """
    Get center-of-mass coordinates for each residue of given type.

    Parameters
    ----------
    struct : parmed.Structure
        Input structure
    residue_name : str
        Residue name to process
    mass_weighted : bool
        If True, compute mass-weighted COM. If False, geometric center.

    Returns
    -------
    NDArray[np.float64]
        COM coordinates, shape (n_residues, 3)

    Examples
    --------
    >>> # Get ethanol COMs for density calculation
    >>> coms = get_residue_com_coords(frame, 'ETA')
    """
    coms = []
    coords = struct.coordinates

    for res in struct.residues:
        if res.name == residue_name:
            indices = [a.idx for a in res.atoms]
            res_coords = coords[indices]

            if mass_weighted:
                masses = np.array([a.mass for a in res.atoms])
                com = np.average(res_coords, axis=0, weights=masses)
            else:
                com = res_coords.mean(axis=0)

            coms.append(com)

    return np.array(coms) if coms else np.empty((0, 3))


def detect_solvent_type(
    struct: parmed.Structure,
    known_solvents: dict[str, Sequence[str]],
) -> str | None:
    """
    Auto-detect solvent mixture from residue composition.

    Parameters
    ----------
    struct : parmed.Structure
        Input structure
    known_solvents : dict[str, Sequence[str]]
        Mapping of solvent names to their residue names.
        E.g., {'MAM': ['MAM', 'WAT'], 'ETA': ['ETA', 'WAT']}

    Returns
    -------
    str | None
        Detected solvent name, or None if no match

    Examples
    --------
    >>> solvents = {'MAM': ['MAM'], 'ETA': ['ETA']}
    >>> solvent_type = detect_solvent_type(frame, solvents)
    """
    # Get organic solvent residues (not water, not ions, not solute)
    solute = get_solute_mask(struct)
    water = get_water_mask(struct)
    ions = get_ion_mask(struct)

    organic_mask = ~(solute | water | ions)
    organic_resnames = set(
        a.residue.name for i, a in enumerate(struct.atoms)
        if organic_mask[i]
    )

    # Match against known solvents
    for solvent_name, residue_names in known_solvents.items():
        # Exclude WAT from comparison (it's always present)
        solvent_organics = set(residue_names) - WATER_RESIDUES
        if organic_resnames == solvent_organics:
            return solvent_name

    return None


# =============================================================================
# Structure Alignment (Biskit magicFit equivalent)
# =============================================================================

def _kabsch_rotation(P: NDArray, Q: NDArray) -> NDArray:
    """
    Compute optimal rotation matrix using Kabsch algorithm.

    Parameters
    ----------
    P : NDArray
        Mobile coordinates (N, 3), centered
    Q : NDArray
        Reference coordinates (N, 3), centered

    Returns
    -------
    NDArray
        Rotation matrix (3, 3)
    """
    # Covariance matrix
    H = P.T @ Q

    # SVD
    U, S, Vt = np.linalg.svd(H)

    # Correct for reflection
    d = np.linalg.det(Vt.T @ U.T)
    correction = np.diag([1, 1, d])

    # Optimal rotation
    R = Vt.T @ correction @ U.T

    return R


def compute_rmsd(
    coords1: NDArray[np.float64],
    coords2: NDArray[np.float64],
) -> float:
    """
    Compute RMSD between two coordinate sets.

    Parameters
    ----------
    coords1 : NDArray[np.float64]
        First coordinate set (N, 3)
    coords2 : NDArray[np.float64]
        Second coordinate set (N, 3)

    Returns
    -------
    float
        RMSD in same units as input (typically Angstroms)
    """
    if coords1.shape != coords2.shape:
        raise ValueError(
            f"Shape mismatch: {coords1.shape} vs {coords2.shape}"
        )

    diff = coords1 - coords2
    return float(np.sqrt(np.mean(np.sum(diff ** 2, axis=1))))


def align_structures(
    mobile: parmed.Structure,
    reference: parmed.Structure,
    mask: NDArray[np.bool_] | None = None,
    in_place: bool = True,
) -> tuple[parmed.Structure, float]:
    """
    Superimpose mobile structure onto reference using Kabsch algorithm.

    Equivalent to Biskit's magicFit().

    Parameters
    ----------
    mobile : parmed.Structure
        Structure to align
    reference : parmed.Structure
        Reference structure to align to
    mask : NDArray[np.bool_] | None
        Boolean mask for atoms to use in fitting.
        If None, uses all atoms (structures must have same atom count).
        If provided, applies to both structures.
    in_place : bool
        If True, modify mobile coordinates in place.
        If False, return a copy.

    Returns
    -------
    tuple[parmed.Structure, float]
        (Aligned structure, RMSD of fitted atoms)

    Examples
    --------
    >>> # Align using backbone atoms
    >>> bb_mask = get_backbone_mask(mobile)
    >>> aligned, rmsd = align_structures(mobile, reference, mask=bb_mask)
    >>> print(f"RMSD: {rmsd:.2f} Å")

    Notes
    -----
    Uses Kabsch algorithm (SVD-based) for optimal rotation.
    For large structures or trajectories, consider using MDAnalysis directly.
    """
    # Get coordinates
    mobile_coords = mobile.coordinates.copy()
    ref_coords = reference.coordinates.copy()

    if mask is not None:
        # Use masked atoms for fitting
        fit_mobile = mobile_coords[mask]
        fit_ref = ref_coords[mask]
    else:
        if len(mobile_coords) != len(ref_coords):
            raise ValueError(
                f"Atom count mismatch: {len(mobile_coords)} vs {len(ref_coords)}. "
                "Provide a mask to specify which atoms to align."
            )
        fit_mobile = mobile_coords
        fit_ref = ref_coords

    # Center both coordinate sets
    mobile_center = fit_mobile.mean(axis=0)
    ref_center = fit_ref.mean(axis=0)

    fit_mobile_centered = fit_mobile - mobile_center
    fit_ref_centered = fit_ref - ref_center

    # Compute optimal rotation
    rotation = _kabsch_rotation(fit_mobile_centered, fit_ref_centered)

    # Apply transformation to ALL atoms:
    # 1. Translate to origin (using fit atoms center)
    # 2. Rotate
    # 3. Translate to reference center
    all_centered = mobile_coords - mobile_center
    all_rotated = all_centered @ rotation.T
    all_transformed = all_rotated + ref_center

    # Compute RMSD on fitted atoms
    fitted_mobile = (fit_mobile - mobile_center) @ rotation.T + ref_center
    rmsd = compute_rmsd(fitted_mobile, fit_ref)

    # Apply to structure
    if in_place:
        result = mobile
    else:
        result = mobile.copy(parmed.Structure)

    result.coordinates = all_transformed

    log.debug(f"Aligned structures with RMSD: {rmsd:.3f} Å")
    return result, rmsd


def align_to_reference(
    mobile: parmed.Structure,
    reference: parmed.Structure,
    selection: str = "backbone",
) -> tuple[parmed.Structure, float]:
    """
    Convenience function to align using common selections.

    Parameters
    ----------
    mobile : parmed.Structure
        Structure to align
    reference : parmed.Structure
        Reference structure
    selection : str
        Selection type: 'backbone', 'ca', 'heavy', 'all', or 'protein'

    Returns
    -------
    tuple[parmed.Structure, float]
        (Aligned structure, RMSD)
    """
    if selection == "backbone":
        mask = get_backbone_mask(mobile)
    elif selection == "ca":
        mask = get_protein_mask(mobile) & get_atom_mask(mobile, "CA")
    elif selection == "heavy":
        mask = get_heavy_atom_mask(mobile)
    elif selection == "protein":
        mask = get_protein_mask(mobile)
    elif selection == "all":
        mask = None
    else:
        raise ValueError(f"Unknown selection: {selection}")

    return align_structures(mobile, reference, mask=mask)
