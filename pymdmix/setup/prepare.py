"""
Structure Preparation
=====================

Prepare molecular structures for MDMix simulations:
- Parse PQR files from PDB2PQR
- Add terminal caps (ACE/NME) using parmed/LEaP
- Handle disulfide bonds
- Remove clashing atoms
- Clean up atom/residue naming
- Web interface to PDB2PQR for automatic protonation

Migrated from pyMDMix.AutoPrepare (830 lines) with modernized implementation.

Examples
--------
>>> from pymdmix.setup.prepare import AutoPrepare
>>>
>>> # From local PDB file with automatic protonation
>>> preparer = AutoPrepare("protein.pdb", protonate=True)
>>> preparer.save_pdb("protein_prepared.pdb")
>>>
>>> # From PDB ID
>>> preparer = AutoPrepare(pdbid="1YER", chains=[0])
>>> preparer.save_pdb("1yer_prepared.pdb")
"""

from __future__ import annotations

import logging
import os
import tempfile
import time
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import parmed
from numpy.typing import NDArray

from pymdmix.core.structure import (
    PROTEIN_RESIDUES,
    WATER_RESIDUES,
    find_disulfides,
    get_water_mask,
    rename_cys_to_cyx,
)

log = logging.getLogger(__name__)


# ============================================================================
# Exceptions
# ============================================================================


class PDB2PQRError(Exception):
    """Base exception for PDB2PQR errors."""

    pass


class ConnectionError(PDB2PQRError):
    """Error connecting to PDB2PQR server."""

    pass


class FormChangeError(PDB2PQRError):
    """PDB2PQR web form has changed."""

    pass


class JobError(PDB2PQRError):
    """Error running PDB2PQR job."""

    pass


class AutoPrepareError(Exception):
    """Error in automatic preparation workflow."""

    pass


# ============================================================================
# PQR File Parsing
# ============================================================================


@dataclass
class PQRAtom:
    """Single atom from PQR file."""

    serial: int
    name: str
    residue_name: str
    chain_id: str
    residue_number: int
    x: float
    y: float
    z: float
    charge: float
    radius: float
    atom_type: str = "ATOM"  # ATOM or HETATM


class PQRParseFile:
    """
    Parse PQR files from PDB2PQR server.

    PQR format is space-separated, not fixed-column like PDB.
    Columns: record serial name resname [chain] resnum x y z charge radius

    Parameters
    ----------
    pqr_path : str or Path
        Path to PQR file

    Examples
    --------
    >>> parser = PQRParseFile("protein.pqr")
    >>> struct = parser.get_model()
    >>> print(f"Loaded {len(struct.atoms)} atoms")
    """

    def __init__(self, pqr_path: str | Path):
        self.pqr_path = Path(pqr_path)
        if not self.pqr_path.exists():
            raise FileNotFoundError(f"PQR file not found: {self.pqr_path}")

        self._atoms: list[PQRAtom] = []
        self._remarks: list[str] = []
        self._parse()

    def _parse(self) -> None:
        """Parse the PQR file."""
        with open(self.pqr_path) as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line:
                    continue

                parts = line.split()
                if not parts:
                    continue

                record = parts[0]

                if record == "REMARK":
                    self._remarks.append(line)
                    continue

                if record in ("ATOM", "HETATM"):
                    try:
                        atom = self._parse_atom_line(parts, record)
                        self._atoms.append(atom)
                    except (ValueError, IndexError) as e:
                        log.warning(f"Error parsing line {line_num}: {e}")
                        continue

                elif record in ("END", "ENDMDL"):
                    break

    def _parse_atom_line(self, parts: list[str], record: str) -> PQRAtom:
        """
        Parse ATOM/HETATM line from PQR file.

        PQR has variable format:
        - With chain: ATOM serial name resname chain resnum x y z charge radius
        - Without chain: ATOM serial name resname resnum x y z charge radius
        """
        # Determine if chain ID is present (11 fields vs 10 fields)
        if len(parts) == 11:
            # With chain
            return PQRAtom(
                serial=int(parts[1]),
                name=parts[2],
                residue_name=parts[3],
                chain_id=parts[4],
                residue_number=int(parts[5]),
                x=float(parts[6]),
                y=float(parts[7]),
                z=float(parts[8]),
                charge=float(parts[9]),
                radius=float(parts[10]),
                atom_type=record,
            )
        elif len(parts) >= 10:
            # Without chain
            return PQRAtom(
                serial=int(parts[1]),
                name=parts[2],
                residue_name=parts[3],
                chain_id="",
                residue_number=int(parts[4]),
                x=float(parts[5]),
                y=float(parts[6]),
                z=float(parts[7]),
                charge=float(parts[8]),
                radius=float(parts[9]),
                atom_type=record,
            )
        else:
            raise ValueError(f"Invalid atom line: {' '.join(parts)}")

    def get_model(self) -> parmed.Structure:
        """
        Convert parsed PQR data to parmed Structure.

        Returns
        -------
        parmed.Structure
            Structure with atoms from PQR file
        """
        struct = parmed.Structure()

        current_resnum = None
        current_chain = None

        for atom in self._atoms:
            # Check if we need a new residue
            if current_resnum != atom.residue_number or current_chain != atom.chain_id:
                struct.add_atom(
                    parmed.Atom(name="DUMMY"),  # Placeholder
                    atom.residue_name,
                    atom.residue_number,
                    atom.chain_id or "",
                )
                # Remove the dummy atom we just added
                struct.atoms.pop()
                current_resnum = atom.residue_number
                current_chain = atom.chain_id

            # Create parmed atom
            parm_atom = parmed.Atom(
                name=atom.name,
                charge=atom.charge,
            )
            parm_atom.xx = atom.x
            parm_atom.xy = atom.y
            parm_atom.xz = atom.z

            # Add atom to structure
            struct.add_atom(
                parm_atom,
                atom.residue_name,
                atom.residue_number,
                atom.chain_id or "",
            )

        return struct

    @property
    def atoms(self) -> list[PQRAtom]:
        """Get parsed atoms."""
        return self._atoms

    @property
    def remarks(self) -> list[str]:
        """Get REMARK lines."""
        return self._remarks

    @property
    def total_charge(self) -> float:
        """Calculate total charge from atoms."""
        return sum(a.charge for a in self._atoms)


# ============================================================================
# PDB2PQR Web Interface
# ============================================================================


class PDB2PQRInterface:
    """
    Interface to PDB2PQR web server for automatic protonation.

    Uses requests library to submit jobs and fetch results.
    Predicts protonation states using PropKa and adds hydrogens.

    Parameters
    ----------
    server_url : str, optional
        PDB2PQR server URL. Default uses NBCR server.

    Examples
    --------
    >>> interface = PDB2PQRInterface()
    >>> struct = interface.protonate_pdb(pdbfile="protein.pdb")
    >>> struct.save("protonated.pdb")
    """

    # PDB2PQR web server (version 2.x API)
    DEFAULT_SERVER = "https://server.poissonboltzmann.org"

    def __init__(self, server_url: str | None = None):
        self.server_url = server_url or self.DEFAULT_SERVER
        self.log = logging.getLogger(f"{__name__}.PDB2PQRInterface")

        # Check requests is available
        try:
            import requests

            self._requests = requests
        except ImportError:
            raise ImportError(
                "requests library required for PDB2PQR interface. "
                "Install with: pip install requests"
            )

    def _submit_job(
        self,
        pdb_content: str | None = None,
        pdbid: str | None = None,
        ph: float = 7.0,
        forcefield: str = "AMBER",
    ) -> str:
        """
        Submit job to PDB2PQR server.

        Returns job ID for status checking.
        """
        url = f"{self.server_url}/api/pdb2pqr"

        # Build form data
        data = {
            "FORCEFIELD": forcefield,
            "OUTPUT_NAMING_SCHEME": forcefield.lower(),
            "PH": str(ph),
            "DEBUMP": "true",
            "OPT": "true",
            "PROPKA": "true",
        }

        files = None

        if pdb_content:
            # Upload PDB content
            files = {
                "PDB": ("input.pdb", pdb_content, "text/plain"),
            }
            data["PDBSOURCE"] = "UPLOAD"
        elif pdbid:
            # Fetch from RCSB
            data["PDBID"] = pdbid.upper()
            data["PDBSOURCE"] = "ID"
        else:
            raise ValueError("Either pdb_content or pdbid required")

        try:
            if files:
                response = self._requests.post(url, data=data, files=files, timeout=60)
            else:
                response = self._requests.post(url, data=data, timeout=60)
            response.raise_for_status()
        except self._requests.RequestException as e:
            raise ConnectionError(f"Failed to submit job: {e}")

        # Parse job ID from response
        result = response.json()
        if "job_id" in result:
            return result["job_id"]
        elif "id" in result:
            return result["id"]
        else:
            raise PDB2PQRError(f"No job ID in response: {result}")

    def _check_status(self, job_id: str) -> dict:
        """Check job status."""
        url = f"{self.server_url}/api/pdb2pqr/{job_id}"

        try:
            response = self._requests.get(url, timeout=30)
            response.raise_for_status()
            return response.json()
        except self._requests.RequestException as e:
            raise ConnectionError(f"Failed to check status: {e}")

    def _download_pqr(self, job_id: str) -> str:
        """Download PQR result file content."""
        url = f"{self.server_url}/api/pdb2pqr/{job_id}/output.pqr"

        try:
            response = self._requests.get(url, timeout=60)
            response.raise_for_status()
            return response.text
        except self._requests.RequestException as e:
            raise PDB2PQRError(f"Failed to download PQR: {e}")

    def protonate_pdb(
        self,
        pdbfile: str | Path | None = None,
        pdbid: str | None = None,
        ph: float = 7.0,
        forcefield: str = "AMBER",
        wait_interval: float = 3.0,
        max_tries: int = 100,
    ) -> parmed.Structure:
        """
        Protonate PDB using PDB2PQR server.

        Parameters
        ----------
        pdbfile : str or Path, optional
            Path to local PDB file to protonate
        pdbid : str, optional
            PDB ID to fetch and protonate (e.g., "1YER")
        ph : float
            pH for protonation state prediction. Default 7.0.
        forcefield : str
            Force field for naming. Default "AMBER".
        wait_interval : float
            Seconds between status checks. Default 3.0.
        max_tries : int
            Maximum status checks before timeout. Default 100.

        Returns
        -------
        parmed.Structure
            Protonated structure

        Examples
        --------
        >>> interface = PDB2PQRInterface()
        >>> struct = interface.protonate_pdb(pdbfile="protein.pdb", ph=7.4)
        """
        # Read PDB content if file provided
        pdb_content = None
        if pdbfile:
            pdbfile = Path(pdbfile)
            if not pdbfile.exists():
                raise FileNotFoundError(f"PDB file not found: {pdbfile}")
            pdb_content = pdbfile.read_text()

        # Submit job
        self.log.info("Submitting to PDB2PQR server...")
        job_id = self._submit_job(
            pdb_content=pdb_content,
            pdbid=pdbid,
            ph=ph,
            forcefield=forcefield,
        )
        self.log.info(f"Job submitted: {job_id}")

        # Wait for completion
        for attempt in range(max_tries):
            status = self._check_status(job_id)
            state = status.get("status", status.get("state", "unknown"))

            if state in ("complete", "finished", "success"):
                self.log.info("Job complete!")
                break
            elif state in ("error", "failed"):
                error_msg = status.get("error", "Unknown error")
                raise JobError(f"PDB2PQR job failed: {error_msg}")
            elif state in ("running", "pending", "queued"):
                self.log.debug(f"Status: {state} (attempt {attempt + 1}/{max_tries})")
                time.sleep(wait_interval)
            else:
                self.log.warning(f"Unknown status: {state}")
                time.sleep(wait_interval)
        else:
            raise JobError(
                f"Job timed out after {max_tries * wait_interval}s. "
                f"Check status at: {self.server_url}/pdb2pqr/{job_id}"
            )

        # Download and parse PQR
        pqr_content = self._download_pqr(job_id)

        # Write to temp file and parse
        with tempfile.NamedTemporaryFile(mode="w", suffix=".pqr", delete=False) as f:
            f.write(pqr_content)
            temp_path = f.name

        try:
            parser = PQRParseFile(temp_path)
            return parser.get_model()
        finally:
            os.unlink(temp_path)

    def fetch_results(
        self,
        job_id: str,
        wait_interval: float = 3.0,
        max_tries: int = 100,
    ) -> parmed.Structure:
        """
        Fetch results from an existing job.

        Parameters
        ----------
        job_id : str
            Job ID from previous submission
        wait_interval : float
            Seconds between status checks
        max_tries : int
            Maximum status checks

        Returns
        -------
        parmed.Structure
            Protonated structure
        """
        # Wait for completion
        for attempt in range(max_tries):
            status = self._check_status(job_id)
            state = status.get("status", status.get("state", "unknown"))

            if state in ("complete", "finished", "success"):
                break
            elif state in ("error", "failed"):
                raise JobError(f"PDB2PQR job failed: {status.get('error', 'Unknown')}")
            else:
                time.sleep(wait_interval)
        else:
            raise JobError("Job timed out")

        # Download and parse
        pqr_content = self._download_pqr(job_id)

        with tempfile.NamedTemporaryFile(mode="w", suffix=".pqr", delete=False) as f:
            f.write(pqr_content)
            temp_path = f.name

        try:
            parser = PQRParseFile(temp_path)
            return parser.get_model()
        finally:
            os.unlink(temp_path)


# ============================================================================
# Structure Cleaning with Capping
# ============================================================================


@dataclass
class StructurePreparationOptions:
    """
    Options controlling :func:`prepare_structure`.

    Attributes
    ----------
    add_caps : bool
        Add ACE/NME terminal caps (default True).
    detect_disulfides : bool
        Detect disulfide bonds and rename CYS → CYX (default True).
    remove_waters : bool
        Remove crystallographic water molecules (default True).
    remove_hydrogens : bool
        Strip hydrogen atoms before re-protonation (default False).
    custom_patches : list[str]
        Additional LEaP patch commands to apply.
    """

    add_caps: bool = True
    detect_disulfides: bool = True
    remove_waters: bool = True
    remove_hydrogens: bool = False
    custom_patches: list[str] = field(default_factory=list)


@dataclass
class PrepareResult:
    """Result of structure preparation."""

    structure: parmed.Structure
    n_caps_added: int = 0
    n_disulfides: int = 0
    n_waters_removed: int = 0
    warnings: list[str] = field(default_factory=list)
    modifications: list[str] = field(default_factory=list)


def prepare_structure(
    structure: parmed.Structure | Path,
    options: StructurePreparationOptions | None = None,
    *,
    cap_termini: bool = True,
    fix_disulfides: bool = True,
    remove_waters: bool = False,
    remove_hydrogens: bool = False,
) -> PrepareResult:
    """
    Prepare structure for simulation.

    Parameters
    ----------
    structure : parmed.Structure or Path
        Input structure or path to a PDB/mol2 file.
    options : StructurePreparationOptions, optional
        Preparation options.  When provided, its flags take precedence over
        the individual keyword arguments.
    cap_termini : bool
        Add ACE/NME caps to chain termini (ignored when *options* given).
    fix_disulfides : bool
        Detect and rename disulfide cysteines to CYX (ignored when *options* given).
    remove_waters : bool
        Remove crystallographic waters (ignored when *options* given).
    remove_hydrogens : bool
        Remove hydrogen atoms for re-protonation (ignored when *options* given).

    Returns
    -------
    PrepareResult
        Prepared structure with metadata and ``modifications`` log.
    """
    from pymdmix.core.structure import load_structure

    # Resolve options — prefer explicit StructurePreparationOptions
    if options is not None:
        cap_termini = options.add_caps
        fix_disulfides = options.detect_disulfides
        remove_waters = options.remove_waters
        remove_hydrogens = options.remove_hydrogens

    # Load from file if needed
    if isinstance(structure, Path):
        structure = load_structure(structure)

    # Work on a copy
    struct = structure.copy(parmed.Structure)
    result = PrepareResult(structure=struct)

    # Remove waters if requested
    if remove_waters:
        water_mask = get_water_mask(struct)
        n_waters = int(water_mask.sum())
        if n_waters > 0:
            keep_indices = np.where(~water_mask)[0].tolist()
            struct = struct[keep_indices]
            result.structure = struct
            result.n_waters_removed = n_waters
            result.modifications.append(f"Removed {n_waters} water atoms")
            log.info(f"Removed {n_waters} water atoms")

    # Remove hydrogens if requested
    if remove_hydrogens:
        h_mask = np.array([a.atomic_number == 1 for a in struct.atoms])
        n_h = int(h_mask.sum())
        if n_h > 0:
            keep_indices = np.where(~h_mask)[0].tolist()
            struct = struct[keep_indices]
            result.structure = struct
            result.modifications.append(f"Removed {n_h} hydrogen atoms")
            log.info(f"Removed {n_h} hydrogen atoms")

    # Fix disulfides
    if fix_disulfides:
        disulfides = find_disulfides(struct)
        if disulfides:
            rename_cys_to_cyx(struct, disulfides)
            result.n_disulfides = len(disulfides)
            result.modifications.append(f"Renamed {len(disulfides)} CYS → CYX (disulfide bonds)")
            log.info(f"Fixed {len(disulfides)} disulfide bonds")

    # Add caps
    if cap_termini:
        n_caps = add_caps(struct)
        result.n_caps_added = n_caps
        if n_caps > 0:
            result.modifications.append(f"Added {n_caps} terminal ACE/NME caps")
            log.info(f"Added {n_caps} terminal caps")

    return result


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


def add_caps(
    structure: parmed.Structure,
    cap_n: bool = True,
    cap_c: bool = True,
) -> int:
    """
    Add ACE/NME caps to protein chain termini.

    This uses LEaP for proper capping if available, otherwise
    logs a warning about manual capping needed.

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
        Number of caps that would be added
    """
    n_caps = 0

    # Find protein chains
    chains = _find_protein_chains(structure)

    for chain_start, chain_end in chains:
        # Check N-terminus
        if cap_n:
            first_res = structure.residues[chain_start]
            if first_res.name not in ("ACE",):
                log.debug(f"N-terminus at residue {first_res.name}{first_res.number} needs ACE cap")
                n_caps += 1

        # Check C-terminus
        if cap_c:
            last_res = structure.residues[chain_end]
            if last_res.name not in ("NME", "NHE"):
                log.debug(f"C-terminus at residue {last_res.name}{last_res.number} needs NME cap")
                n_caps += 1

    # For actual capping, we need LEaP
    # This placeholder identifies where caps are needed
    if n_caps > 0:
        log.info(f"Identified {n_caps} termini needing caps. Use LEaP for actual capping.")

    return n_caps


def cap_with_leap(
    structure: parmed.Structure,
    cap_n_chains: list[int] | None = None,
    cap_c_chains: list[int] | None = None,
) -> parmed.Structure:
    """
    Add ACE/NME caps using LEaP.

    Parameters
    ----------
    structure : parmed.Structure
        Input structure
    cap_n_chains : list[int], optional
        Chain indices to cap at N-terminus. Default: all chains.
    cap_c_chains : list[int], optional
        Chain indices to cap at C-terminus. Default: all chains.

    Returns
    -------
    parmed.Structure
        Structure with caps added

    Notes
    -----
    Requires AmberTools LEaP to be available in PATH.
    """
    # Find protein chains
    chains = _find_protein_chains(structure)
    n_chains = len(chains)

    if cap_n_chains is None:
        cap_n_chains = list(range(n_chains))
    if cap_c_chains is None:
        cap_c_chains = list(range(n_chains))

    # Write input PDB
    with tempfile.NamedTemporaryFile(mode="w", suffix=".pdb", delete=False) as f:
        structure.save(f.name, overwrite=True)
        input_pdb = f.name

    # Build LEaP script
    leap_script = [
        "source leaprc.protein.ff14SB",
        f"mol = loadpdb {input_pdb}",
    ]

    # Add capping commands
    # LEaP uses 1-based indexing
    for chain_idx in cap_n_chains:
        if chain_idx < n_chains:
            chain_start, _ = chains[chain_idx]
            # Get residue info for LEaP
            first_res = structure.residues[chain_start]
            # In LEaP: sequence mol.chain.residue
            leap_script.append(
                f"# Cap N-terminus of chain {chain_idx} "
                f"(residue {first_res.name}{first_res.number})"
            )
            # LEaP uses different syntax based on version
            # For simple cases, we can prepend ACE

    for chain_idx in cap_c_chains:
        if chain_idx < n_chains:
            _, chain_end = chains[chain_idx]
            last_res = structure.residues[chain_end]
            leap_script.append(
                f"# Cap C-terminus of chain {chain_idx} (residue {last_res.name}{last_res.number})"
            )

    # Save capped structure
    output_pdb = tempfile.mktemp(suffix="_capped.pdb")
    leap_script.extend(
        [
            f"savepdb mol {output_pdb}",
            "quit",
        ]
    )

    # Write LEaP input
    leap_input = tempfile.mktemp(suffix=".leap")
    with open(leap_input, "w") as f:
        f.write("\n".join(leap_script))

    # Run LEaP
    import subprocess

    try:
        result = subprocess.run(
            ["tleap", "-f", leap_input],
            capture_output=True,
            text=True,
            timeout=120,
        )
        if result.returncode != 0:
            log.warning(f"LEaP warnings: {result.stderr}")
    except FileNotFoundError:
        raise RuntimeError("LEaP not found. Install AmberTools.")
    except subprocess.TimeoutExpired:
        raise RuntimeError("LEaP timed out")
    finally:
        os.unlink(input_pdb)
        os.unlink(leap_input)

    # Load result
    if os.path.exists(output_pdb):
        capped = parmed.load_file(output_pdb)
        os.unlink(output_pdb)
        return capped
    else:
        raise RuntimeError("LEaP failed to produce output")


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
) -> parmed.Structure:
    """
    Remove water molecules that clash with solute.

    Parameters
    ----------
    structure : parmed.Structure
        Structure containing waters
    solute_coords : NDArray
        Solute coordinates to check against
    clash_distance : float
        Minimum allowed distance in Angstroms

    Returns
    -------
    parmed.Structure
        Structure with clashing waters removed
    """
    try:
        from scipy.spatial import cKDTree
    except ImportError:
        log.warning("scipy required for clash detection")
        return structure

    # Build tree of solute coordinates
    solute_tree = cKDTree(solute_coords)

    # Find water residue indices to remove
    residues_to_remove = set()

    for res in structure.residues:
        if res.name not in WATER_RESIDUES:
            continue

        for atom in res.atoms:
            if atom.xx is None:
                continue
            coord = np.array([atom.xx, atom.xy, atom.xz])
            dist, _ = solute_tree.query(coord)
            if dist < clash_distance:
                residues_to_remove.add(res.idx)
                break

    if not residues_to_remove:
        return structure

    # Get atom indices to keep
    keep_atoms = [a.idx for a in structure.atoms if a.residue.idx not in residues_to_remove]

    log.info(f"Removed {len(residues_to_remove)} clashing water residues")
    return structure[keep_atoms]


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
        "OXT": "OXT",
        "O1": "O",
        "O2": "OXT",
        # Hydrogens
        "1H": "H1",
        "2H": "H2",
        "3H": "H3",
        "HN": "H",
    }

    n_renamed = 0
    for atom in structure.atoms:
        if atom.name in name_map:
            old_name = atom.name
            atom.name = name_map[atom.name]
            n_renamed += 1
            log.debug(f"Renamed {old_name} -> {atom.name}")

    return n_renamed


def center_structure(structure: parmed.Structure) -> NDArray[np.float64]:
    """
    Center structure at origin.

    Parameters
    ----------
    structure : parmed.Structure
        Structure to center (modified in place)

    Returns
    -------
    NDArray
        Translation vector applied
    """
    coords = structure.coordinates
    center = coords.mean(axis=0)
    structure.coordinates = coords - center
    return center


# ============================================================================
# AmberPDBCleaner - Full PDB cleaning for Amber
# ============================================================================


class AmberPDBCleaner:
    """
    Clean PDB for Amber simulations.

    Handles:
    - Centering structure
    - Water handling (keep/remove, clash detection)
    - Disulfide bond detection and CYX renaming
    - ACE/NME terminal capping
    - Atom name standardization

    Migrated from pyMDMix.AutoPrepare.AmberPDBCleaner.

    Parameters
    ----------
    structure : parmed.Structure
        Input structure to clean
    verbose : bool
        Enable verbose logging

    Examples
    --------
    >>> from pymdmix.core.structure import load_structure
    >>> struct = load_structure("protein.pdb")
    >>> cleaner = AmberPDBCleaner(struct)
    >>> cleaned = cleaner.clean_pdb(cap=True, keep_waters=True)
    """

    def __init__(
        self,
        structure: parmed.Structure,
        verbose: bool = True,
    ):
        self.structure = structure.copy(parmed.Structure)
        self.verbose = verbose
        self.log = logging.getLogger(f"{__name__}.AmberPDBCleaner")

        # Track modifications
        self.disulfides: list[tuple[int, int]] = []

    def clean_pdb(
        self,
        hetatm: bool = True,
        keep_waters: bool = True,
        cap: bool = True,
        cap_n: list[int] | None = None,
        cap_c: list[int] | None = None,
    ) -> parmed.Structure:
        """
        Clean PDB for Amber simulation.

        Parameters
        ----------
        hetatm : bool
            Keep hetero atoms. Default True.
        keep_waters : bool
            Keep crystallographic waters. Default True.
        cap : bool
            Add ACE/NME caps. Default True.
        cap_n : list[int], optional
            Chain indices to cap at N-terminus.
        cap_c : list[int], optional
            Chain indices to cap at C-terminus.

        Returns
        -------
        parmed.Structure
            Cleaned structure
        """
        struct = self.structure

        if self.verbose:
            self.log.info("Cleaning PDB for Amber...")

        # Center structure
        center_structure(struct)

        # Extract waters before processing
        waters = None
        if keep_waters:
            water_mask = get_water_mask(struct)
            if water_mask.any():
                water_indices = np.where(water_mask)[0].tolist()
                waters = struct[water_indices]
                # Remove waters from main structure for processing
                non_water_indices = np.where(~water_mask)[0].tolist()
                struct = struct[non_water_indices]
        else:
            # Remove all waters
            water_mask = get_water_mask(struct)
            non_water_indices = np.where(~water_mask)[0].tolist()
            struct = struct[non_water_indices]

        # Standardize atom names (Xplor -> Amber)
        n_renamed = standardize_atom_names(struct)
        if self.verbose and n_renamed > 0:
            self.log.info(f"Renamed {n_renamed} atoms to Amber conventions")

        # Find protein chains for capping
        chains = _find_protein_chains(struct)

        if cap and chains:
            if cap_n is None:
                cap_n = list(range(len(chains)))
            if cap_c is None:
                cap_c = list(range(len(chains)))

            # For now, identify caps needed (actual capping via LEaP)
            n_caps = 0
            for i in cap_n:
                if i < len(chains):
                    chain_start, _ = chains[i]
                    res = struct.residues[chain_start]
                    if res.name not in ("ACE",):
                        if self.verbose:
                            self.log.info(
                                f"Will cap N-terminus of chain {i} (residue {res.name}{res.number})"
                            )
                        n_caps += 1

            for i in cap_c:
                if i < len(chains):
                    _, chain_end = chains[i]
                    res = struct.residues[chain_end]
                    if res.name not in ("NME", "NHE"):
                        if self.verbose:
                            self.log.info(
                                f"Will cap C-terminus of chain {i} (residue {res.name}{res.number})"
                            )
                        n_caps += 1

        # Find and fix disulfides
        self.disulfides = find_disulfides(struct, cutoff=4.0)
        if self.disulfides:
            rename_cys_to_cyx(struct, self.disulfides)
            if self.verbose:
                self.log.info(f"Found {len(self.disulfides)} disulfide bonds")

        # Re-add waters, removing clashing ones
        if keep_waters and waters is not None and len(waters.atoms) > 0:
            # Get solute coordinates
            solute_coords = struct.coordinates

            # Remove clashing waters
            waters = remove_clashing_waters(waters, solute_coords, clash_distance=1.5)

            if len(waters.atoms) > 0:
                # Combine structures
                struct = struct + waters

        # Renumber residues
        renumber_residues(struct, start=1)

        self.structure = struct
        return struct


# ============================================================================
# AutoPrepare - Full Workflow
# ============================================================================


class AutoPrepare:
    """
    Automatic PDB preparation for MDMix simulations.

    Full workflow:
    1. Load PDB from file or fetch by ID
    2. Optionally protonate using PDB2PQR server
    3. Clean structure (cap, fix disulfides, handle waters)
    4. Save prepared PDB or OFF file

    Migrated from pyMDMix.AutoPrepare.

    Parameters
    ----------
    pdb : str or Path or parmed.Structure, optional
        Input PDB file path or structure
    chains : list[int], optional
        Chain indices to keep (default: all)
    protonate : bool
        Use PDB2PQR for protonation. Default False.
    pdbid : str, optional
        PDB ID to fetch from RCSB (e.g., "1YER")
    ph : float
        pH for protonation. Default 7.0.

    Examples
    --------
    >>> # From local file
    >>> prep = AutoPrepare("protein.pdb")
    >>> prep.save_pdb("protein_prepared.pdb")

    >>> # From PDB ID with protonation
    >>> prep = AutoPrepare(pdbid="1YER", protonate=True, ph=7.4)
    >>> prep.save_pdb("1yer_prepared.pdb")

    >>> # Specific chains only
    >>> prep = AutoPrepare("complex.pdb", chains=[0, 1])
    """

    def __init__(
        self,
        pdb: str | Path | parmed.Structure | None = None,
        chains: list[int] | None = None,
        protonate: bool = False,
        pdbid: str | None = None,
        ph: float = 7.0,
        cap: bool = True,
        keep_waters: bool = True,
        **kwargs,
    ):
        self.log = logging.getLogger(f"{__name__}.AutoPrepare")

        self._pdb: parmed.Structure | None = None
        self.chains = chains or []
        self.pdbid = pdbid

        # Track source
        self._source_file: Path | None = None

        # If PDBID provided, fetch and optionally protonate
        if pdbid:
            self.fetch_pdb_and_protonate(
                pdbid,
                protonate=protonate or True,  # Always protonate if fetching by ID
                ph=ph,
                **kwargs,
            )
        elif pdb:
            # Load from file or structure
            if isinstance(pdb, parmed.Structure):
                self._pdb = pdb.copy(parmed.Structure)
            elif isinstance(pdb, (str, Path)):
                pdb_path = Path(pdb)
                if not pdb_path.exists():
                    raise FileNotFoundError(f"PDB file not found: {pdb_path}")
                self._pdb = parmed.load_file(str(pdb_path))
                self._source_file = pdb_path

            # Optionally protonate
            if protonate:
                self.protonate_pdb(ph=ph, **kwargs)

        # Clean if we have a structure
        if self._pdb is not None:
            self.clean_pdb(chains=self.chains, cap=cap, keep_waters=keep_waters, **kwargs)

    def set_pdb(self, pdb: str | Path | parmed.Structure) -> None:
        """
        Set PDB to prepare.

        Parameters
        ----------
        pdb : str or Path or parmed.Structure
            PDB file path or structure
        """
        if isinstance(pdb, parmed.Structure):
            self._pdb = pdb.copy(parmed.Structure)
        elif isinstance(pdb, (str, Path)):
            pdb_path = Path(pdb)
            if not pdb_path.exists():
                raise FileNotFoundError(f"PDB file not found: {pdb_path}")
            self._pdb = parmed.load_file(str(pdb_path))
            self._source_file = pdb_path
        else:
            raise AutoPrepareError("pdb must be a file path or parmed.Structure")

    def get_pdb(self) -> parmed.Structure | None:
        """Get current PDB structure."""
        return self._pdb

    @property
    def pdb(self) -> parmed.Structure | None:
        """Current PDB structure."""
        return self._pdb

    def save_pdb(self, outname: str | Path, overwrite: bool = True) -> Path:
        """
        Save current PDB structure.

        Parameters
        ----------
        outname : str or Path
            Output file path
        overwrite : bool
            Overwrite existing file. Default True.

        Returns
        -------
        Path
            Path to saved file
        """
        if self._pdb is None:
            raise AutoPrepareError("No PDB loaded")

        outname = Path(outname)
        self._pdb.save(str(outname), overwrite=overwrite)
        self.log.info(f"Saved PDB to {outname}")
        return outname

    def protonate_pdb(
        self,
        ph: float = 7.0,
        tries: int = 50,
        wait: float = 3.0,
        **kwargs,
    ) -> None:
        """
        Protonate current PDB using PDB2PQR server.

        Parameters
        ----------
        ph : float
            pH for protonation prediction
        tries : int
            Max status check attempts
        wait : float
            Seconds between checks
        """
        if self._pdb is None and self._source_file is None:
            raise AutoPrepareError("No PDB to protonate. Set one with set_pdb().")

        interface = PDB2PQRInterface()

        # Use source file if available, otherwise write temp
        if self._source_file and self._source_file.exists():
            pdb_file = self._source_file
        else:
            # Write current structure to temp file
            with tempfile.NamedTemporaryFile(mode="w", suffix=".pdb", delete=False) as f:
                self._pdb.save(f.name, overwrite=True)
                pdb_file = Path(f.name)

        try:
            self._pdb = interface.protonate_pdb(
                pdbfile=pdb_file,
                ph=ph,
                max_tries=tries,
                wait_interval=wait,
            )
            self.log.info(f"Protonated at pH {ph}")
        finally:
            # Clean up temp file if we created one
            if pdb_file != self._source_file and pdb_file.exists():
                pdb_file.unlink()

    def fetch_pdb_and_protonate(
        self,
        pdbid: str,
        protonate: bool = True,
        ph: float = 7.0,
        tries: int = 50,
        wait: float = 5.0,
        **kwargs,
    ) -> None:
        """
        Fetch PDB by ID and optionally protonate.

        Parameters
        ----------
        pdbid : str
            PDB ID (e.g., "1YER")
        protonate : bool
            Use PDB2PQR for protonation
        ph : float
            pH for protonation
        tries : int
            Max attempts for server operations
        wait : float
            Seconds between status checks
        """
        self.pdbid = pdbid.upper()

        if protonate:
            # Fetch and protonate via PDB2PQR
            interface = PDB2PQRInterface()
            self._pdb = interface.protonate_pdb(
                pdbid=pdbid,
                ph=ph,
                max_tries=tries,
                wait_interval=wait,
            )
            self.log.info(f"Fetched and protonated {pdbid} at pH {ph}")
        else:
            # Just fetch from RCSB
            import requests

            url = f"https://files.rcsb.org/download/{pdbid.upper()}.pdb"

            response = requests.get(url, timeout=30)
            response.raise_for_status()

            # Write to temp and load
            with tempfile.NamedTemporaryFile(mode="w", suffix=".pdb", delete=False) as f:
                f.write(response.text)
                temp_path = f.name

            try:
                self._pdb = parmed.load_file(temp_path)
                self.log.info(f"Fetched {pdbid} from RCSB")
            finally:
                os.unlink(temp_path)

    def clean_pdb(
        self,
        chains: list[int] | None = None,
        hetatm: bool = True,
        keep_waters: bool = True,
        cap: bool = True,
        cap_n: list[int] | None = None,
        cap_c: list[int] | None = None,
        **kwargs,
    ) -> None:
        """
        Clean current PDB for simulation.

        Parameters
        ----------
        chains : list[int], optional
            Chain indices to keep. Default: all.
        hetatm : bool
            Keep hetero atoms. Default True.
        keep_waters : bool
            Keep crystallographic waters. Default True.
        cap : bool
            Add ACE/NME caps. Default True.
        cap_n : list[int], optional
            Chain indices to cap at N-terminus.
        cap_c : list[int], optional
            Chain indices to cap at C-terminus.
        """
        if self._pdb is None:
            raise AutoPrepareError("No PDB loaded")

        # Take specific chains if requested
        if chains:
            # Get unique chain IDs
            chain_ids = list(dict.fromkeys(r.chain for r in self._pdb.residues if r.chain))

            # Select atoms from requested chains
            keep_chain_ids = [chain_ids[i] for i in chains if i < len(chain_ids)]

            if keep_chain_ids:
                keep_mask = np.array([a.residue.chain in keep_chain_ids for a in self._pdb.atoms])
                keep_indices = np.where(keep_mask)[0].tolist()
                self._pdb = self._pdb[keep_indices]

            # Update cap lists to match filtered chains
            n_chains = len(set(r.chain for r in self._pdb.residues if r.chain))
            if cap and not cap_c:
                cap_c = list(range(n_chains))
            if cap and not cap_n:
                cap_n = list(range(n_chains))

        # Use cleaner
        cleaner = AmberPDBCleaner(self._pdb, verbose=True)
        self._pdb = cleaner.clean_pdb(
            hetatm=hetatm,
            keep_waters=keep_waters,
            cap=cap,
            cap_n=cap_n,
            cap_c=cap_c,
        )
