"""
Amber OFF (Object File Format) Manager
=======================================

Provides reading and parsing capabilities for Amber OFF/lib files.
These files contain parameterized units (molecules) with atomic information,
connectivity, coordinates, and box dimensions.

Examples
--------
>>> from pymdmix.io.off_manager import OFFManager
>>> off = OFFManager.from_file("solvent.off")
>>> print(off.get_units())
['ETA', 'ETAWAT20', 'WAT']
>>> coords = off.get_coords("ETA")
>>> volume = off.get_volume("ETAWAT20")
"""

from __future__ import annotations

import logging
import tempfile
from collections.abc import Iterator
from dataclasses import dataclass, field
from functools import reduce
from pathlib import Path

import numpy as np

log = logging.getLogger(__name__)


class OFFManagerError(Exception):
    """Base exception for OFFManager errors."""

    pass


class OFFSectionError(OFFManagerError):
    """Raised when a requested section is not found in the OFF file."""

    pass


class OFFUnitNotFoundError(OFFManagerError):
    """Raised when a requested unit is not found in the OFF file."""

    pass


@dataclass
class Atom:
    """
    Container for atomic information from OFF file.

    Attributes
    ----------
    id : int
        Atom ID (1-indexed in OFF file)
    name : str
        Atom name
    atom_type : str
        Amber atom type
    element : int
        Element atomic number
    charge : float
        Partial charge
    """

    id: int
    name: str
    atom_type: str
    element: int
    charge: float

    def __repr__(self) -> str:
        return f"Atom({self.name!r}, type={self.atom_type!r}, charge={self.charge:.4f})"

    def __eq__(self, other: object) -> bool:
        """Compare by name."""
        if isinstance(other, Atom):
            return self.name == other.name
        if isinstance(other, str):
            return self.name == other
        return NotImplemented


@dataclass
class Residue:
    """
    Container for residue unit information from OFF file.

    Attributes
    ----------
    name : str
        Residue name
    atoms : list[Atom]
        List of atoms in the residue
    connectivity : list[tuple[int, int]]
        Bond connectivity as atom index pairs
    xyz : np.ndarray
        Coordinates (N x 3)
    """

    name: str
    atoms: list[Atom] = field(default_factory=list)
    connectivity: list[tuple[int, int]] = field(default_factory=list)
    xyz: np.ndarray | None = None

    def __post_init__(self):
        """Calculate derived properties."""
        # Build lookup maps
        self._atom_ids: dict[int, Atom] = {a.id: a for a in self.atoms}
        self._atom_names: dict[str, Atom] = {a.name: a for a in self.atoms}

    @property
    def charge(self) -> float:
        """Total charge of the residue."""
        if not self.atoms:
            return 0.0
        return round(sum(a.charge for a in self.atoms), 4)

    @property
    def num_atoms(self) -> int:
        """Number of atoms in the residue."""
        return len(self.atoms)

    def get_atom_by_id(self, atom_id: int) -> Atom | None:
        """Get atom by ID."""
        return self._atom_ids.get(atom_id)

    def get_atom_by_name(self, name: str) -> Atom | None:
        """Get atom by name."""
        return self._atom_names.get(name)

    def center(self) -> np.ndarray | None:
        """Get center of mass of the residue."""
        if self.xyz is None or len(self.xyz) == 0:
            return None
        return self.xyz.mean(axis=0)

    def set_xyz(self, xyz: np.ndarray) -> bool:
        """
        Set new coordinates.

        Parameters
        ----------
        xyz : np.ndarray
            New coordinates array

        Returns
        -------
        bool
            True if successful, False if shape mismatch
        """
        if self.xyz is not None and xyz.shape != self.xyz.shape:
            return False
        self.xyz = xyz
        return True

    def __repr__(self) -> str:
        return f"Residue({self.name!r}, atoms={len(self.atoms)})"

    def __eq__(self, other: object) -> bool:
        """Compare by name."""
        if isinstance(other, Residue):
            return self.name == other.name
        if isinstance(other, str):
            return self.name == other
        return NotImplemented


class OFFManager:
    """
    Manager for Amber OFF (Object File Format) files.

    Provides reading and parsing capabilities for OFF/lib files used by
    tLEaP and other Amber tools.

    Parameters
    ----------
    content : str
        Content of the OFF file as a string

    Attributes
    ----------
    SECTIONS : list[str]
        Known section names in OFF format

    Examples
    --------
    >>> off = OFFManager.from_file("solvent.off")
    >>> units = off.get_units()
    >>> residue = off.get_residue("ETA")
    >>> print(residue.charge)
    """

    SECTIONS = [
        "atoms",
        "atomspertinfo",
        "boundbox",
        "childsequence",
        "connect",
        "connectivity",
        "hierarchy",
        "name",
        "positions",
        "residueconnect",
        "residues",
        "residuesPdbSequenceNumber",
        "solventcap",
        "velocities",
    ]

    def __init__(self, content: str):
        """
        Initialize OFFManager with file content.

        Parameters
        ----------
        content : str
            Content of the OFF file as a string
        """
        self._content = content
        self._tmpfile: Path | None = None

    @classmethod
    def from_file(cls, path: str | Path) -> OFFManager:
        """
        Create OFFManager from a file path.

        Parameters
        ----------
        path : str | Path
            Path to the OFF file

        Returns
        -------
        OFFManager
            Initialized manager

        Raises
        ------
        OFFManagerError
            If file does not exist
        """
        path = Path(path)
        if not path.exists():
            raise OFFManagerError(f"Object File {path} not found")
        content = path.read_text()
        return cls(content)

    @classmethod
    def from_string(cls, content: str) -> OFFManager:
        """
        Create OFFManager from a string.

        Parameters
        ----------
        content : str
            OFF file content as string

        Returns
        -------
        OFFManager
            Initialized manager
        """
        return cls(content)

    @property
    def content(self) -> str:
        """Raw content of the OFF file."""
        return self._content

    def _iter_lines(self) -> Iterator[str]:
        """Return file content as an iterable of lines."""
        return iter(self._content.split("\n"))

    def get_residue(self, unit: str, skip_h: bool = False) -> Residue:
        """
        Fetch residue unit and return a Residue instance with atomic information.

        Parameters
        ----------
        unit : str
            Residue/unit name in the OFF file
        skip_h : bool, optional
            Skip hydrogen atoms, by default False

        Returns
        -------
        Residue
            Residue instance with atoms, connectivity, and coordinates
        """
        atoms = self.get_atoms(unit, skip_h=skip_h)
        connectivity = self.get_connectivity(unit)
        xyz = self.get_coords(unit)

        residue = Residue(name=unit, atoms=atoms, connectivity=connectivity, xyz=xyz)
        return residue

    def get_coords(self, unit: str) -> np.ndarray:
        """
        Fetch position coordinates for a unit.

        This reads the !entry.UNIT.unit.positions section.

        Parameters
        ----------
        unit : str
            Unit name

        Returns
        -------
        np.ndarray
            Coordinates array of shape (N, 3)
        """
        positions = self.read_off_section(unit, "positions")
        xyz = np.array([line.split() for line in positions], dtype=np.float32)
        return xyz

    def get_connectivity(self, unit: str) -> list[tuple[int, int]]:
        """
        Fetch connectivity table for a unit.

        This reads the !entry.UNIT.unit.connectivity section.

        Parameters
        ----------
        unit : str
            Unit name

        Returns
        -------
        list[tuple[int, int]]
            List of bonded atom pairs (bidirectional)
        """
        try:
            connect_info = self.read_off_section(unit, "connectivity")
        except OFFSectionError:
            return []

        connectivity = []
        for line in connect_info:
            parts = line.split()
            if len(parts) >= 2:
                pair = (int(parts[0]), int(parts[1]))
                connectivity.append(pair)
                connectivity.append((pair[1], pair[0]))  # Add reverse direction
        return connectivity

    def get_atoms(self, unit: str, skip_h: bool = False) -> list[Atom]:
        """
        Fetch atomic information for a unit.

        The atoms section format has columns:
        name, type, typex, reession, residuePdbSequenceNumber, flags, sequence, elmnt, charge

        Parameters
        ----------
        unit : str
            Unit name
        skip_h : bool, optional
            Skip hydrogen atoms (element=1), by default False

        Returns
        -------
        list[Atom]
            List of Atom instances
        """
        atom_list = []
        atom_info = self.read_off_section(unit, "atoms")

        for line in atom_info:
            parts = line.split()
            if len(parts) < 9:
                continue

            # Column indices based on OFF format:
            # 0=name, 1=type, 2=typex, 3=reession, 4=residuePdbSeq,
            # 5=flags, 6=sequence, 7=elmnt, 8=charge
            element = int(parts[7])
            if skip_h and element == 1:
                continue

            atom_id = int(parts[6])  # sequence column
            name = parts[0].strip('"')
            atom_type = parts[1].strip('"')
            charge = float(parts[8])

            atom = Atom(
                id=atom_id,
                name=name,
                atom_type=atom_type,
                element=element,
                charge=charge,
            )
            atom_list.append(atom)

        return atom_list

    def get_units(self) -> list[str]:
        """
        Get list of unit names in the OFF file.

        Returns
        -------
        list[str]
            Unit names
        """
        lines = self._iter_lines()
        # Skip first line (file header)
        next(lines)

        units = []
        for line in lines:
            line = line.strip()
            if "!" in line:
                break
            if '"' in line:
                units.append(line.split('"')[1])

        return units

    def has_unit(self, unit_name: str) -> bool:
        """
        Check if the OFF file contains a specific unit.

        Parameters
        ----------
        unit_name : str
            Unit name to check

        Returns
        -------
        bool
            True if unit exists
        """
        return unit_name in self.get_units()

    def is_parameter(self, unit: str) -> bool:
        """
        Check if a unit is a parameter unit (not a molecule/residue).

        Parameters
        ----------
        unit : str
            Unit name

        Returns
        -------
        bool
            True if it's a parameter unit
        """
        matches = self._find(f"!entry.{unit}.")
        if matches:
            return "parm" in matches[0]
        return False

    def _find(self, expr: str, return_lines: bool = True) -> list[str] | int | bool:
        """
        Find expression in OFF file using fnmatch wildcards.

        Parameters
        ----------
        expr : str
            Expression to search for (wildcards allowed)
        return_lines : bool, optional
            Return matching lines if True, else count

        Returns
        -------
        list[str] | int | bool
            Matching lines, count, or False if no match
        """
        import fnmatch

        lines = list(self._iter_lines())
        matches = fnmatch.filter(lines, f"*{expr}*")

        if not matches:
            return False
        if return_lines:
            return matches
        return len(matches)

    def get_residue_list(self, unit: str, unique: bool = True) -> list[str]:
        """
        Get list of residue names within a unit.

        Parameters
        ----------
        unit : str
            Unit name
        unique : bool, optional
            Return unique names only, by default True

        Returns
        -------
        list[str]
            Residue names
        """
        section = self.read_off_section(unit, "residues")
        res_names = []
        for line in section:
            parts = line.split()
            if parts:
                res_names.append(parts[0].replace('"', ""))

        if unique:
            return list(set(res_names))
        return res_names

    def read_off_section(self, unit: str, section: str, with_header: bool = False) -> list[str]:
        """
        Read a complete section from the OFF file for a given unit.

        Parameters
        ----------
        unit : str
            Unit name
        section : str
            Section name (e.g., 'residues', 'atoms', 'positions')
        with_header : bool, optional
            Include section header line, by default False

        Returns
        -------
        list[str]
            Section content lines

        Raises
        ------
        OFFSectionError
            If section not found
        """
        search = f"!entry.{unit}.unit.{section}"
        lines = self._iter_lines()

        # Find section start
        line = next(lines, None)
        while line is not None and search not in line:
            line = next(lines, None)

        if line is None:
            raise OFFSectionError(f"Section {section} for unit {unit} not in file.")

        out = []
        if with_header:
            out.append(line)

        # Read until next section
        line = next(lines, None)
        while line is not None and "!entry" not in line:
            out.append(line.strip())
            line = next(lines, None)

        return out

    def get_num_res(self, unit: str, residue: str | None = None) -> int:
        """
        Count residues in a unit.

        Parameters
        ----------
        unit : str
            Unit name
        residue : str | None, optional
            If given, count only this residue type

        Returns
        -------
        int
            Residue count
        """
        res_list = self.get_residue_list(unit, unique=False)
        if residue is None:
            return len(res_list)
        return res_list.count(residue)

    def get_num_atoms(self, unit: str, residue: str, atom_name: str) -> int:
        """
        Count atoms with a specific name in a unit.

        Parameters
        ----------
        unit : str
            Unit name
        residue : str
            Residue name
        atom_name : str
            Atom name to count

        Returns
        -------
        int
            Atom count
        """
        n_res = self.get_num_res(unit, residue)
        section = self.read_off_section(residue, "atoms")
        n = sum(1 for line in section if atom_name in line)
        return n * n_res

    def get_box_dimensions(self, unit: str) -> list[float]:
        """
        Get box dimensions from the OFF file.

        Parameters
        ----------
        unit : str
            Unit name

        Returns
        -------
        list[float]
            Box dimensions [x, y, z] (or [angle, x, y, z] depending on format)
        """
        box_section = self.read_off_section(unit, "boundbox")
        # Skip first two lines (box type info), read dimensions
        return [float(x) for x in box_section[2:]]

    def get_volume(self, unit: str) -> float | None:
        """
        Calculate volume from box dimensions.

        Parameters
        ----------
        unit : str
            Unit name

        Returns
        -------
        float | None
            Box volume, or None if no box dimensions
        """
        try:
            dims = self.get_box_dimensions(unit)
            if not dims:
                return None
            return reduce(lambda x, y: x * y, dims)
        except OFFSectionError:
            return None

    def write(self, outpath: str | Path) -> Path:
        """
        Write OFF content to a file.

        Parameters
        ----------
        outpath : str | Path
            Output file path

        Returns
        -------
        Path
            Absolute path to the written file
        """
        outpath = Path(outpath)
        outpath.write_text(self._content)
        return outpath.absolute()

    def write_tmp(self) -> Path:
        """
        Write OFF content to a temporary file.

        The path is stored in _tmpfile for later cleanup.

        Returns
        -------
        Path
            Path to temporary file
        """
        fd, tmp_path = tempfile.mkstemp(suffix=".off")
        self._tmpfile = Path(tmp_path)
        with open(fd, "w") as f:
            f.write(self._content)
        return self._tmpfile

    def clean_tmp(self) -> bool:
        """
        Remove temporary file created by write_tmp().

        Returns
        -------
        bool
            True if file was removed, False otherwise
        """
        if self._tmpfile and self._tmpfile.exists():
            self._tmpfile.unlink()
            self._tmpfile = None
            return True
        return False

    def __repr__(self) -> str:
        units = self.get_units()
        return f"OFFManager(units={units})"
