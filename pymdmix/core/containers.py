"""
Data Container Classes
======================

Simple container classes for atomic and residue information.

Migrated from the original pyMDMix/containers.py.

Note: The Probe class is now in pymdmix.core.solvent.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass
class Atom:
    """
    Container for atomic information from OFF/MOL2 files.

    Attributes
    ----------
    id : int
        Atom ID (1-indexed).
    name : str
        Atom name.
    type : str
        Atom type (AMBER type).
    element : int | str
        Element (atomic number or symbol).
    charge : float
        Partial charge.

    Examples
    --------
    >>> atom = Atom(id=1, name='CA', type='CT', element=6, charge=-0.1234)
    >>> atom.name
    'CA'
    """

    id: int
    name: str
    type: str
    element: int | str
    charge: float

    def __repr__(self) -> str:
        return self.name

    def __str__(self) -> str:
        return (
            f"ID: {self.id}\n"
            f"ATOM NAME: {self.name}\n"
            f"ATOM TYPE: {self.type}\n"
            f"ELEMENT: {self.element}\n"
            f"CHARGE: {self.charge:.4f}"
        )

    def __eq__(self, other: object) -> bool:
        """Compare by name."""
        if isinstance(other, Atom):
            return self.name == other.name
        if isinstance(other, str):
            return self.name == other
        return NotImplemented

    def __hash__(self) -> int:
        return hash(self.name)


@dataclass
class Residue:
    """
    Container for residue information from OFF/MOL2 files.

    Attributes
    ----------
    name : str
        Residue name.
    atoms : list[Atom]
        List of atoms in the residue.
    connectivity : tuple | list
        Connectivity information (bonds).
    xyz : np.ndarray
        Coordinates (N x 3 array).

    Examples
    --------
    >>> atoms = [Atom(1, 'C1', 'CT', 6, 0.0), Atom(2, 'C2', 'CT', 6, 0.0)]
    >>> xyz = np.array([[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]])
    >>> residue = Residue(name='ETA', atoms=atoms, connectivity=[(0, 1)], xyz=xyz)
    >>> residue.charge
    0.0
    """

    name: str
    atoms: list[Atom]
    connectivity: tuple | list = field(default_factory=list)
    xyz: np.ndarray = field(default_factory=lambda: np.array([]))

    # Computed fields
    _charge: float = field(init=False, repr=False, default=0.0)
    _atids: dict[int, Atom] = field(init=False, repr=False, default_factory=dict)
    _atnames: dict[str, Atom] = field(init=False, repr=False, default_factory=dict)

    def __post_init__(self):
        """Calculate derived properties."""
        # Calculate total charge
        if self.atoms:
            self._charge = round(sum(a.charge for a in self.atoms), 4)

        # Build lookup maps
        self._atids = {a.id: a for a in self.atoms}
        self._atnames = {a.name: a for a in self.atoms}

    @property
    def charge(self) -> float:
        """Total charge of the residue."""
        return self._charge

    @property
    def atids(self) -> dict[int, Atom]:
        """Dictionary mapping atom IDs to Atom instances."""
        return self._atids

    @property
    def atnames(self) -> dict[str, Atom]:
        """Dictionary mapping atom names to Atom instances."""
        return self._atnames

    def __repr__(self) -> str:
        return self.name

    def __str__(self) -> str:
        return f"RESIDUE NAME: {self.name}\nATOMS: {self.atoms}"

    def __eq__(self, other: object) -> bool:
        """Compare by name."""
        if isinstance(other, Residue):
            return self.name == other.name
        if isinstance(other, str):
            return self.name == other
        return NotImplemented

    def __hash__(self) -> int:
        return hash(self.name)

    def center(self) -> np.ndarray:
        """
        Get geometric center of residue.

        Returns
        -------
        np.ndarray
            Center coordinates (x, y, z).
        """
        if len(self.xyz) == 0:
            return np.array([0.0, 0.0, 0.0])  # type: ignore[no-any-return]
        return self.xyz.mean(axis=0)  # type: ignore[no-any-return]

    def set_xyz(self, xyz: np.ndarray) -> bool:
        """
        Set new coordinates.

        Parameters
        ----------
        xyz : np.ndarray
            New coordinates (must match shape).

        Returns
        -------
        bool
            True if successful.
        """
        if xyz.shape == self.xyz.shape:
            self.xyz = xyz
            return True
        return False

    def get_atom(self, name: str) -> Atom | None:
        """
        Get atom by name.

        Parameters
        ----------
        name : str
            Atom name.

        Returns
        -------
        Atom | None
            Atom if found, None otherwise.
        """
        return self._atnames.get(name)

    def get_atom_by_id(self, atom_id: int) -> Atom | None:
        """
        Get atom by ID.

        Parameters
        ----------
        atom_id : int
            Atom ID.

        Returns
        -------
        Atom | None
            Atom if found, None otherwise.
        """
        return self._atids.get(atom_id)


if __name__ == "__main__":
    # Test
    atoms = [
        Atom(id=1, name="C1", type="CT", element=6, charge=0.123),
        Atom(id=2, name="O1", type="OH", element=8, charge=-0.456),
        Atom(id=3, name="H1", type="HO", element=1, charge=0.333),
    ]

    xyz = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.5, 0.0, 0.0],
            [2.0, 0.5, 0.0],
        ]
    )

    residue = Residue(name="ETA", atoms=atoms, connectivity=[(0, 1), (1, 2)], xyz=xyz)

    print(f"Residue: {residue}")
    print(f"Total charge: {residue.charge}")
    print(f"Center: {residue.center()}")
    print(f"Atom C1: {residue.get_atom('C1')}")
    print(f"Atom by ID 2: {residue.get_atom_by_id(2)}")
