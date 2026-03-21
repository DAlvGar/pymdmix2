"""
Solvent Mixture Definitions
===========================

Defines solvent mixtures used in MDMix simulations.
Each solvent has residue definitions and probe atoms for analysis.

Examples
--------
>>> library = SolventLibrary()
>>> ethanol = library.get("ETA")
>>> print(ethanol.probes)
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
import json
import logging

log = logging.getLogger(__name__)


@dataclass
class Probe:
    """
    Probe atom(s) to track in density analysis.

    Attributes
    ----------
    name : str
        Probe identifier (e.g., "OH", "CT", "WAT")
    residue : str
        Parent residue name
    atoms : list[str]
        Atom names to track
    description : str
        Human-readable description
    """
    name: str
    residue: str
    atoms: list[str]
    description: str = ""

    @property
    def selection(self) -> str:
        """MDAnalysis selection string for this probe."""
        atoms_str = " ".join(self.atoms)
        return f"resname {self.residue} and name {atoms_str}"

    def __repr__(self) -> str:
        return f"Probe({self.name!r}, residue={self.residue!r}, atoms={self.atoms})"


@dataclass
class SolventResidue:
    """
    Residue component of a solvent mixture.

    Attributes
    ----------
    name : str
        Residue name (as in topology)
    count : int
        Number of this residue in solvent box
    description : str
        Human-readable description
    """
    name: str
    count: int = 0
    description: str = ""

    def __repr__(self) -> str:
        return f"SolventResidue({self.name!r}, count={self.count})"


@dataclass
class Solvent:
    """
    Solvent mixture definition.

    Attributes
    ----------
    name : str
        Solvent identifier (e.g., "ETA", "MAM")
    description : str
        Human-readable description
    residues : list[SolventResidue]
        Residue components
    probes : list[Probe]
        Probe definitions for analysis
    off_file : Path | None
        LEaP object file for solvation

    Examples
    --------
    >>> solvent = Solvent(
    ...     name="ETA",
    ...     description="20% Ethanol / 80% Water",
    ...     residues=[
    ...         SolventResidue("ETA", 200),
    ...         SolventResidue("WAT", 800),
    ...     ],
    ...     probes=[
    ...         Probe("OH", "ETA", ["O1"], "Ethanol hydroxyl"),
    ...         Probe("CT", "ETA", ["C2"], "Ethanol methyl"),
    ...     ],
    ... )
    """
    name: str
    description: str = ""
    full_name: str = ""
    residues: list[SolventResidue] = field(default_factory=list)
    probes: list[Probe] = field(default_factory=list)
    off_file: Path | None = None
    box_unit: str = ""
    water_model: str = "TIP3P"

    def __post_init__(self):
        """Set full_name to name if not provided."""
        if not self.full_name:
            self.full_name = self.name

    def get_probe(self, name: str) -> Probe | None:
        """Get probe by name."""
        for p in self.probes:
            if p.name == name:
                return p
        return None

    def get_residue_names(self) -> list[str]:
        """Get list of residue names in this solvent."""
        return [r.name for r in self.residues]

    def get_probe_names(self) -> list[str]:
        """Get list of probe names."""
        return [p.name for p in self.probes]

    @classmethod
    def from_dict(cls, data: dict) -> Solvent:
        """Create Solvent from dictionary."""
        residues = [
            SolventResidue(**r) if isinstance(r, dict) else r
            for r in data.get('residues', [])
        ]
        probes = [
            Probe(**p) if isinstance(p, dict) else p
            for p in data.get('probes', [])
        ]

        off_file = data.get('off_file')
        if off_file:
            off_file = Path(off_file)

        return cls(
            name=data['name'],
            description=data.get('description', ''),
            full_name=data.get('full_name', data['name']),
            residues=residues,
            probes=probes,
            off_file=off_file,
            box_unit=data.get('box_unit', ''),
            water_model=data.get('water_model', 'TIP3P'),
        )

    @classmethod
    def from_json(cls, path: str | Path) -> Solvent:
        """Load solvent from JSON file."""
        with open(path) as f:
            data = json.load(f)

        solvent = cls.from_dict(data)

        # Resolve off_file path relative to JSON file
        if solvent.off_file and not solvent.off_file.is_absolute():
            solvent.off_file = Path(path).parent / solvent.off_file

        return solvent

    def to_dict(self) -> dict:
        """Convert to dictionary."""
        return {
            'name': self.name,
            'description': self.description,
            'residues': [
                {'name': r.name, 'count': r.count, 'description': r.description}
                for r in self.residues
            ],
            'probes': [
                {'name': p.name, 'residue': p.residue,
                 'atoms': p.atoms, 'description': p.description}
                for p in self.probes
            ],
            'off_file': str(self.off_file) if self.off_file else None,
        }

    def to_json(self, path: str | Path) -> None:
        """Save solvent to JSON file."""
        with open(path, 'w') as f:
            json.dump(self.to_dict(), f, indent=2)

    def __repr__(self) -> str:
        return (
            f"Solvent({self.name!r}, "
            f"residues={len(self.residues)}, probes={len(self.probes)})"
        )


class SolventLibrary:
    """
    Manager for available solvent mixtures.

    Loads solvents from JSON files in the library directory.

    Parameters
    ----------
    library_path : Path | None
        Path to solvent library directory.
        Defaults to package data/solvents directory.

    Examples
    --------
    >>> library = SolventLibrary()
    >>> print(library.list_solvents())
    ['ETA', 'MAM', 'WAT', ...]
    >>> ethanol = library.get("ETA")
    """

    def __init__(self, library_path: Path | None = None):
        self.library_path = library_path or self._default_path()
        self._solvents: dict[str, Solvent] = {}
        self._load_library()

    def _default_path(self) -> Path:
        """Get default library path (package data)."""
        return Path(__file__).parent.parent / 'data' / 'solvents'

    def _load_library(self) -> None:
        """Load all solvents from library directory."""
        if not self.library_path.exists():
            log.warning(f"Solvent library path does not exist: {self.library_path}")
            return

        for json_file in self.library_path.glob('*.json'):
            try:
                solvent = Solvent.from_json(json_file)
                self._solvents[solvent.name] = solvent
                log.debug(f"Loaded solvent: {solvent.name}")
            except Exception as e:
                log.warning(f"Failed to load solvent from {json_file}: {e}")

    def get(self, name: str) -> Solvent | None:
        """
        Get solvent by name.

        Parameters
        ----------
        name : str
            Solvent name

        Returns
        -------
        Solvent | None
            Solvent if found, None otherwise
        """
        return self._solvents.get(name)

    def list_solvents(self) -> list[str]:
        """Get list of available solvent names."""
        return sorted(self._solvents.keys())

    def add(self, solvent: Solvent) -> None:
        """Add or update a solvent in the library."""
        self._solvents[solvent.name] = solvent

    def save_all(self) -> None:
        """Save all solvents to library directory."""
        self.library_path.mkdir(parents=True, exist_ok=True)
        for name, solvent in self._solvents.items():
            path = self.library_path / f"{name.lower()}.json"
            solvent.to_json(path)
            log.debug(f"Saved solvent {name} to {path}")

    def __len__(self) -> int:
        return len(self._solvents)

    def __iter__(self):
        return iter(self._solvents.values())

    def __repr__(self) -> str:
        return f"SolventLibrary(n_solvents={len(self)})"


# Pre-defined standard solvents
def create_standard_solvents() -> list[Solvent]:
    """Create standard MDMix solvent definitions."""
    return [
        Solvent(
            name="WAT",
            description="Pure water",
            residues=[SolventResidue("WAT", 1000)],
            probes=[Probe("O", "WAT", ["O"], "Water oxygen")],
        ),
        Solvent(
            name="ETA",
            description="20% Ethanol / 80% Water",
            residues=[
                SolventResidue("ETA", 200, "Ethanol"),
                SolventResidue("WAT", 800, "Water"),
            ],
            probes=[
                Probe("OH", "ETA", ["O1"], "Ethanol hydroxyl oxygen"),
                Probe("CT", "ETA", ["C2"], "Ethanol methyl carbon"),
                Probe("WAT", "WAT", ["O"], "Water oxygen"),
            ],
        ),
        Solvent(
            name="MAM",
            description="20% Methylamine / 80% Water",
            residues=[
                SolventResidue("MAM", 200, "Methylamine"),
                SolventResidue("WAT", 800, "Water"),
            ],
            probes=[
                Probe("N", "MAM", ["N"], "Amine nitrogen"),
                Probe("CT", "MAM", ["C"], "Methyl carbon"),
                Probe("WAT", "WAT", ["O"], "Water oxygen"),
            ],
        ),
        Solvent(
            name="ACN",
            description="20% Acetonitrile / 80% Water",
            residues=[
                SolventResidue("ACN", 200, "Acetonitrile"),
                SolventResidue("WAT", 800, "Water"),
            ],
            probes=[
                Probe("N", "ACN", ["N1"], "Nitrile nitrogen"),
                Probe("CT", "ACN", ["C1"], "Methyl carbon"),
                Probe("WAT", "WAT", ["O"], "Water oxygen"),
            ],
        ),
    ]
