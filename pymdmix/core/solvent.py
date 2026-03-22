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

import json
import logging
from dataclasses import dataclass, field
from pathlib import Path

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
    probe_types : list[str]
        Chemical types for this probe (e.g., ["Don", "Acc"], ["Hyd"], ["Wat"])
    probability : float | None
        Expected probability per probe atom (derived from volume)
    """

    name: str
    residue: str
    atoms: list[str]
    description: str = ""
    probe_types: list[str] = field(default_factory=list)
    probability: float | None = None

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
    volume : float | None
        Box volume in Å³
    is_ionic : bool
        Whether the solvent contains ionic species
    total_charge : float
        Total charge of the solvent box
    corrections : dict[str, float]
        Correction factors per probe

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
    ...         Probe("OH", "ETA", ["O1"], "Ethanol hydroxyl", ["Don", "Acc"]),
    ...         Probe("CT", "ETA", ["C2"], "Ethanol methyl", ["Hyd"]),
    ...     ],
    ...     volume=7988.43,
    ...     is_ionic=False,
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
    volume: float | None = None
    is_ionic: bool = False
    total_charge: float = 0.0
    corrections: dict[str, float] = field(default_factory=dict)

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

    def get_types_map(self) -> dict[str, str]:
        """
        Get mapping of probe names to chemical types.

        Returns a dict like {'OH': 'Don,Acc', 'CT': 'Hyd', 'WAT': 'Wat'}
        """
        return {p.name: ",".join(p.probe_types) for p in self.probes if p.probe_types}

    def calculate_probability(self, probe_name: str) -> float | None:
        """
        Calculate expected probability for a probe.

        Probability = 1 / volume (per probe atom in standard units).
        This is used for Boltzmann energy calculations.

        Parameters
        ----------
        probe_name : str
            Name of the probe

        Returns
        -------
        float | None
            Probability value, or None if volume not set or probe not found
        """
        if self.volume is None or self.volume <= 0:
            return None

        probe = self.get_probe(probe_name)
        if probe is None:
            return None

        # Return stored probability if available
        if probe.probability is not None:
            return probe.probability

        # Calculate from volume (simplified - actual calculation may vary)
        return 1.0 / self.volume

    @classmethod
    def from_dict(cls, data: dict) -> Solvent:
        """Create Solvent from dictionary."""
        residues = [
            SolventResidue(**r) if isinstance(r, dict) else r for r in data.get("residues", [])
        ]

        # Handle probes with new fields
        probes = []
        for p in data.get("probes", []):
            if isinstance(p, dict):
                probe = Probe(
                    name=p["name"],
                    residue=p["residue"],
                    atoms=p["atoms"],
                    description=p.get("description", ""),
                    probe_types=p.get("probe_types", []),
                    probability=p.get("probability"),
                )
                probes.append(probe)
            else:
                probes.append(p)

        off_file = data.get("off_file")
        if off_file:
            off_file = Path(off_file)

        return cls(
            name=data["name"],
            description=data.get("description", ""),
            full_name=data.get("full_name", data["name"]),
            residues=residues,
            probes=probes,
            off_file=off_file,
            box_unit=data.get("box_unit", ""),
            water_model=data.get("water_model", "TIP3P"),
            volume=data.get("volume"),
            is_ionic=data.get("is_ionic", False),
            total_charge=data.get("total_charge", 0.0),
            corrections=data.get("corrections", {}),
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
        result = {
            "name": self.name,
            "full_name": self.full_name,
            "description": self.description,
            "off_file": str(self.off_file) if self.off_file else None,
            "box_unit": self.box_unit,
            "water_model": self.water_model,
            "volume": self.volume,
            "is_ionic": self.is_ionic,
            "total_charge": self.total_charge,
            "corrections": self.corrections,
            "residues": [
                {"name": r.name, "count": r.count, "description": r.description}
                for r in self.residues
            ],
            "probes": [
                {
                    "name": p.name,
                    "residue": p.residue,
                    "atoms": p.atoms,
                    "description": p.description,
                    "probe_types": p.probe_types,
                    "probability": p.probability,
                }
                for p in self.probes
            ],
        }
        return result

    def to_json(self, path: str | Path) -> None:
        """Save solvent to JSON file."""
        with open(path, "w") as f:
            json.dump(self.to_dict(), f, indent=2)

    def write_off(self, path: str | Path) -> int:
        """
        Write the solvent OFF file to disk.

        This preserves the legacy `Solvent.writeOff()` behavior used by the
        original Amber solvation workflow, where the solvent object could emit
        a temporary OFF file for LEaP loading.

        Parameters
        ----------
        path : str | Path
            Destination OFF file path.

        Returns
        -------
        int
            Number of characters written.

        Raises
        ------
        ValueError
            If this solvent has no associated OFF file.
        FileNotFoundError
            If the configured OFF file does not exist.
        """
        if self.off_file is None:
            raise ValueError(f"Solvent {self.name} has no OFF file configured")
        if not self.off_file.exists():
            raise FileNotFoundError(f"Solvent OFF file not found: {self.off_file}")

        content = self.off_file.read_text()
        return Path(path).write_text(content)

    @classmethod
    def from_file(cls, path: str | Path) -> Solvent:
        """
        Load a Solvent from a file, dispatching by extension.

        Supports:
        - ``.json`` — native JSON format (``from_json``)
        - ``.cfg`` / ``.ini`` — legacy pyMDMix INI format

        Parameters
        ----------
        path : str or Path
            Path to solvent definition file

        Returns
        -------
        Solvent
            Loaded solvent definition

        Examples
        --------
        >>> solvent = Solvent.from_file("my_solvent.json")
        >>> solvent = Solvent.from_file("my_solvent.cfg")
        """
        path = Path(path)
        if not path.exists():
            raise FileNotFoundError(f"Solvent file not found: {path}")

        if path.suffix == ".json":
            return cls.from_json(path)
        elif path.suffix in (".cfg", ".ini"):
            return cls._from_cfg(path)
        else:
            # Attempt JSON by default
            return cls.from_json(path)

    @classmethod
    def _from_cfg(cls, path: Path) -> Solvent:
        """
        Load Solvent from legacy pyMDMix ``.cfg`` INI format.

        Expected sections: ``[SOLVENT]``, ``[PROBES]``, ``[RESIDUES]``
        """
        import configparser

        parser = configparser.ConfigParser()
        parser.read(str(path))

        # --- [SOLVENT] section ---
        solv_sec = dict(parser["SOLVENT"]) if parser.has_section("SOLVENT") else {}
        name = solv_sec.get("name", path.stem.upper())
        description = solv_sec.get("description", "")
        full_name = solv_sec.get("full_name", name)
        off_file_str = solv_sec.get("off_file")
        box_unit = solv_sec.get("box_unit", "")
        volume = float(solv_sec["volume"]) if "volume" in solv_sec else None
        is_ionic = solv_sec.get("is_ionic", "false").lower() in ("true", "1", "yes")
        water_model = solv_sec.get("water_model", "TIP3P")

        # --- [RESIDUES] section ---
        residues: list[SolventResidue] = []
        if parser.has_section("RESIDUES"):
            for res_name, count_str in parser.items("RESIDUES"):
                if res_name.lower() == "default":
                    continue
                try:
                    residues.append(SolventResidue(name=res_name.upper(), count=int(count_str)))
                except ValueError:
                    pass

        # --- [PROBES] / per-probe sections ---
        probes: list[Probe] = []
        probe_names = [k for k in solv_sec if k.startswith("probe")] + (
            [p.strip() for p in solv_sec.get("probes", "").split(",") if p.strip()]
        )
        # Deduplicate preserving order
        seen: set[str] = set()
        unique_probes = [p for p in probe_names if not (p in seen or seen.add(p))]  # type: ignore[func-returns-value]

        for probe_name in unique_probes:
            sec = probe_name.upper()
            if parser.has_section(sec):
                pd = dict(parser[sec])
                atoms_raw = pd.get("atoms", pd.get("atom", ""))
                atoms = [a.strip() for a in atoms_raw.split(",") if a.strip()]
                types_raw = pd.get("types", pd.get("type", ""))
                probe_types = [t.strip() for t in types_raw.split(",") if t.strip()]
                probes.append(
                    Probe(
                        name=sec,
                        residue=pd.get("residue", ""),
                        atoms=atoms,
                        description=pd.get("description", ""),
                        probe_types=probe_types,
                        probability=float(pd["probability"]) if "probability" in pd else None,
                    )
                )

        off_file = Path(off_file_str) if off_file_str else None
        if off_file and not off_file.is_absolute():
            off_file = path.parent / off_file

        return cls(
            name=name,
            description=description,
            full_name=full_name,
            residues=residues,
            probes=probes,
            off_file=off_file,
            box_unit=box_unit,
            water_model=water_model,
            volume=volume,
            is_ionic=is_ionic,
        )

    def __repr__(self) -> str:
        return f"Solvent({self.name!r}, residues={len(self.residues)}, probes={len(self.probes)})"


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
        return Path(__file__).parent.parent / "data" / "solvents"

    def _load_library(self) -> None:
        """Load all solvents from library directory."""
        if not self.library_path.exists():
            log.warning(f"Solvent library path does not exist: {self.library_path}")
            return

        for json_file in self.library_path.glob("*.json"):
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
            probes=[Probe("O", "WAT", ["O"], "Water oxygen", ["Wat"], 0.00408)],
            volume=6617.51,
            is_ionic=False,
        ),
        Solvent(
            name="ETA",
            description="20% Ethanol / 80% Water",
            residues=[
                SolventResidue("ETA", 200, "Ethanol"),
                SolventResidue("WAT", 800, "Water"),
            ],
            probes=[
                Probe("OH", "ETA", ["O1"], "Ethanol hydroxyl oxygen", ["Don", "Acc"], 0.000266),
                Probe("CT", "ETA", ["C1"], "Ethanol methyl carbon", ["Hyd"], 0.000266),
                Probe("WAT", "WAT", ["O"], "Water oxygen", ["Wat"], 0.00327),
            ],
            volume=7988.43,
            is_ionic=False,
        ),
        Solvent(
            name="MAM",
            description="20% Acetamide / 80% Water",
            residues=[
                SolventResidue("MAM", 200, "Acetamide"),
                SolventResidue("WAT", 800, "Water"),
            ],
            probes=[
                Probe("N", "MAM", ["N1"], "Amine nitrogen", ["Don"], 0.000288),
                Probe("O", "MAM", ["O1"], "Carbonyl oxygen", ["Acc"], 0.000288),
                Probe("CT", "MAM", ["C2"], "Methyl carbon", ["Hyd"], 0.000288),
                Probe("WAT", "WAT", ["O"], "Water oxygen", ["Wat"], 0.00334),
            ],
            volume=7817.91,
            is_ionic=False,
        ),
        Solvent(
            name="ANT",
            description="20% Acetonitrile / 80% Water",
            residues=[
                SolventResidue("ANT", 200, "Acetonitrile"),
                SolventResidue("WAT", 800, "Water"),
            ],
            probes=[
                Probe("N", "ANT", ["N1"], "Nitrile nitrogen", ["Acc"], 0.000271),
                Probe("C", "ANT", ["C3"], "Methyl carbon", ["Hyd"], 0.000271),
                Probe("WAT", "WAT", ["O"], "Water oxygen", ["Wat"], 0.00333),
            ],
            volume=7834.51,
            is_ionic=False,
        ),
    ]
