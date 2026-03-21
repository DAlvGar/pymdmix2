"""
Project Management
==================

A Project is the top-level container for an MDMix study.
It manages multiple replicas across different solvents.

Examples
--------
>>> from pymdmix.project import Project, Config
>>>
>>> # Create new project
>>> project = Project.create(
...     name="kinase_study",
...     path=Path("./kinase_study"),
...     config=Config(),
... )
>>>
>>> # Add replicas
>>> project.add_replicas("ETA", n_replicas=3)
>>> project.add_replicas("MAM", n_replicas=3)
>>>
>>> # Save and reload
>>> project.save()
>>> project = Project.load(Path("./kinase_study"))
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterator, Any
import json
import logging
from datetime import datetime

from pymdmix.project.config import Config
from pymdmix.project.replica import Replica, ReplicaState, create_replica

log = logging.getLogger(__name__)


@dataclass
class Project:
    """
    MDMix project container.

    Attributes
    ----------
    name : str
        Project name
    path : Path
        Project directory
    config : Config
        Project configuration
    replicas : list[Replica]
        List of replicas
    pdb_file : str | None
        Input PDB filename
    created_at : str
        Creation timestamp
    """
    name: str
    path: Path
    config: Config = field(default_factory=Config)
    replicas: list[Replica] = field(default_factory=list)
    pdb_file: str | None = None
    created_at: str = field(default_factory=lambda: datetime.now().isoformat())
    modified_at: str = field(default_factory=lambda: datetime.now().isoformat())

    def __post_init__(self):
        """Normalize path."""
        if isinstance(self.path, str):
            self.path = Path(self.path)

    @property
    def replicas_path(self) -> Path:
        """Path to replicas directory."""
        return self.path / "replicas"

    @property
    def input_path(self) -> Path:
        """Path to input files directory."""
        return self.path / "input"

    @property
    def analysis_path(self) -> Path:
        """Path to analysis output directory."""
        return self.path / "analysis"

    @property
    def n_replicas(self) -> int:
        """Total number of replicas."""
        return len(self.replicas)

    @property
    def solvents(self) -> list[str]:
        """List of unique solvents used."""
        return list(set(r.solvent for r in self.replicas))

    def get_replicas_by_solvent(self, solvent: str) -> list[Replica]:
        """Get all replicas for a specific solvent."""
        return [r for r in self.replicas if r.solvent == solvent]

    def get_replicas_by_state(self, state: ReplicaState) -> list[Replica]:
        """Get all replicas in a specific state."""
        return [r for r in self.replicas if r.state == state]

    def get_replica(self, name: str) -> Replica | None:
        """Get replica by name."""
        for r in self.replicas:
            if r.name == name:
                return r
        return None

    def add_replica(self, replica: Replica) -> None:
        """Add a replica to the project."""
        if self.get_replica(replica.name):
            raise ValueError(f"Replica {replica.name} already exists")

        self.replicas.append(replica)
        self._update_modified()
        log.info(f"Added replica: {replica.name}")

    def add_replicas(
        self,
        solvent: str,
        n_replicas: int = 3,
        name_template: str = "{solvent}_{i}",
    ) -> list[Replica]:
        """
        Add multiple replicas for a solvent.

        Parameters
        ----------
        solvent : str
            Solvent name
        n_replicas : int
            Number of replicas to create
        name_template : str
            Name template with {solvent} and {i} placeholders

        Returns
        -------
        list[Replica]
            Created replicas
        """
        created = []

        # Find existing replicas for this solvent to get next index
        existing = self.get_replicas_by_solvent(solvent)
        start_idx = len(existing) + 1

        for i in range(start_idx, start_idx + n_replicas):
            name = name_template.format(solvent=solvent, i=i)

            replica = create_replica(
                name=name,
                solvent=solvent,
                base_path=self.replicas_path,
            )

            self.add_replica(replica)
            created.append(replica)

        return created

    def remove_replica(self, name: str) -> bool:
        """
        Remove a replica from the project.

        Note: Does not delete files.
        """
        for i, r in enumerate(self.replicas):
            if r.name == name:
                del self.replicas[i]
                self._update_modified()
                log.info(f"Removed replica: {name}")
                return True
        return False

    def __iter__(self) -> Iterator[Replica]:
        """Iterate over replicas."""
        return iter(self.replicas)

    def __len__(self) -> int:
        """Number of replicas."""
        return len(self.replicas)

    def _update_modified(self) -> None:
        """Update modified timestamp."""
        self.modified_at = datetime.now().isoformat()

    def create_directories(self) -> None:
        """Create project directory structure."""
        self.path.mkdir(parents=True, exist_ok=True)
        self.replicas_path.mkdir(exist_ok=True)
        self.input_path.mkdir(exist_ok=True)
        self.analysis_path.mkdir(exist_ok=True)

        log.info(f"Created project directories: {self.path}")

    def status(self) -> dict[str, Any]:
        """
        Get project status summary.

        Returns
        -------
        dict
            Status information
        """
        state_counts = {}
        for state in ReplicaState:
            count = len(self.get_replicas_by_state(state))
            if count > 0:
                state_counts[state.name] = count

        return {
            "name": self.name,
            "path": str(self.path),
            "n_replicas": self.n_replicas,
            "solvents": self.solvents,
            "states": state_counts,
            "created_at": self.created_at,
            "modified_at": self.modified_at,
        }

    def print_status(self) -> str:
        """Get formatted status string."""
        status = self.status()

        lines = [
            f"Project: {status['name']}",
            f"Path: {status['path']}",
            f"Replicas: {status['n_replicas']}",
            f"Solvents: {', '.join(status['solvents']) or 'None'}",
            "",
            "State Summary:",
        ]

        for state, count in status['states'].items():
            lines.append(f"  {state}: {count}")

        return "\n".join(lines)

    def to_dict(self) -> dict[str, Any]:
        """Convert to dictionary for serialization."""
        return {
            "name": self.name,
            "path": str(self.path),
            "config": self.config.to_dict(),
            "replicas": [r.to_dict() for r in self.replicas],
            "pdb_file": self.pdb_file,
            "created_at": self.created_at,
            "modified_at": self.modified_at,
        }

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> Project:
        """Create Project from dictionary."""
        config = Config.from_dict(data.pop("config", {}))
        replicas_data = data.pop("replicas", [])
        replicas = [Replica.from_dict(r) for r in replicas_data]

        if "path" in data:
            data["path"] = Path(data["path"])

        return cls(config=config, replicas=replicas, **data)

    def save(self, path: Path | None = None) -> None:
        """
        Save project to JSON file.

        Parameters
        ----------
        path : Path | None
            Output path. Defaults to {self.path}/project.json
        """
        if path is None:
            path = self.path / "project.json"

        # Ensure directory exists
        path.parent.mkdir(parents=True, exist_ok=True)

        with open(path, 'w') as f:
            json.dump(self.to_dict(), f, indent=2)

        # Also save individual replica files
        for replica in self.replicas:
            if replica.path:
                replica.save()

        log.info(f"Saved project to {path}")

    @classmethod
    def load(cls, path: Path) -> Project:
        """
        Load project from directory or JSON file.

        Parameters
        ----------
        path : Path
            Path to project directory or project.json file

        Returns
        -------
        Project
            Loaded project
        """
        if path.is_dir():
            json_path = path / "project.json"
        else:
            json_path = path

        if not json_path.exists():
            raise FileNotFoundError(f"Project file not found: {json_path}")

        with open(json_path) as f:
            data = json.load(f)

        return cls.from_dict(data)

    @classmethod
    def create(
        cls,
        name: str,
        path: Path | None = None,
        config: Config | None = None,
    ) -> Project:
        """
        Create a new project with directory structure.

        Parameters
        ----------
        name : str
            Project name
        path : Path | None
            Project directory. Defaults to ./{name}
        config : Config | None
            Project configuration. Defaults to Config()

        Returns
        -------
        Project
            New project with created directories
        """
        if path is None:
            path = Path.cwd() / name

        if config is None:
            config = Config()

        project = cls(name=name, path=path, config=config)
        project.create_directories()
        project.save()

        log.info(f"Created project: {name}")
        return project

    def __repr__(self) -> str:
        return f"Project(name={self.name!r}, replicas={self.n_replicas})"
