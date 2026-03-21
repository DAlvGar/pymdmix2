"""
Replica Management
==================

A Replica represents a single MD simulation run.
Multiple replicas with different random seeds improve sampling.

Replicas track:
- Setup state (topology, coordinates)
- Simulation progress
- Analysis status

Examples
--------
>>> from pymdmix.project import Replica, ReplicaState
>>>
>>> replica = Replica(
...     name="ETA_1",
...     solvent="ETA",
...     path=Path("replicas/ETA_1"),
... )
>>>
>>> # Check state
>>> if replica.state == ReplicaState.READY:
...     replica.submit()
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum, auto
from pathlib import Path
from typing import Any
import json
import logging
from datetime import datetime

log = logging.getLogger(__name__)


class ReplicaState(Enum):
    """
    State of a replica in the workflow.

    States progress: CREATED → SETUP → READY → RUNNING → COMPLETE → ANALYZED
    """
    CREATED = auto()    # Just created, no files
    SETUP = auto()      # Topology/coordinates generated
    READY = auto()      # Ready to submit
    RUNNING = auto()    # Simulation in progress
    COMPLETE = auto()   # Simulation finished
    ALIGNED = auto()    # Trajectory aligned
    ANALYZED = auto()   # Analysis complete
    ERROR = auto()      # Error state


@dataclass
class Replica:
    """
    A single MD simulation replica.

    Attributes
    ----------
    name : str
        Replica identifier (e.g., "ETA_1")
    solvent : str
        Solvent mixture name
    path : Path | None
        Directory containing replica files
    state : ReplicaState
        Current workflow state
    topology : str | None
        Topology filename (relative to path)
    coordinates : str | None
        Coordinate filename
    trajectory : str | None
        Trajectory filename
    reference : str | None
        Reference PDB for alignment
    seed : int
        Random seed for simulation
    n_steps : int
        Total steps completed
    """
    name: str
    solvent: str
    path: Path | None = None
    state: ReplicaState = ReplicaState.CREATED

    # Files (relative to path)
    topology: str | None = None
    coordinates: str | None = None
    trajectory: str | None = None
    reference: str | None = None

    # Simulation info
    seed: int = 0
    n_steps: int = 0

    # Timestamps
    created_at: str = field(default_factory=lambda: datetime.now().isoformat())
    modified_at: str = field(default_factory=lambda: datetime.now().isoformat())

    # Metadata
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self):
        """Normalize path."""
        if self.path and isinstance(self.path, str):
            self.path = Path(self.path)

        if isinstance(self.state, str):
            self.state = ReplicaState[self.state]

    @property
    def topology_path(self) -> Path | None:
        """Full path to topology file."""
        if self.path and self.topology:
            return self.path / self.topology
        return None

    @property
    def coordinates_path(self) -> Path | None:
        """Full path to coordinates file."""
        if self.path and self.coordinates:
            return self.path / self.coordinates
        return None

    @property
    def trajectory_path(self) -> Path | None:
        """Full path to trajectory file."""
        if self.path and self.trajectory:
            return self.path / self.trajectory
        return None

    @property
    def reference_path(self) -> Path | None:
        """Full path to reference PDB."""
        if self.path and self.reference:
            return self.path / self.reference
        return None

    def exists(self) -> bool:
        """Check if replica directory exists."""
        return self.path is not None and self.path.exists()

    def has_topology(self) -> bool:
        """Check if topology file exists."""
        return self.topology_path is not None and self.topology_path.exists()

    def has_trajectory(self) -> bool:
        """Check if trajectory file exists."""
        return self.trajectory_path is not None and self.trajectory_path.exists()

    def create_directory(self) -> None:
        """Create replica directory structure."""
        if self.path is None:
            raise ValueError("Replica path not set")

        self.path.mkdir(parents=True, exist_ok=True)

        # Create subdirectories
        (self.path / "input").mkdir(exist_ok=True)
        (self.path / "output").mkdir(exist_ok=True)
        (self.path / "analysis").mkdir(exist_ok=True)

        self._update_modified()
        log.info(f"Created replica directory: {self.path}")

    def set_state(self, state: ReplicaState) -> None:
        """Update replica state."""
        old_state = self.state
        self.state = state
        self._update_modified()
        log.info(f"Replica {self.name}: {old_state.name} → {state.name}")

    def _update_modified(self) -> None:
        """Update modified timestamp."""
        self.modified_at = datetime.now().isoformat()

    def to_dict(self) -> dict[str, Any]:
        """Convert to dictionary for serialization."""
        return {
            "name": self.name,
            "solvent": self.solvent,
            "path": str(self.path) if self.path else None,
            "state": self.state.name,
            "topology": self.topology,
            "coordinates": self.coordinates,
            "trajectory": self.trajectory,
            "reference": self.reference,
            "seed": self.seed,
            "n_steps": self.n_steps,
            "created_at": self.created_at,
            "modified_at": self.modified_at,
            "metadata": self.metadata,
        }

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> Replica:
        """Create Replica from dictionary."""
        # Handle state enum
        if "state" in data and isinstance(data["state"], str):
            data["state"] = ReplicaState[data["state"]]

        # Handle path
        if "path" in data and data["path"]:
            data["path"] = Path(data["path"])

        return cls(**data)

    def save(self, path: Path | None = None) -> None:
        """
        Save replica state to JSON file.

        Parameters
        ----------
        path : Path | None
            Output path. Defaults to {self.path}/replica.json
        """
        if path is None:
            if self.path is None:
                raise ValueError("No path specified and replica path not set")
            path = self.path / "replica.json"

        with open(path, 'w') as f:
            json.dump(self.to_dict(), f, indent=2)

        log.debug(f"Saved replica to {path}")

    @classmethod
    def load(cls, path: Path) -> Replica:
        """
        Load replica from JSON file.

        Parameters
        ----------
        path : Path
            Path to replica.json file

        Returns
        -------
        Replica
            Loaded replica
        """
        with open(path) as f:
            data = json.load(f)

        # Set path from file location if not in data
        if data.get("path") is None:
            data["path"] = path.parent

        return cls.from_dict(data)

    def __repr__(self) -> str:
        return (
            f"Replica(name={self.name!r}, solvent={self.solvent!r}, "
            f"state={self.state.name})"
        )


def create_replica(
    name: str,
    solvent: str,
    base_path: Path,
    seed: int | None = None,
) -> Replica:
    """
    Create a new replica with directory structure.

    Parameters
    ----------
    name : str
        Replica name
    solvent : str
        Solvent mixture name
    base_path : Path
        Base directory for replicas
    seed : int | None
        Random seed (auto-generated if None)

    Returns
    -------
    Replica
        New replica with created directory
    """
    import random

    if seed is None:
        seed = random.randint(1, 999999)

    replica = Replica(
        name=name,
        solvent=solvent,
        path=base_path / name,
        seed=seed,
    )

    replica.create_directory()
    replica.save()

    return replica
