"""
Project Management
==================

Classes for managing MDMix projects, replicas, and configuration.

A Project contains multiple Systems, each with multiple Replicas.
Replicas track simulation state and analysis progress.

Examples
--------
>>> from pymdmix.project import Project, Replica, Config
>>>
>>> # Create new project
>>> project = Project.create("my_project", config=Config())
>>>
>>> # Add replicas
>>> replica = Replica(name="ETA_1", solvent="ETA")
>>> project.add_replica(replica)
"""

from pymdmix.project.config import Config, MDSettings
from pymdmix.project.project import Project
from pymdmix.project.replica import Replica, ReplicaState

__all__ = [
    "Config",
    "MDSettings",
    "Replica",
    "ReplicaState",
    "Project",
]
