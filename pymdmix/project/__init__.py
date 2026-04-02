"""
Project Management
==================

Classes for managing MDMix projects, replicas, and configuration.

A Project contains multiple Systems, each with multiple Replicas.
Replicas track simulation state and analysis progress.

Examples
--------
>>> from pymdmix.project import Project, Replica, Config, MDSettings
>>>
>>> # Create new project
>>> project = Project.create("my_project", config=Config())
>>>
>>> # Add replicas with MD settings
>>> settings = MDSettings(nanos=40, temperature=300.0)
>>> replica = Replica(name="ETA_1", solvent="ETA", settings=settings)
>>> project.add_replica(replica)
"""

from pymdmix.project.config import Config
from pymdmix.project.project import Project
from pymdmix.project.replica import MDSettings, Replica, ReplicaState

__all__ = [
    "Config",
    "MDSettings",
    "Replica",
    "ReplicaState",
    "Project",
]
