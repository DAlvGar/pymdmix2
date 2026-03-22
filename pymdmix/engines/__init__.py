"""
MD Engines
==========

Interfaces to molecular dynamics engines for running simulations.

Supported engines:
- Amber (pmemd, sander) - primary
- OpenMM (optional) - GPU-accelerated alternative
- NAMD (optional) - parallel MD engine
- GROMACS - high-performance MD

Also includes queue script generation for HPC clusters.
"""

from pymdmix.engines.amber import AmberEngine, AmberInput
from pymdmix.engines.executor import (
    AsyncExecutor,
    Executor,
    Job,
    JobResult,
    run_command,
)
from pymdmix.engines.gromacs import GromacsCheck, GromacsConfig, GromacsEngine
from pymdmix.engines.namd import NAMDConfig, NAMDEngine
from pymdmix.engines.openmm import (
    OpenMMConfig,
    OpenMMEngine,
    apply_harmonic_positional_restraints,
    get_backbone_indices,
    get_heavy_atom_indices,
    set_context_from_restart,
)
from pymdmix.engines.queue import QueueConfig, generate_queue_script

__all__ = [
    # Amber
    "AmberEngine",
    "AmberInput",
    # OpenMM
    "OpenMMEngine",
    "OpenMMConfig",
    # NAMD
    "NAMDEngine",
    "NAMDConfig",
    # GROMACS
    "GromacsEngine",
    "GromacsConfig",
    "GromacsCheck",
    # Queue
    "QueueConfig",
    "generate_queue_script",
    # Executor
    "Executor",
    "AsyncExecutor",
    "Job",
    "JobResult",
    "run_command",
    # OpenMM utilities
    "apply_harmonic_positional_restraints",
    "set_context_from_restart",
    "get_backbone_indices",
    "get_heavy_atom_indices",
]
