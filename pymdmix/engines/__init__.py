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
from pymdmix.engines.queue import QueueConfig, generate_queue_script
from pymdmix.engines.openmm import OpenMMEngine, OpenMMConfig
from pymdmix.engines.namd import NAMDEngine, NAMDConfig
from pymdmix.engines.gromacs import GromacsEngine, GromacsConfig, GromacsCheck

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
]
