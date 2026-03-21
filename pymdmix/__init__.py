"""
pyMDMix - Molecular Dynamics with Solvent Mixtures
===================================================

A toolkit for setting up, running, and analyzing molecular dynamics
simulations with organic solvent mixtures for drug discovery.

Example
-------
>>> from pymdmix import open_trajectory, Grid
>>> traj = open_trajectory("system.prmtop", "trajectory.nc")
>>> for frame in traj:
...     # Process frame coordinates
...     pass
"""

__version__ = "0.3.0"
__author__ = "Daniel Alvarez-Garcia"

from pymdmix.core.grid import Grid
from pymdmix.core.trajectory import open_trajectory, Frame
from pymdmix.core.structure import load_structure

__all__ = [
    "Grid",
    "open_trajectory",
    "Frame",
    "load_structure",
    "__version__",
]
