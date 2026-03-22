"""
Analysis Actions
================

Plugin-based analysis framework for MD trajectory analysis.

Available actions:
- density: Calculate probe density grids
- residence: Calculate residence times at hotspots
- hotspots: Detect binding hotspots from density grids
- align: Trajectory alignment
- energy: Free energy conversion and Boltzmann averaging

Examples
--------
>>> from pymdmix.analysis import DensityAction, run_action
>>> action = DensityAction()
>>> result = action.run(replica, trajectory, probes=["OH", "CT"])

>>> from pymdmix.analysis import boltzmann_average
>>> avg_grid = boltzmann_average(["rep1.dx", "rep2.dx"])
"""

from pymdmix.analysis.align import AlignAction, AlignmentResult, align_replica, align_trajectory
from pymdmix.analysis.base import (
    Action,
    ActionResult,
    get_action,
    list_actions,
    register_action,
    run_action,
)
from pymdmix.analysis.density import DensityAction, calculate_density
from pymdmix.analysis.energy import (
    EnergyAction,
    EnergyResult,
    boltzmann_average,
    calculate_expected_density,
    density_to_free_energy,
    normalize_grid,
    replica_average,
)
from pymdmix.analysis.hotspots import Hotspot, HotspotAction, detect_hotspots
from pymdmix.analysis.residence import ResidenceAction, calculate_residence

__all__ = [
    # Base
    "Action",
    "ActionResult",
    "register_action",
    "get_action",
    "list_actions",
    "run_action",
    # Density
    "DensityAction",
    "calculate_density",
    # Residence
    "ResidenceAction",
    "calculate_residence",
    # Hotspots
    "HotspotAction",
    "Hotspot",
    "detect_hotspots",
    # Alignment
    "AlignAction",
    "align_trajectory",
    "align_replica",
    "AlignmentResult",
    # Energy
    "EnergyAction",
    "density_to_free_energy",
    "boltzmann_average",
    "calculate_expected_density",
    "normalize_grid",
    "replica_average",
    "EnergyResult",
]
