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

from pymdmix.analysis.base import (
    Action,
    ActionResult,
    register_action,
    get_action,
    list_actions,
    run_action,
)
from pymdmix.analysis.density import DensityAction, calculate_density
from pymdmix.analysis.residence import ResidenceAction, calculate_residence
from pymdmix.analysis.hotspots import HotspotAction, Hotspot, detect_hotspots
from pymdmix.analysis.align import align_trajectory, align_replica, AlignmentResult
from pymdmix.analysis.energy import (
    density_to_free_energy,
    boltzmann_average,
    calculate_expected_density,
    normalize_grid,
    replica_average,
    EnergyResult,
)

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
    "align_trajectory",
    "align_replica",
    "AlignmentResult",
    # Energy
    "density_to_free_energy",
    "boltzmann_average",
    "calculate_expected_density",
    "normalize_grid",
    "replica_average",
    "EnergyResult",
]
