"""
MD Settings
===========

Configuration for molecular dynamics simulations.

Replaces legacy MDSettings class with a modern dataclass-based approach.
Settings can come from:
1. Defaults (hardcoded)
2. Config files (TOML/CFG)
3. Constructor arguments (highest priority)

Examples
--------
>>> settings = MDSettings(solvent='ETA', nanos=40)
>>> settings.nanos
40
>>> settings.temperature
300.0
"""

from __future__ import annotations

from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Literal
import logging

# tomllib is Python 3.11+, use tomli for 3.10
try:
    import tomllib
except ImportError:
    import tomli as tomllib

log = logging.getLogger(__name__)


# Restraint mode type
RestraintMode = Literal["FREE", "BB", "HA", "CUSTOM"]

# MD program type
MDProgram = Literal["AMBER", "NAMD", "OPENMM", "GROMACS"]

# Ensemble type
Ensemble = Literal["NVT", "NPT"]


@dataclass
class MDSettings:
    """
    MD simulation parameters.
    
    Attributes
    ----------
    solvent : str
        Solvent mixture name (e.g., 'ETA', 'WAT', 'MAM')
    name : str | None
        Settings name. Auto-generated if not provided.
    nanos : int
        Production length in nanoseconds
    temperature : float
        Simulation temperature in Kelvin
    timestep : float
        Integration timestep in femtoseconds
    restraint_mode : RestraintMode
        Restraining scheme: FREE, BB (backbone), HA (heavy atoms), CUSTOM
    restraint_force : float
        Restraining force in kcal/mol·Å²
    restraint_mask : str
        Amber mask for restrained atoms (if CUSTOM mode)
    align_mask : str
        Amber mask for trajectory alignment
    md_program : MDProgram
        Simulation program to use
    production_ensemble : Ensemble
        Ensemble for production (NVT or NPT)
    trajectory_frequency : int
        Trajectory writing frequency (steps)
    minimization_steps : int
        Number of minimization steps
    heating_steps : int
        Heating equilibration steps
    npt_equilibration_steps : int
        NPT equilibration steps
    production_steps : int
        Steps per production file
    heating_initial_temp : float
        Starting temperature for heating
    npt_pressure : float
        Target pressure for NPT ensemble (bar)
    use_netcdf : bool
        Write trajectory in NetCDF format
    wrap_coordinates : bool
        Apply coordinate wrapping (IWRAP)
    """
    
    # Required
    solvent: str
    
    # Optional with defaults
    name: str | None = None
    nanos: int = 20
    temperature: float = 300.0
    timestep: float = 2.0
    restraint_mode: RestraintMode = "FREE"
    restraint_force: float = 0.0
    restraint_mask: str = ""
    align_mask: str = ""
    
    # Program settings
    md_program: MDProgram = "AMBER"
    production_ensemble: Ensemble = "NVT"
    
    # Step counts
    trajectory_frequency: int = 500
    minimization_steps: int = 5000
    heating_steps: int = 100000
    npt_equilibration_steps: int = 500000
    production_steps: int = 500000
    
    # Other parameters
    heating_initial_temp: float = 100.0
    npt_pressure: float = 1.0
    use_netcdf: bool = True
    wrap_coordinates: bool = True
    
    # Folder names
    md_folder: str = "md"
    eq_folder: str = "eq"
    min_folder: str = "min"
    align_folder: str = "align"
    density_folder: str = "dgrids"
    energy_folder: str = "egrids"
    
    # File templates
    md_file_template: str = "md{step}.{extension}"
    eq_file_template: str = "eq{step}.{extension}"
    
    # Force fields
    force_fields: list[str] = field(
        default_factory=lambda: ["leaprc.protein.ff14SB", "leaprc.gaff"]
    )
    
    def __post_init__(self):
        """Compute derived attributes."""
        # Auto-generate name if not provided
        if self.name is None:
            self.name = f"Settings_{self.solvent}_{self.restraint_mode}"
        
        # Validate restraint mode
        if self.restraint_mode not in ("FREE", "BB", "HA", "CUSTOM"):
            raise ValueError(f"Invalid restraint_mode: {self.restraint_mode}")
    
    @property
    def has_restraints(self) -> bool:
        """Check if restraints are enabled."""
        return self.restraint_mode != "FREE"
    
    @property
    def n_trajectory_files(self) -> int:
        """Calculate expected number of trajectory files."""
        # Convert ns to fs, divide by timestep to get total steps
        # Then divide by steps per file
        total_steps = (self.nanos * 1e6) / self.timestep
        return int(total_steps / self.production_steps)
    
    @property
    def n_snapshots(self) -> int:
        """Calculate total number of snapshots."""
        return int(self.nanos * (self.production_steps / self.trajectory_frequency))
    
    @property
    def snapshots_per_ns(self) -> int:
        """Snapshots per nanosecond."""
        return int(self.production_steps / self.trajectory_frequency)
    
    @property
    def trajectory_extension(self) -> str:
        """Get trajectory file extension based on format."""
        return "nc" if self.use_netcdf else "x"
    
    def to_dict(self) -> dict:
        """Convert to dictionary."""
        return asdict(self)
    
    @classmethod
    def from_dict(cls, data: dict) -> "MDSettings":
        """Create from dictionary."""
        # Filter to only known fields
        known = {f.name for f in cls.__dataclass_fields__.values()}
        filtered = {k: v for k, v in data.items() if k in known}
        return cls(**filtered)
    
    @classmethod
    def from_toml(cls, path: str | Path, section: str = "mdsettings") -> "MDSettings":
        """
        Load settings from TOML file.
        
        Parameters
        ----------
        path : str or Path
            Path to TOML file
        section : str
            Section name in TOML file
            
        Returns
        -------
        MDSettings
            Settings loaded from file
        """
        path = Path(path)
        if not path.exists():
            raise FileNotFoundError(f"Settings file not found: {path}")
        
        with open(path, "rb") as f:
            data = tomllib.load(f)
        
        if section not in data:
            raise KeyError(f"Section '{section}' not found in {path}")
        
        return cls.from_dict(data[section])
    
    def summary(self) -> str:
        """Get human-readable summary."""
        lines = [
            f"MD Settings: {self.name}",
            "-" * 40,
            f"Solvent:      {self.solvent}",
            f"Duration:     {self.nanos} ns",
            f"Temperature:  {self.temperature} K",
            f"Timestep:     {self.timestep} fs",
            f"Program:      {self.md_program}",
            f"Ensemble:     {self.production_ensemble}",
            f"Restraints:   {self.restraint_mode}",
        ]
        
        if self.has_restraints:
            lines.append(f"Restr. Force: {self.restraint_force} kcal/mol·Å²")
            if self.restraint_mask:
                lines.append(f"Restr. Mask:  {self.restraint_mask}")
        
        lines.extend([
            "-" * 40,
            f"Traj. files:  {self.n_trajectory_files}",
            f"Snapshots:    {self.n_snapshots}",
        ])
        
        return "\n".join(lines)
    
    def __str__(self) -> str:
        return self.summary()
    
    def __repr__(self) -> str:
        return f"MDSettings(solvent={self.solvent!r}, name={self.name!r})"


def load_settings(
    path: str | Path,
    solvent: str | None = None,
    **overrides,
) -> MDSettings:
    """
    Load settings from file with optional overrides.
    
    Parameters
    ----------
    path : str or Path
        Path to settings file (TOML)
    solvent : str | None
        Solvent name (required if not in file)
    **overrides
        Additional parameters to override
        
    Returns
    -------
    MDSettings
        Loaded settings
    """
    settings_dict = {}
    
    path = Path(path)
    if path.exists():
        with open(path, "rb") as f:
            data = tomllib.load(f)
        if "mdsettings" in data:
            settings_dict.update(data["mdsettings"])
    
    # Apply overrides
    settings_dict.update(overrides)
    
    # Solvent is required
    if solvent:
        settings_dict["solvent"] = solvent
    
    if "solvent" not in settings_dict:
        raise ValueError("Solvent name is required")
    
    return MDSettings.from_dict(settings_dict)
