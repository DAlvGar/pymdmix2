"""
Configuration Management
========================

Global and per-simulation configuration for MDMix.

Replaces the old settings.py, SettingsParser.py, MDSettings.py, and Parsers.py
with a cleaner, type-safe approach using dataclasses and JSON/YAML.

Examples
--------
>>> from pymdmix.project import Config, MDSettings
>>>
>>> # Load config
>>> config = Config.from_file("mdmix.yaml")
>>>
>>> # Access settings
>>> print(config.amber_home)
>>> print(config.md_settings.nsteps)
"""

from __future__ import annotations

import json
import logging
import os
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

log = logging.getLogger(__name__)

# Try to import YAML support
try:
    import yaml

    HAS_YAML = True
except ImportError:
    HAS_YAML = False
    yaml = None


@dataclass
class MDSettings:
    """
    Molecular dynamics simulation settings.

    Attributes
    ----------
    nsteps : int
        Number of MD steps per production run
    timestep : float
        Integration timestep in femtoseconds
    temperature : float
        Target temperature in Kelvin
    pressure : float
        Target pressure in bar (for NPT)
    ntwx : int
        Trajectory write frequency (steps)
    ntpr : int
        Energy output frequency (steps)
    ntwr : int
        Restart file frequency (steps)
    ntr : int
        Restraint flag (0=off, 1=on)
    restraint_wt : float
        Restraint weight in kcal/mol/Å²
    restraint_mask : str
        Amber mask for restrained atoms
    """

    # Production settings
    nsteps: int = 500000
    timestep: float = 2.0
    temperature: float = 300.0
    pressure: float = 1.0

    # Output frequencies
    ntwx: int = 500
    ntpr: int = 500
    ntwr: int = 5000

    # Restraints
    ntr: int = 0
    restraint_wt: float = 5.0
    restraint_mask: str = "@CA"

    # Cutoffs
    cut: float = 10.0

    # Thermostat/Barostat
    ntt: int = 3  # Langevin thermostat
    gamma_ln: float = 2.0
    ntp: int = 1  # Isotropic pressure
    barostat: int = 2  # Monte Carlo barostat

    # Shake
    ntc: int = 2  # Shake on hydrogens
    ntf: int = 2  # No force calc for shaken bonds

    def to_dict(self) -> dict[str, Any]:
        """Convert to dictionary."""
        return asdict(self)

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> MDSettings:
        """Create from dictionary."""
        # Filter to only known fields
        known_fields = {f.name for f in cls.__dataclass_fields__.values()}
        filtered = {k: v for k, v in data.items() if k in known_fields}
        return cls(**filtered)

    def to_amber_mdin(self) -> str:
        """Generate Amber mdin content."""
        return f"""Production MD
&cntrl
  imin=0, irest=1, ntx=5,
  nstlim={self.nsteps}, dt={self.timestep / 1000:.4f},
  temp0={self.temperature}, tempi={self.temperature},
  ntt={self.ntt}, gamma_ln={self.gamma_ln},
  ntp={self.ntp}, barostat={self.barostat}, pres0={self.pressure},
  ntc={self.ntc}, ntf={self.ntf}, cut={self.cut},
  ntwx={self.ntwx}, ntpr={self.ntpr}, ntwr={self.ntwr},
  ntxo=2, ioutfm=1,
  ntr={self.ntr},
&end
"""


@dataclass
class Config:
    """
    Global MDMix configuration.

    Handles paths, environment variables, and default settings.

    Attributes
    ----------
    amber_home : Path | None
        AMBER installation directory
    work_dir : Path
        Working directory for projects
    md_settings : MDSettings
        Default MD simulation settings
    default_solvent : str
        Default solvent for new systems
    n_replicas : int
        Default number of replicas
    queue_system : str
        Default queue system (slurm, sge, pbs, local)
    """

    amber_home: Path | None = None
    work_dir: Path = field(default_factory=Path.cwd)
    md_settings: MDSettings = field(default_factory=MDSettings)
    default_solvent: str = "ETA"
    n_replicas: int = 3
    queue_system: str = "slurm"

    # Executables (relative to amber_home/bin or absolute)
    tleap_exe: str = "tleap"
    cpptraj_exe: str = "cpptraj"
    pmemd_exe: str = "pmemd.cuda"

    # Force fields
    protein_ff: str = "leaprc.protein.ff19SB"
    water_ff: str = "leaprc.water.opc"

    def __post_init__(self):
        """Resolve paths and environment variables."""
        # Get AMBERHOME from environment if not set
        if self.amber_home is None:
            env_amber = os.environ.get("AMBERHOME")
            if env_amber:
                self.amber_home = Path(env_amber)

        # Convert work_dir to Path if string
        if isinstance(self.work_dir, str):
            self.work_dir = Path(self.work_dir)

        # Convert md_settings from dict if needed
        if isinstance(self.md_settings, dict):
            self.md_settings = MDSettings.from_dict(self.md_settings)

    @property
    def amber_bin(self) -> Path | None:
        """Path to AMBER bin directory."""
        if self.amber_home:
            bin_path = self.amber_home / "bin"
            if bin_path.exists():
                return bin_path
        return None

    def get_executable(self, name: str) -> Path | None:
        """Get path to an executable."""
        exe_name = getattr(self, f"{name}_exe", name)

        # Check if absolute path
        if Path(exe_name).is_absolute():
            return Path(exe_name) if Path(exe_name).exists() else None

        # Check in AMBER bin
        if self.amber_bin:
            exe_path = self.amber_bin / exe_name
            if exe_path.exists():
                return exe_path

        # Check in PATH (would need shutil.which)
        return None

    def validate(self) -> list[str]:
        """
        Validate configuration.

        Returns
        -------
        list[str]
            List of validation errors
        """
        errors = []

        if self.amber_home is None:
            errors.append("AMBERHOME not set")
        elif not self.amber_home.exists():
            errors.append(f"AMBERHOME does not exist: {self.amber_home}")

        if self.amber_bin and not self.amber_bin.exists():
            errors.append(f"AMBER bin directory not found: {self.amber_bin}")

        return errors

    def to_dict(self) -> dict[str, Any]:
        """Convert to dictionary for serialization."""
        return {
            "amber_home": str(self.amber_home) if self.amber_home else None,
            "work_dir": str(self.work_dir),
            "md_settings": self.md_settings.to_dict(),
            "default_solvent": self.default_solvent,
            "n_replicas": self.n_replicas,
            "queue_system": self.queue_system,
            "tleap_exe": self.tleap_exe,
            "cpptraj_exe": self.cpptraj_exe,
            "pmemd_exe": self.pmemd_exe,
            "protein_ff": self.protein_ff,
            "water_ff": self.water_ff,
        }

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> Config:
        """Create Config from dictionary."""
        # Handle nested MDSettings
        md_data = data.pop("md_settings", {})
        md_settings = MDSettings.from_dict(md_data) if md_data else MDSettings()

        # Handle Path fields
        if "amber_home" in data and data["amber_home"]:
            data["amber_home"] = Path(data["amber_home"])
        if "work_dir" in data:
            data["work_dir"] = Path(data["work_dir"])

        return cls(md_settings=md_settings, **data)

    def to_json(self, path: str | Path) -> None:
        """Save configuration to JSON file."""
        with open(path, "w") as f:
            json.dump(self.to_dict(), f, indent=2)
        log.info(f"Saved config to {path}")

    @classmethod
    def from_json(cls, path: str | Path) -> Config:
        """Load configuration from JSON file."""
        with open(path) as f:
            data = json.load(f)
        return cls.from_dict(data)

    def to_yaml(self, path: str | Path) -> None:
        """Save configuration to YAML file."""
        if not HAS_YAML:
            raise ImportError("PyYAML required for YAML support")

        with open(path, "w") as f:
            yaml.dump(self.to_dict(), f, default_flow_style=False)
        log.info(f"Saved config to {path}")

    @classmethod
    def from_yaml(cls, path: str | Path) -> Config:
        """Load configuration from YAML file."""
        if not HAS_YAML:
            raise ImportError("PyYAML required for YAML support")

        with open(path) as f:
            data = yaml.safe_load(f)
        return cls.from_dict(data)

    @classmethod
    def from_file(cls, path: str | Path) -> Config:
        """Load configuration from file (auto-detect format)."""
        path = Path(path)

        if path.suffix in (".yaml", ".yml"):
            return cls.from_yaml(path)
        elif path.suffix == ".json":
            return cls.from_json(path)
        else:
            raise ValueError(f"Unknown config format: {path.suffix}")

    def save(self, path: str | Path) -> None:
        """Save configuration to file (auto-detect format)."""
        path = Path(path)

        if path.suffix in (".yaml", ".yml"):
            self.to_yaml(path)
        else:
            self.to_json(path)


def get_default_config() -> Config:
    """
    Get default configuration, reading from environment.

    Returns
    -------
    Config
        Configuration with defaults and environment variables
    """
    return Config()
