"""
Amber MD Engine
===============

Interface for running Amber molecular dynamics simulations.

Includes:
- LeapSession: Interactive tleap subprocess control
- LeapBuilder: System building (solvation, neutralization)
- AmberChecker: Verify simulation completion
- AmberWriter: Generate input files
- AmberEngine: MD input generation

Examples
--------
>>> from pymdmix.engines.amber import LeapSession, AmberEngine
>>>
>>> # Use tleap interactively
>>> with LeapSession() as leap:
...     leap.command("source leaprc.protein.ff19SB")
...     leap.load_off("system.off")
...     leap.save_amber_parm("mol", "system.prmtop", "system.rst7")
>>>
>>> # Generate MD input
>>> engine = AmberEngine(exe="pmemd.cuda")
>>> mdin = engine.production_input(nsteps=500000)
"""

from __future__ import annotations

import fnmatch
import logging
import os
import re
import shutil
import subprocess
import tempfile
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from pymdmix.project.replica import Replica

log = logging.getLogger(__name__)


# =============================================================================
# Exceptions
# =============================================================================

class AmberError(Exception):
    """Base exception for Amber operations."""
    pass


class LeapError(AmberError):
    """Error in tleap execution."""
    pass


class AmberCheckError(AmberError):
    """Error checking simulation status."""
    pass


class AmberWriterError(AmberError):
    """Error writing input files."""
    pass


# =============================================================================
# Environment Detection
# =============================================================================

def get_amber_home() -> Path | None:
    """Get AMBERHOME from environment."""
    amber_home = os.environ.get("AMBERHOME")
    if amber_home:
        return Path(amber_home)
    return None


def get_amber_exe() -> Path | None:
    """Get Amber executables directory."""
    amber_home = get_amber_home()
    if not amber_home:
        return None
        
    # Check common locations
    for subdir in ("exe", "bin"):
        exe_path = amber_home / subdir
        if exe_path.exists():
            return exe_path
            
    return None


def find_forcefield(name: str) -> Path | None:
    """
    Find forcefield file by name.
    
    Searches in AMBERHOME/dat/leap/cmd and AMBERHOME/dat/leap/parm.
    
    Parameters
    ----------
    name : str
        Forcefield name (e.g., "leaprc.protein.ff19SB")
        
    Returns
    -------
    Path | None
        Path to forcefield file, or None if not found.
    """
    # Check if it's already a valid path
    if Path(name).exists():
        return Path(name).resolve()
        
    amber_home = get_amber_home()
    if not amber_home:
        return None
        
    search_paths = [
        amber_home / "dat" / "leap" / "cmd",
        amber_home / "dat" / "leap" / "parm",
    ]
    
    for search_path in search_paths:
        if not search_path.exists():
            continue
            
        for root, _, files in os.walk(search_path):
            matches = fnmatch.filter(files, f"*{name}*")
            if matches:
                if name in matches:
                    return Path(root) / name
                return Path(root) / matches[0]
                
    return None


# =============================================================================
# LeapSession - Interactive tleap Control
# =============================================================================

class LeapSession:
    """
    Interactive tleap subprocess session.
    
    Provides a context manager for running tleap commands and
    capturing output.
    
    Examples
    --------
    >>> with LeapSession() as leap:
    ...     leap.source("leaprc.protein.ff19SB")
    ...     leap.load_off("system.off")
    ...     output = leap.command("desc mol")
    ...     leap.save_amber_parm("mol", "out.prmtop", "out.rst7")
    """
    
    def __init__(
        self,
        extra_ff: list[str] | None = None,
        cwd: Path | str | None = None,
    ):
        """
        Initialize leap session.
        
        Parameters
        ----------
        extra_ff : list[str] | None
            Additional forcefield files to load.
        cwd : Path | str | None
            Working directory for tleap.
        """
        self.extra_ff = extra_ff or []
        self.cwd = Path(cwd) if cwd else None
        self._process: subprocess.Popen | None = None
        self._terminator = 'mdmix_terminator_command = "mdmix_terminator_command"'
        self._terminator_check = "desc mdmix_terminator_command"
        
    def __enter__(self) -> LeapSession:
        self.start()
        return self
        
    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        self.close()
        
    def start(self) -> None:
        """Start tleap subprocess."""
        exe_path = get_amber_exe()
        if not exe_path:
            raise LeapError("AMBERHOME not set or tleap not found")
            
        tleap = exe_path / "tleap"
        if not tleap.exists():
            raise LeapError(f"tleap not found at {tleap}")
            
        self._process = subprocess.Popen(
            [str(tleap), "-f", "-"],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            cwd=self.cwd,
            text=True,
        )
        
        # Set up terminator for output detection
        self._write(self._terminator)
        
        # Load extra forcefields
        for ff in self.extra_ff:
            if "frcmod" in ff:
                self.command(f"loadAmberParams {ff}")
            else:
                self.command(f"source {ff}")
                
        log.debug("LeapSession started")
        
    def close(self) -> None:
        """Close tleap subprocess."""
        if self._process:
            try:
                self._process.terminate()
                self._process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                self._process.kill()
            self._process = None
            
        log.debug("LeapSession closed")
        
    def _write(self, cmd: str) -> None:
        """Write command to tleap stdin."""
        if not self._process or not self._process.stdin:
            raise LeapError("LeapSession not started")
        self._process.stdin.write(cmd + "\n")
        self._process.stdin.flush()
        
    def _read_output(self) -> list[str]:
        """Read output until terminator."""
        if not self._process or not self._process.stdout:
            return []
            
        output = []
        while True:
            line = self._process.stdout.readline().strip()
            
            if "mdmix_terminator_command" in line:
                break
            if self._process.poll() is not None:
                break
            if "Fatal Error!" in line:
                break
            if "Exiting LEaP" in line:
                break
                
            if line:
                output.append(line)
                log.debug(f"tleap: {line}")
                
        return output
        
    def command(self, cmd: str) -> list[str]:
        """
        Execute a tleap command.
        
        Parameters
        ----------
        cmd : str
            tleap command string.
            
        Returns
        -------
        list[str]
            Output lines from command.
        """
        log.debug(f"tleap command: {cmd}")
        self._write(cmd)
        self._write(self._terminator_check)
        return self._read_output()
        
    def source(self, ff_name: str) -> list[str]:
        """Source a forcefield file."""
        ff_path = find_forcefield(ff_name)
        if ff_path:
            return self.command(f"source {ff_path}")
        return self.command(f"source {ff_name}")
        
    def load_off(self, off_file: Path | str) -> list[str]:
        """Load an OFF file."""
        return self.command(f"loadOff {off_file}")
        
    def load_pdb(self, pdb_file: Path | str, unit_name: str = "mol") -> list[str]:
        """Load a PDB file."""
        return self.command(f"{unit_name} = loadPdb {pdb_file}")
        
    def save_amber_parm(
        self,
        unit: str,
        top: Path | str,
        crd: Path | str,
    ) -> bool:
        """
        Save Amber topology and coordinates.
        
        Parameters
        ----------
        unit : str
            Unit name in tleap.
        top : Path | str
            Output topology file.
        crd : Path | str
            Output coordinate file.
            
        Returns
        -------
        bool
            True if files were created successfully.
        """
        output = self.command(f"saveAmberParm {unit} {top} {crd}")
        time.sleep(0.5)  # Give filesystem time to sync
        
        top_path = Path(top)
        crd_path = Path(crd)
        
        if top_path.exists() and crd_path.exists():
            if top_path.stat().st_size > 0 and crd_path.stat().st_size > 0:
                log.info(f"Saved topology and coordinates for {unit}")
                return True
                
        log.error(f"Failed to save {top} / {crd}")
        log.error(f"tleap output: {output}")
        return False
        
    def charge(self, unit: str) -> float:
        """Get charge of a unit."""
        output = self.command(f"charge {unit}")
        if output:
            try:
                return float(output[0].split()[-1])
            except (ValueError, IndexError):
                pass
        return 0.0
        
    def solvate_box(
        self,
        unit: str,
        solvent_box: str,
        buffer: float = 12.0,
        cubic: bool = False,
    ) -> list[str]:
        """
        Solvate a unit with a solvent box.
        
        Parameters
        ----------
        unit : str
            Unit to solvate.
        solvent_box : str
            Solvent box name.
        buffer : float
            Buffer distance in Angstroms.
        cubic : bool
            Use cubic (True) or octahedral (False) box.
        """
        cmd = "solvateBox" if cubic else "solvateOct"
        return self.command(f"{cmd} {unit} {solvent_box} {buffer:.2f} iso 1")
        
    def add_ions(
        self,
        unit: str,
        ion: str = "Na+",
        count: int = 0,
    ) -> list[str]:
        """Add ions to neutralize or reach specific count."""
        return self.command(f"addIons {unit} {ion} {count}")
        
    def neutralize_nacl(self, unit: str) -> bool:
        """Neutralize with Na+ and Cl-."""
        self.add_ions(unit, "Na+", 0)
        self.add_ions(unit, "Cl-", 0)
        
        charge = self.charge(unit)
        if abs(charge) < 0.01:
            log.info(f"Neutralized {unit} with NaCl")
            return True
            
        log.warning(f"Could not fully neutralize {unit}, charge: {charge:.3f}")
        return False


# =============================================================================
# LeapBuilder - System Building
# =============================================================================

class LeapBuilder:
    """
    Build solvated systems using tleap.
    
    Handles forcefield loading, solvation, and neutralization.
    
    Examples
    --------
    >>> builder = LeapBuilder()
    >>> builder.add_ff("leaprc.protein.ff19SB")
    >>> builder.load_off("protein.off")
    >>> builder.solvate("mol", "ETA", buffer=12.0)
    >>> builder.save("mol", "solvated.prmtop", "solvated.rst7")
    """
    
    def __init__(
        self,
        replica: Replica | None = None,
        forcefields: list[str] | None = None,
    ):
        """
        Initialize builder.
        
        Parameters
        ----------
        replica : Replica | None
            Replica to build system for.
        forcefields : list[str] | None
            Initial forcefield list.
        """
        self.replica = replica
        self.forcefields: list[str] = []
        self.leap: LeapSession | None = None
        
        # Load default forcefields
        defaults = ["leaprc.protein.ff19SB", "leaprc.gaff2", "leaprc.water.tip3p"]
        for ff in defaults:
            self.add_ff(ff)
            
        if forcefields:
            for ff in forcefields:
                self.add_ff(ff)
                
    def add_ff(self, ff_name: str) -> bool:
        """
        Add a forcefield.
        
        Parameters
        ----------
        ff_name : str
            Forcefield name or path.
            
        Returns
        -------
        bool
            True if forcefield was found.
        """
        ff_path = find_forcefield(ff_name)
        if ff_path:
            if str(ff_path) not in self.forcefields:
                self.forcefields.append(str(ff_path))
                log.debug(f"Added forcefield: {ff_path}")
            return True
            
        # Try adding as-is
        if ff_name not in self.forcefields:
            self.forcefields.append(ff_name)
            log.warning(f"Forcefield not found, using as-is: {ff_name}")
        return False
        
    def init_leap(self) -> LeapSession:
        """Initialize tleap session."""
        if self.leap is None:
            self.leap = LeapSession(extra_ff=self.forcefields)
            self.leap.start()
        return self.leap
        
    def close(self) -> None:
        """Close leap session."""
        if self.leap:
            self.leap.close()
            self.leap = None
            
    def load_off(self, off_file: Path | str) -> bool:
        """Load an OFF file."""
        leap = self.init_leap()
        leap.load_off(off_file)
        return True
        
    def solvate(
        self,
        unit: str,
        solvent: str,
        buffer: float = 12.0,
        cubic: bool = False,
    ) -> bool:
        """
        Solvate a unit with organic solvent mixture.
        
        Parameters
        ----------
        unit : str
            Unit name to solvate.
        solvent : str
            Solvent name (e.g., "ETA", "MAM").
        buffer : float
            Box buffer distance.
        cubic : bool
            Use cubic box.
            
        Returns
        -------
        bool
            True if solvation succeeded.
        """
        from pymdmix.core.solvent import SolventLibrary
        
        leap = self.init_leap()
        library = SolventLibrary()
        
        try:
            solv = library.get(solvent)
        except KeyError:
            log.error(f"Unknown solvent: {solvent}")
            return False
            
        # Write solvent OFF to temp file and load
        with tempfile.NamedTemporaryFile(suffix=".off", delete=False) as tmp:
            tmp_path = Path(tmp.name)
            
        try:
            solv.write_off(tmp_path)
            leap.load_off(tmp_path)
            
            # Solvate
            box_unit = solv.box_unit if hasattr(solv, "box_unit") else f"{solvent}BOX"
            leap.solvate_box(unit, box_unit, buffer, cubic)
            
            # Neutralize
            if hasattr(solv, "is_ionic") and solv.is_ionic():
                # Handle ionic solvents specially
                self._neutralize_ionic(unit, solv)
            else:
                leap.neutralize_nacl(unit)
                
            return True
            
        finally:
            tmp_path.unlink(missing_ok=True)
            
    def _neutralize_ionic(self, unit: str, solvent: Any) -> bool:
        """Neutralize ionic solvent box."""
        if not self.leap:
            return False
            
        charge = self.leap.charge(unit)
        if abs(charge) < 0.01:
            return True
            
        # Try standard NaCl neutralization
        return self.leap.neutralize_nacl(unit)
        
    def save(
        self,
        unit: str,
        top: Path | str,
        crd: Path | str,
    ) -> bool:
        """Save topology and coordinates."""
        if not self.leap:
            self.init_leap()
        return self.leap.save_amber_parm(unit, top, crd)
        
    def create_off(
        self,
        pdb: Path | str,
        output: Path | str,
        unit_name: str = "mol",
        extra_ff: list[str] | None = None,
    ) -> bool:
        """
        Create OFF file from PDB.
        
        Parameters
        ----------
        pdb : Path | str
            Input PDB file.
        output : Path | str
            Output OFF file.
        unit_name : str
            Unit name.
        extra_ff : list[str] | None
            Additional forcefields.
            
        Returns
        -------
        bool
            True if successful.
        """
        leap = self.init_leap()
        
        if extra_ff:
            for ff in extra_ff:
                if "frcmod" in ff:
                    leap.command(f"loadAmberParams {ff}")
                else:
                    leap.source(ff)
                    
        leap.load_pdb(pdb, unit_name)
        leap.command(f"saveOff {unit_name} {output}")
        
        return Path(output).exists()


# =============================================================================
# AmberChecker - Verify Simulation Status
# =============================================================================

class AmberChecker:
    """
    Check Amber simulation completion status.
    
    Verifies minimization, equilibration, and production stages.
    
    Examples
    --------
    >>> checker = AmberChecker(replica)
    >>> if checker.check_production([1, 2, 3]):
    ...     print("First 3 ns complete")
    """
    
    def __init__(self, replica: Replica, warn: bool = True):
        """
        Initialize checker.
        
        Parameters
        ----------
        replica : Replica
            Replica to check.
        warn : bool
            Log warnings for incomplete stages.
        """
        self.replica = replica
        self.warn = warn
        
    def check_minimization(self) -> bool:
        """Check if minimization is complete."""
        min_path = self.replica.min_path
        if not min_path or not min_path.exists():
            if self.warn:
                log.warning(f"Minimization folder not found: {min_path}")
            return False
            
        # Check for restart file
        rst_files = list(min_path.glob("*.rst7")) + list(min_path.glob("*.rst"))
        if rst_files:
            return True
            
        # Check output file for completion
        out_files = list(min_path.glob("*.out"))
        for out_file in out_files:
            content = out_file.read_text()
            if "FINAL RESULTS" in content or "Final Performance" in content:
                return True
                
        if self.warn:
            log.warning("Minimization not complete")
        return False
        
    def check_equilibration(
        self,
        step_selection: list[int] | None = None,
    ) -> bool:
        """
        Check if equilibration is complete.
        
        Parameters
        ----------
        step_selection : list[int] | None
            Specific equilibration steps to check.
        """
        eq_path = self.replica.eq_path
        if not eq_path or not eq_path.exists():
            if self.warn:
                log.warning(f"Equilibration folder not found: {eq_path}")
            return False
            
        if step_selection is None:
            step_selection = [1, 2]  # Default: eq1, eq2
            
        for step in step_selection:
            # Look for output/restart files
            patterns = [f"eq{step}.out", f"eq{step:02d}.out", f"eq_{step}.out"]
            found = False
            
            for pattern in patterns:
                matches = list(eq_path.glob(pattern))
                if matches:
                    content = matches[0].read_text()
                    if "FINAL RESULTS" in content or "Final Performance" in content:
                        found = True
                        break
                        
            if not found:
                if self.warn:
                    log.warning(f"Equilibration step {step} not complete")
                return False
                
        return True
        
    def check_production(
        self,
        step_selection: list[int] | None = None,
    ) -> bool:
        """
        Check if production MD is complete.
        
        Parameters
        ----------
        step_selection : list[int] | None
            Specific production steps to check.
        """
        md_path = self.replica.md_path
        if not md_path or not md_path.exists():
            if self.warn:
                log.warning(f"Production folder not found: {md_path}")
            return False
            
        if step_selection is None:
            if self.replica.settings:
                step_selection = list(range(1, self.replica.settings.n_traj_files + 1))
            else:
                step_selection = [1]
                
        extensions = self.replica.check_production_extension(step_selection)
        
        for step in step_selection:
            if extensions.get(step) is None:
                if self.warn:
                    log.warning(f"Production step {step} trajectory not found")
                return False
                
        return True
        
    def get_sim_volume(self, step: int = 1) -> float | None:
        """
        Get simulation box volume from output file.
        
        Parameters
        ----------
        step : int
            Production step number.
            
        Returns
        -------
        float | None
            Volume in Å³, or None if not found.
        """
        md_path = self.replica.md_path
        if not md_path:
            return None
            
        # Find output file
        template = self.replica.md_output_template
        out_name = template.format(step=step, extension="out")
        out_file = md_path / out_name
        
        if not out_file.exists():
            return None
            
        # Parse volume from output
        volume_pattern = re.compile(r"VOLUME\s+=\s+([\d.]+)")
        
        content = out_file.read_text()
        matches = volume_pattern.findall(content)
        
        if matches:
            return float(matches[-1])  # Last value
            
        return None


# =============================================================================
# AmberWriter - Generate Input Files
# =============================================================================

class AmberWriter:
    """
    Generate Amber input files for a replica.
    
    Creates minimization, equilibration, and production inputs
    along with submission scripts.
    
    Examples
    --------
    >>> writer = AmberWriter(replica)
    >>> writer.write_replica_input()
    >>> writer.write_commands()
    """
    
    def __init__(self, replica: Replica):
        """
        Initialize writer.
        
        Parameters
        ----------
        replica : Replica
            Replica to write inputs for.
        """
        self.replica = replica
        self.engine = AmberEngine()
        
    def write_replica_input(self) -> None:
        """Write all MD input files."""
        if not self.replica.path:
            raise AmberWriterError("Replica path not set")
            
        settings = self.replica.settings
        if not settings:
            raise AmberWriterError("Replica settings not configured")
            
        # Get restraint mask
        restr_mask = self._get_restraint_mask()
        
        # Minimization
        min_path = self.replica.min_path
        if min_path:
            min_path.mkdir(exist_ok=True)
            
            mdin = self.engine.minimization_input(
                restraint_wt=settings.restraint_force if settings.has_restraints else 0,
                restraint_mask=restr_mask,
            )
            (min_path / "min.in").write_text(mdin)
            
        # Equilibration
        eq_path = self.replica.eq_path
        if eq_path:
            eq_path.mkdir(exist_ok=True)
            
            # Heating (NVT)
            heat = self.engine.heating_input(
                temp_end=settings.temperature,
                restraint_wt=settings.restraint_force if settings.has_restraints else 5.0,
                restraint_mask=restr_mask or "@CA",
            )
            (eq_path / "eq1.in").write_text(heat)
            
            # Equilibration (NPT)
            equil = self.engine.equilibration_input(
                temp=settings.temperature,
                restraint_wt=settings.restraint_force / 2 if settings.has_restraints else 1.0,
                restraint_mask=restr_mask or "@CA",
            )
            (eq_path / "eq2.in").write_text(equil)
            
        # Production
        md_path = self.replica.md_path
        if md_path:
            md_path.mkdir(exist_ok=True)
            
            prod = self.engine.production_input(
                nsteps=settings.prod_steps,
                dt=settings.timestep / 1000,  # fs to ps
                temp=settings.temperature,
                ntwx=settings.traj_frequency,
            )
            (md_path / "prod.in").write_text(prod)
            
        log.info(f"Wrote input files for replica {self.replica.name}")
        
    def write_commands(self, output: Path | str = "COMMANDS.sh") -> None:
        """
        Write shell script with all run commands.
        
        Parameters
        ----------
        output : Path | str
            Output script path.
        """
        if not self.replica.path:
            raise AmberWriterError("Replica path not set")
            
        settings = self.replica.settings
        if not settings:
            raise AmberWriterError("Replica settings not configured")
            
        lines = [
            "#!/bin/bash",
            f"# Commands for replica {self.replica.name}",
            "",
            "# Set paths",
            f"TOP={self.replica.topology}",
            f"CRD={self.replica.coordinates}",
            f"REF={self.replica.reference or self.replica.coordinates}",
            "",
        ]
        
        # Minimization
        lines.extend([
            "# Minimization",
            f"cd {self.replica.min_folder}",
            self.engine.run_command(
                Path("..") / self.replica.topology,
                Path("..") / self.replica.coordinates,
                Path("min.in"),
                "min",
            ),
            "cd ..",
            "",
        ])
        
        # Equilibration
        lines.extend([
            "# Equilibration",
            f"cd {self.replica.eq_folder}",
        ])
        
        # eq1 (heating)
        cmd = self.engine.run_command(
            Path("..") / self.replica.topology,
            Path("..") / self.replica.min_folder / "min.rst7",
            Path("eq1.in"),
            "eq1",
        )
        lines.append(cmd)
        
        # eq2 (equilibration)
        cmd = self.engine.run_command(
            Path("..") / self.replica.topology,
            Path("eq1.rst7"),
            Path("eq2.in"),
            "eq2",
        )
        lines.append(cmd)
        lines.extend(["cd ..", ""])
        
        # Production
        lines.extend([
            "# Production",
            f"cd {self.replica.md_folder}",
        ])
        
        for step in range(1, settings.n_traj_files + 1):
            if step == 1:
                prev_rst = Path("..") / self.replica.eq_folder / "eq2.rst7"
            else:
                prev_name = self.replica.md_output_template.format(step=step-1, extension="rst7")
                prev_rst = Path(prev_name)
                
            out_prefix = self.replica.md_output_template.format(step=step, extension="").rstrip(".")
            
            cmd = self.engine.run_command(
                Path("..") / self.replica.topology,
                prev_rst,
                Path("prod.in"),
                out_prefix,
            )
            lines.append(f"# Step {step}")
            lines.append(cmd)
            
        lines.extend(["cd ..", ""])
        
        script = "\n".join(lines)
        script_path = self.replica.path / output
        script_path.write_text(script)
        script_path.chmod(0o755)
        
        log.info(f"Wrote commands script: {script_path}")
        
    def _get_restraint_mask(self) -> str:
        """Get Amber restraint mask for replica."""
        settings = self.replica.settings
        if not settings:
            return ""
            
        if settings.restraint_mask:
            return settings.restraint_mask
            
        mode = settings.restraint_mode.upper()
        
        if mode == "FREE":
            return ""
        elif mode == "BB":
            return "@CA,C,N,O"
        elif mode == "HA":
            return "!@H="
        else:
            return "@CA"


# =============================================================================
# AmberInput - Input File Generation
# =============================================================================

@dataclass
class AmberInput:
    """
    Amber input file (.in) content.

    Attributes
    ----------
    title : str
        Job title
    cntrl : dict[str, Any]
        Control namelist parameters
    restraints : str | None
        Restraint specification
    """
    title: str = "Amber MD"
    cntrl: dict[str, Any] = field(default_factory=dict)
    restraints: str | None = None

    def to_string(self) -> str:
        """Generate Amber input file content."""
        lines = [self.title, "&cntrl"]

        # Format control parameters
        for key, value in self.cntrl.items():
            if isinstance(value, bool):
                value = 1 if value else 0
            lines.append(f"  {key}={value},")

        lines.append("&end")

        if self.restraints:
            lines.append(self.restraints)

        return "\n".join(lines) + "\n"

    def __str__(self) -> str:
        return self.to_string()


# =============================================================================
# AmberEngine - MD Input Generation
# =============================================================================

@dataclass
class AmberEngine:
    """
    Amber MD engine interface.

    Generates input files for Amber simulations.

    Attributes
    ----------
    exe : str
        Executable name (pmemd, pmemd.cuda, sander)
    """
    exe: str = "pmemd.cuda"

    def minimization_input(
        self,
        maxcyc: int = 5000,
        ncyc: int = 2500,
        cut: float = 10.0,
        restraint_wt: float = 0.0,
        restraint_mask: str | None = None,
    ) -> str:
        """
        Generate minimization input.

        Parameters
        ----------
        maxcyc : int
            Maximum minimization cycles
        ncyc : int
            Cycles of steepest descent before conjugate gradient
        cut : float
            Nonbonded cutoff
        restraint_wt : float
            Restraint weight (0 = no restraints)
        restraint_mask : str
            Amber mask for restrained atoms

        Returns
        -------
        str
            Amber input file content
        """
        cntrl = {
            "imin": 1,
            "maxcyc": maxcyc,
            "ncyc": ncyc,
            "cut": cut,
            "ntb": 1,
            "ntr": 1 if restraint_wt > 0 else 0,
        }

        if restraint_wt > 0:
            cntrl["restraint_wt"] = restraint_wt
            cntrl["restraintmask"] = f"'{restraint_mask or '@CA'}'"

        inp = AmberInput(title="Minimization", cntrl=cntrl)
        return inp.to_string()

    def heating_input(
        self,
        nsteps: int = 25000,
        dt: float = 0.002,
        temp_start: float = 0.0,
        temp_end: float = 300.0,
        cut: float = 10.0,
        restraint_wt: float = 5.0,
        restraint_mask: str = "@CA",
    ) -> str:
        """
        Generate NVT heating input.
        """
        cntrl = {
            "imin": 0,
            "irest": 0,
            "ntx": 1,
            "nstlim": nsteps,
            "dt": dt,
            "ntc": 2,
            "ntf": 2,
            "cut": cut,
            "ntb": 1,
            "ntp": 0,
            "ntt": 3,
            "gamma_ln": 2.0,
            "tempi": temp_start,
            "temp0": temp_end,
            "ntwx": 500,
            "ntwr": 5000,
            "ntpr": 500,
            "ioutfm": 1,
            "ntxo": 2,
            "ntr": 1 if restraint_wt > 0 else 0,
        }

        if restraint_wt > 0:
            cntrl["restraint_wt"] = restraint_wt
            cntrl["restraintmask"] = f"'{restraint_mask}'"

        inp = AmberInput(title="NVT Heating", cntrl=cntrl)
        return inp.to_string()

    def equilibration_input(
        self,
        nsteps: int = 50000,
        dt: float = 0.002,
        temp: float = 300.0,
        pressure: float = 1.0,
        cut: float = 10.0,
        restraint_wt: float = 1.0,
        restraint_mask: str = "@CA",
    ) -> str:
        """
        Generate NPT equilibration input.
        """
        cntrl = {
            "imin": 0,
            "irest": 1,
            "ntx": 5,
            "nstlim": nsteps,
            "dt": dt,
            "ntc": 2,
            "ntf": 2,
            "cut": cut,
            "ntb": 2,
            "ntp": 1,
            "barostat": 2,
            "pres0": pressure,
            "ntt": 3,
            "gamma_ln": 2.0,
            "temp0": temp,
            "ntwx": 500,
            "ntwr": 5000,
            "ntpr": 500,
            "ioutfm": 1,
            "ntxo": 2,
            "ntr": 1 if restraint_wt > 0 else 0,
        }

        if restraint_wt > 0:
            cntrl["restraint_wt"] = restraint_wt
            cntrl["restraintmask"] = f"'{restraint_mask}'"

        inp = AmberInput(title="NPT Equilibration", cntrl=cntrl)
        return inp.to_string()

    def production_input(
        self,
        nsteps: int = 500000,
        dt: float = 0.002,
        temp: float = 300.0,
        pressure: float = 1.0,
        cut: float = 10.0,
        ntwx: int = 500,
        ntpr: int = 500,
    ) -> str:
        """
        Generate production MD input.
        """
        cntrl = {
            "imin": 0,
            "irest": 1,
            "ntx": 5,
            "nstlim": nsteps,
            "dt": dt,
            "ntc": 2,
            "ntf": 2,
            "cut": cut,
            "ntb": 2,
            "ntp": 1,
            "barostat": 2,
            "pres0": pressure,
            "ntt": 3,
            "gamma_ln": 2.0,
            "temp0": temp,
            "ntwx": ntwx,
            "ntwr": 10000,
            "ntpr": ntpr,
            "ioutfm": 1,
            "ntxo": 2,
            "ntr": 0,
        }

        inp = AmberInput(title="Production MD", cntrl=cntrl)
        return inp.to_string()

    def run_command(
        self,
        topology: Path,
        coordinates: Path,
        input_file: Path,
        output_prefix: str,
        restart: Path | None = None,
    ) -> str:
        """
        Generate command to run Amber.
        """
        ref = coordinates if restart is None else restart

        cmd = (
            f"{self.exe} -O "
            f"-i {input_file} "
            f"-o {output_prefix}.out "
            f"-p {topology} "
            f"-c {coordinates} "
            f"-r {output_prefix}.rst7 "
            f"-x {output_prefix}.nc "
            f"-ref {ref} "
            f"-inf {output_prefix}.mdinfo"
        )

        return cmd


# =============================================================================
# Utility Functions
# =============================================================================

def ambpdb(
    topology: Path | str,
    coordinates: Path | str,
    output: Path | str,
) -> bool:
    """
    Convert Amber topology/coordinates to PDB using ambpdb.
    
    Parameters
    ----------
    topology : Path | str
        Amber topology file.
    coordinates : Path | str
        Amber coordinate file.
    output : Path | str
        Output PDB file.
        
    Returns
    -------
    bool
        True if successful.
    """
    exe_path = get_amber_exe()
    if not exe_path:
        log.error("AMBERHOME not set")
        return False
        
    ambpdb_exe = exe_path / "ambpdb"
    if not ambpdb_exe.exists():
        # Try cpptraj as fallback
        ambpdb_exe = exe_path / "cpptraj"
        
    if not ambpdb_exe.exists():
        log.error("ambpdb/cpptraj not found")
        return False
        
    try:
        with open(coordinates) as crd_file:
            result = subprocess.run(
                [str(ambpdb_exe), "-p", str(topology)],
                stdin=crd_file,
                capture_output=True,
                text=True,
            )
            
        if result.returncode == 0:
            Path(output).write_text(result.stdout)
            return True
        else:
            log.error(f"ambpdb failed: {result.stderr}")
            return False
            
    except Exception as e:
        log.error(f"Error running ambpdb: {e}")
        return False


def generate_ss_leap_commands(
    unit_name: str,
    ss_pairs: list[tuple[int, int]],
) -> list[str]:
    """
    Generate tleap commands for disulfide bonds.

    Equivalent to Biskit's AmberParmBuilder.__fLines(ss_bond, ss).

    Parameters
    ----------
    unit_name : str
        Leap unit name (e.g., 'mol')
    ss_pairs : list[tuple[int, int]]
        List of (residue_idx1, residue_idx2) pairs for disulfides.
        Indices should be 1-based (leap convention).

    Returns
    -------
    list[str]
        Leap commands to create disulfide bonds

    Examples
    --------
    >>> ss_pairs = [(10, 25), (33, 48)]
    >>> cmds = generate_ss_leap_commands('protein', ss_pairs)
    >>> # ['bond protein.10.SG protein.25.SG', 'bond protein.33.SG protein.48.SG']
    """
    commands = []
    for res1, res2 in ss_pairs:
        cmd = f"bond {unit_name}.{res1}.SG {unit_name}.{res2}.SG"
        commands.append(cmd)
        log.debug(f"SS bond command: {cmd}")
    return commands
