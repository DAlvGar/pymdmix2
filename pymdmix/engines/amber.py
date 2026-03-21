"""
Amber MD Engine
===============

Interface for running Amber molecular dynamics simulations.

Generates input files for:
- Minimization
- Equilibration (NVT, NPT)
- Production MD

Examples
--------
>>> from pymdmix.engines import AmberEngine
>>>
>>> engine = AmberEngine(pmemd_exe="pmemd.cuda")
>>>
>>> # Generate production input
>>> mdin = engine.production_input(nsteps=500000)
>>> Path("prod.in").write_text(mdin)
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any
import logging

log = logging.getLogger(__name__)


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
        params = []
        for key, value in self.cntrl.items():
            if isinstance(value, bool):
                value = 1 if value else 0
            elif isinstance(value, float):
                value = f"{value:.4f}" if value != int(value) else str(int(value))
            params.append(f"  {key}={value}")

        lines.extend(params)
        lines.append("&end")

        if self.restraints:
            lines.extend(["", self.restraints])

        return ",\n".join(lines[:2]) + ",\n" + ",\n".join(lines[2:-1]) + "\n" + lines[-1] + "\n"

    def __str__(self) -> str:
        return self.to_string()


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
        return self._format_mdin(inp)

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

        Parameters
        ----------
        nsteps : int
            Number of MD steps
        dt : float
            Timestep in ps
        temp_start : float
            Starting temperature
        temp_end : float
            Target temperature
        cut : float
            Nonbonded cutoff
        restraint_wt : float
            Restraint weight
        restraint_mask : str
            Restraint mask

        Returns
        -------
        str
            Amber input file content
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
            "ntp": 0,  # NVT
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
        return self._format_mdin(inp)

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

        Parameters
        ----------
        nsteps : int
            Number of MD steps
        dt : float
            Timestep in ps
        temp : float
            Temperature
        pressure : float
            Pressure in bar
        cut : float
            Nonbonded cutoff
        restraint_wt : float
            Restraint weight
        restraint_mask : str
            Restraint mask

        Returns
        -------
        str
            Amber input file content
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
        return self._format_mdin(inp)

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

        Parameters
        ----------
        nsteps : int
            Number of MD steps
        dt : float
            Timestep in ps
        temp : float
            Temperature
        pressure : float
            Pressure in bar
        cut : float
            Nonbonded cutoff
        ntwx : int
            Trajectory write frequency
        ntpr : int
            Energy output frequency

        Returns
        -------
        str
            Amber input file content
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
        return self._format_mdin(inp)

    def _format_mdin(self, inp: AmberInput) -> str:
        """Format AmberInput to proper mdin format."""
        lines = [inp.title, "&cntrl"]

        # Format parameters
        for key, value in inp.cntrl.items():
            if isinstance(value, bool):
                value = 1 if value else 0
            lines.append(f"  {key}={value},")

        lines.append("&end")

        if inp.restraints:
            lines.append(inp.restraints)

        return "\n".join(lines) + "\n"

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

        Parameters
        ----------
        topology : Path
            Topology file
        coordinates : Path
            Coordinate/restart file
        input_file : Path
            Input file
        output_prefix : str
            Prefix for output files
        restart : Path
            Restart file (for continuation)

        Returns
        -------
        str
            Shell command
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
