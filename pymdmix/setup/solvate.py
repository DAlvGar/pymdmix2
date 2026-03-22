"""
Solvation
=========

Solvate structures with organic solvent mixtures using LEaP.

Examples
--------
>>> from pymdmix.setup import solvate_structure
>>> from pymdmix.core import load_structure, SolventLibrary
>>>
>>> struct = load_structure("protein_prepared.pdb")
>>> library = SolventLibrary()
>>> solvent = library.get("ETA")
>>>
>>> solvated = solvate_structure(struct, solvent, box_buffer=12.0)
"""

from __future__ import annotations

import logging
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path

import parmed

from pymdmix.core.solvent import Solvent

log = logging.getLogger(__name__)


@dataclass
class SolvationOptions:
    """
    High-level solvation options passed to :func:`solvate_structure`.

    Attributes
    ----------
    box_buffer : float
        Minimum distance (Å) between solute and box edge (default 12.0).
    box_shape : str
        Box geometry: ``"truncated_octahedron"`` (default), ``"octahedron"``,
        or ``"box"``.
    neutralize : bool
        Add counter-ions to neutralise the system (default True).
    ion_concentration : float
        Additional NaCl concentration in mol/L (default 0.15).
    cation : str
        Cation residue name for LEaP (default ``"Na+"``).
    anion : str
        Anion residue name for LEaP (default ``"Cl-"``).
    extra_forcefields : list[str]
        Additional force-field files to source in LEaP.
    """

    box_buffer: float = 12.0
    box_shape: str = "truncated_octahedron"
    neutralize: bool = True
    ion_concentration: float = 0.15
    cation: str = "Na+"
    anion: str = "Cl-"
    extra_forcefields: list[str] = field(default_factory=list)


@dataclass
class SolvateResult:
    """Result of solvation."""

    topology: Path | None = None
    coordinates: Path | None = None
    n_solvent_residues: int = 0
    n_ions_added: int = 0
    box_dimensions: tuple[float, float, float] | None = None
    leap_log: str = ""
    success: bool = False
    error: str | None = None

    def save_coordinates(self, path: str | Path) -> None:
        """
        Copy the coordinate file to *path*.

        Parameters
        ----------
        path : str or Path
            Destination file path.

        Raises
        ------
        ValueError
            If no coordinates were produced (solvation failed).
        """
        if self.coordinates is None:
            raise ValueError("No coordinates available — solvation may have failed")
        shutil.copy2(self.coordinates, Path(path))

    def save_topology(self, path: str | Path) -> None:
        """
        Copy the topology file to *path*.

        Parameters
        ----------
        path : str or Path
            Destination file path.

        Raises
        ------
        ValueError
            If no topology was produced (solvation failed).
        """
        if self.topology is None:
            raise ValueError("No topology available — solvation may have failed")
        shutil.copy2(self.topology, Path(path))


@dataclass
class BoxConfig:
    """Solvation box configuration."""

    shape: str = "truncated_octahedron"  # or "box", "octahedron"
    buffer: float = 12.0  # Angstroms from solute to box edge
    closeness: float = 1.0  # Minimum distance between solute and solvent

    def to_leap_command(self, unit_name: str, solvent_box: str) -> str:
        """Generate LEaP solvation command."""
        if self.shape == "truncated_octahedron":
            return f"solvateoct {unit_name} {solvent_box} {self.buffer} {self.closeness}"
        elif self.shape == "octahedron":
            return f"solvateoct {unit_name} {solvent_box} {self.buffer} {self.closeness}"
        else:  # box
            return f"solvatebox {unit_name} {solvent_box} {self.buffer} {self.closeness}"


@dataclass
class IonConfig:
    """Ion addition configuration."""

    neutralize: bool = True
    add_ions: bool = True
    cation: str = "Na+"
    anion: str = "Cl-"
    concentration: float = 0.15  # Molar

    def to_leap_commands(self, unit_name: str) -> list[str]:
        """Generate LEaP ion commands."""
        commands = []

        if self.neutralize:
            # First neutralize
            commands.append(f"addions {unit_name} {self.cation} 0")
            commands.append(f"addions {unit_name} {self.anion} 0")

        if self.add_ions and self.concentration > 0:
            # Add to target concentration
            # This is approximate - LEaP uses count, not concentration
            commands.append(f"# Target concentration: {self.concentration} M")
            # Would need box volume to calculate exact ion count

        return commands


def generate_leap_script(
    pdb_path: Path,
    solvent: Solvent,
    output_prefix: str,
    box_config: BoxConfig | None = None,
    ion_config: IonConfig | None = None,
    force_fields: list[str] | None = None,
    extra_commands: list[str] | None = None,
) -> str:
    """
    Generate LEaP script for solvation.

    Parameters
    ----------
    pdb_path : Path
        Input PDB file
    solvent : Solvent
        Solvent mixture definition
    output_prefix : str
        Prefix for output files
    box_config : BoxConfig
        Box configuration
    ion_config : IonConfig
        Ion configuration
    force_fields : list[str]
        Force field files to load
    extra_commands : list[str]
        Additional LEaP commands

    Returns
    -------
    str
        LEaP script content
    """
    if box_config is None:
        box_config = BoxConfig()
    if ion_config is None:
        ion_config = IonConfig()
    if force_fields is None:
        force_fields = [
            "leaprc.protein.ff19SB",
            "leaprc.water.opc",
            "leaprc.gaff2",
        ]

    lines = [
        "# LEaP script generated by pyMDMix",
        "",
        "# Load force fields",
    ]

    for ff in force_fields:
        lines.append(f"source {ff}")

    lines.extend(
        [
            "",
            "# Load solvent parameters",
        ]
    )

    if solvent.off_file and solvent.off_file.exists():
        lines.append(f"loadoff {solvent.off_file}")

    lines.extend(
        [
            "",
            "# Load structure",
            f"sys = loadpdb {pdb_path}",
            "",
            "# Solvate",
        ]
    )

    # Determine solvent box name
    solvent_box = f"{solvent.name}BOX"
    lines.append(box_config.to_leap_command("sys", solvent_box))

    # Add ions
    lines.extend(["", "# Add ions"])
    lines.extend(ion_config.to_leap_commands("sys"))

    # Extra commands
    if extra_commands:
        lines.extend(["", "# Additional commands"])
        lines.extend(extra_commands)

    # Save outputs
    lines.extend(
        [
            "",
            "# Check and save",
            "check sys",
            f"saveamberparm sys {output_prefix}.prmtop {output_prefix}.inpcrd",
            f"savepdb sys {output_prefix}.pdb",
            "",
            "quit",
        ]
    )

    return "\n".join(lines)


def run_leap(
    script: str,
    leap_exe: str = "tleap",
    work_dir: Path | None = None,
) -> tuple[bool, str]:
    """
    Run LEaP with a script.

    Parameters
    ----------
    script : str
        LEaP script content
    leap_exe : str
        Path to tleap executable
    work_dir : Path
        Working directory

    Returns
    -------
    tuple[bool, str]
        (success, log_output)
    """
    if work_dir is None:
        work_dir = Path.cwd()

    # Write script to temp file
    script_path = work_dir / "leap.in"
    script_path.write_text(script)

    try:
        result = subprocess.run(
            [leap_exe, "-f", str(script_path)],
            cwd=work_dir,
            capture_output=True,
            text=True,
            timeout=300,
        )

        log_output = result.stdout + result.stderr
        success = result.returncode == 0

        # Check for errors in output
        if "FATAL" in log_output or "Error" in log_output:
            success = False

        return success, log_output

    except FileNotFoundError:
        return False, f"LEaP executable not found: {leap_exe}"
    except subprocess.TimeoutExpired:
        return False, "LEaP timed out"
    except Exception as e:
        return False, f"LEaP error: {e}"


def solvate_structure(
    structure: parmed.Structure | Path,
    solvent: Solvent,
    output_dir: Path | None = None,
    output_prefix: str = "system",
    options: SolvationOptions | None = None,
    box_config: BoxConfig | None = None,
    ion_config: IonConfig | None = None,
    leap_exe: str = "tleap",
) -> SolvateResult:
    """
    Solvate a structure with a solvent mixture.

    Parameters
    ----------
    structure : parmed.Structure or Path
        Input structure or PDB path
    solvent : Solvent
        Solvent mixture
    output_dir : Path
        Output directory
    output_prefix : str
        Prefix for output files
    box_config : BoxConfig
        Box configuration
    ion_config : IonConfig
        Ion configuration
    leap_exe : str
        Path to tleap

    Returns
    -------
    SolvateResult
        Solvation result with file paths
    """
    # Apply SolvationOptions to BoxConfig / IonConfig when provided
    if options is not None:
        if box_config is None:
            box_config = BoxConfig(
                shape=options.box_shape,
                buffer=options.box_buffer,
            )
        if ion_config is None:
            ion_config = IonConfig(
                neutralize=options.neutralize,
                add_ions=options.ion_concentration > 0,
                cation=options.cation,
                anion=options.anion,
                concentration=options.ion_concentration,
            )

    output_dir = Path(output_dir) if output_dir is not None else Path.cwd()
    output_dir.mkdir(parents=True, exist_ok=True)

    result = SolvateResult()

    # Save structure to PDB if needed
    if isinstance(structure, parmed.Structure):
        pdb_path = output_dir / "input.pdb"
        structure.save(str(pdb_path), overwrite=True)
    else:
        pdb_path = Path(structure)

    if not pdb_path.exists():
        result.error = f"Input PDB not found: {pdb_path}"
        return result

    # Generate LEaP script
    script = generate_leap_script(
        pdb_path=pdb_path,
        solvent=solvent,
        output_prefix=output_prefix,
        box_config=box_config,
        ion_config=ion_config,
    )

    log.info(f"Running LEaP for {solvent.name} solvation")

    # Run LEaP
    success, leap_log = run_leap(script, leap_exe, output_dir)
    result.leap_log = leap_log

    if not success:
        result.error = "LEaP failed - check leap_log for details"
        log.error(f"LEaP failed: {leap_log[:500]}")
        return result

    # Check outputs
    topology = output_dir / f"{output_prefix}.prmtop"
    coordinates = output_dir / f"{output_prefix}.inpcrd"

    if topology.exists() and coordinates.exists():
        result.success = True
        result.topology = topology
        result.coordinates = coordinates
        log.info(f"Solvation complete: {topology}")
    else:
        result.error = "Output files not created"

    return result
