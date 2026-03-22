"""
Replica Management
==================

A Replica represents a single MD simulation run.
Multiple replicas with different random seeds improve sampling.

Replicas track:
- Setup state (topology, coordinates)
- Simulation progress
- Analysis status

Examples
--------
>>> from pymdmix.project import Replica, ReplicaState
>>>
>>> replica = Replica(
...     name="ETA_1",
...     solvent="ETA",
...     path=Path("replicas/ETA_1"),
... )
>>>
>>> # Check state
>>> if replica.state == ReplicaState.READY:
...     replica.submit()
"""

from __future__ import annotations

import logging
import re
import shutil
from dataclasses import dataclass, field
from datetime import datetime
from enum import Enum, auto
from pathlib import Path
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    import parmed

    from pymdmix.core.grid import Grid
    from pymdmix.core.trajectory import TrajectoryReader

log = logging.getLogger(__name__)


class ChainedTrajectoryReader:
    """Chain multiple trajectory files into a single reader interface."""

    def __init__(
        self,
        topology: Path,
        trajectories: list[Path],
        frame_step: int = 1,
    ):
        from pymdmix.core.trajectory import FrameSliceReader, open_trajectory

        if not trajectories:
            raise ReplicaError("No trajectory files selected")

        self._readers = [
            FrameSliceReader(open_trajectory(topology, traj), step=frame_step)
            for traj in trajectories
        ]

    @property
    def n_frames(self) -> int:
        return sum(reader.n_frames for reader in self._readers)

    @property
    def n_atoms(self) -> int:
        return int(self._readers[0].n_atoms)

    def __len__(self) -> int:
        return self.n_frames

    def __iter__(self):
        for reader in self._readers:
            yield from reader


class ReplicaState(Enum):
    """
    State of a replica in the workflow.

    States progress: CREATED → SETUP → READY → RUNNING → COMPLETE → ANALYZED
    """

    CREATED = auto()  # Just created, no files
    SETUP = auto()  # Topology/coordinates generated
    READY = auto()  # Ready to submit
    RUNNING = auto()  # Simulation in progress
    COMPLETE = auto()  # Simulation finished
    ALIGNED = auto()  # Trajectory aligned
    ANALYZED = auto()  # Analysis complete
    ERROR = auto()  # Error state


class ReplicaError(Exception):
    """Base exception for replica operations."""

    pass


class BadFileError(ReplicaError):
    """Invalid or missing file."""

    pass


# Default folder names
DEFAULT_FOLDERS = {
    "min": "min",
    "eq": "eq",
    "md": "md",
    "align": "align",
    "density": "density",
    "energy": "energy",
}

# Known trajectory extensions
TRAJECTORY_EXTENSIONS = ("nc", "netcdf", "dcd", "mdcrd", "crd", "xtc", "trr")


@dataclass
class MDSettings:
    """
    MD simulation settings for a replica.

    Attributes
    ----------
    nanos : int
        Simulation length in nanoseconds
    temperature : float
        Temperature in Kelvin
    timestep : float
        Integration timestep in femtoseconds
    prod_steps : int
        Steps per production output file
    traj_frequency : int
        Trajectory save frequency (steps)
    restraint_mode : str
        Restraint mode: FREE, BB (backbone), HA (heavy atoms), CUSTOM
    restraint_force : float
        Restraint force constant (kcal/mol/Å²)
    restraint_mask : str
        Amber mask for custom restraints
    align_mask : str
        Amber mask for trajectory alignment
    """

    nanos: int = 20
    temperature: float = 300.0
    timestep: float = 2.0
    prod_steps: int = 500000  # Steps per output file (1ns at 2fs)
    traj_frequency: int = 1000
    restraint_mode: str = "FREE"
    restraint_force: float = 0.0
    restraint_mask: str = ""
    align_mask: str = ""
    md_program: str = "AMBER"

    @property
    def has_restraints(self) -> bool:
        """Check if restraints are enabled."""
        return self.restraint_mode != "FREE" and self.restraint_force > 0

    @property
    def n_traj_files(self) -> int:
        """Expected number of trajectory files."""
        # nanos * 1e6 fs / timestep / prod_steps
        return int((self.nanos * 1e6 / self.timestep) / self.prod_steps)

    @property
    def n_snapshots(self) -> int:
        """Total expected snapshots."""
        return self.n_traj_files * (self.prod_steps // self.traj_frequency)


@dataclass
class Replica:
    """
    A single MD simulation replica.

    Attributes
    ----------
    name : str
        Replica identifier (e.g., "ETA_1")
    solvent : str
        Solvent mixture name
    path : Path | None
        Directory containing replica files
    state : ReplicaState
        Current workflow state
    topology : str | None
        Topology filename (relative to path)
    coordinates : str | None
        Coordinate filename
    trajectory : str | None
        Trajectory filename
    reference : str | None
        Reference PDB for alignment
    seed : int
        Random seed for simulation
    settings : MDSettings | None
        MD simulation settings
    """

    name: str
    solvent: str
    path: Path | None = None
    state: ReplicaState = ReplicaState.CREATED

    # Files (relative to path)
    topology: str | None = None
    coordinates: str | None = None
    trajectory: str | None = None
    reference: str | None = None
    pdb: str | None = None

    # Simulation info
    seed: int = 0
    n_steps: int = 0

    # MD Settings
    settings: MDSettings | None = field(default_factory=MDSettings)

    # Folder names
    min_folder: str = "min"
    eq_folder: str = "eq"
    md_folder: str = "md"
    align_folder: str = "align"
    density_folder: str = "density"
    energy_folder: str = "energy"

    # Output templates
    md_output_template: str = "md_{step:03d}.{extension}"
    eq_output_template: str = "eq_{step:02d}.{extension}"

    # Timestamps
    created_at: str = field(default_factory=lambda: datetime.now().isoformat())
    modified_at: str = field(default_factory=lambda: datetime.now().isoformat())

    # Metadata
    metadata: dict[str, Any] = field(default_factory=dict)

    # Internal state
    _folders_created: bool = field(default=False, repr=False)
    _md_input_written: bool = field(default=False, repr=False)
    _grids: list | None = field(default=None, repr=False)

    def __post_init__(self):
        """Normalize path and settings."""
        if self.path and isinstance(self.path, str):
            self.path = Path(self.path)

        if isinstance(self.state, str):
            self.state = ReplicaState[self.state]

        if self.settings is None:
            self.settings = MDSettings()
        elif isinstance(self.settings, dict):
            self.settings = MDSettings(**self.settings)

    # =========================================================================
    # Path Properties
    # =========================================================================

    @property
    def topology_path(self) -> Path | None:
        """Full path to topology file."""
        if self.path and self.topology:
            return self.path / self.topology
        return None

    @property
    def coordinates_path(self) -> Path | None:
        """Full path to coordinates file."""
        if self.path and self.coordinates:
            return self.path / self.coordinates
        return None

    @property
    def trajectory_path(self) -> Path | None:
        """Full path to trajectory file."""
        if self.path and self.trajectory:
            return self.path / self.trajectory
        return None

    @property
    def reference_path(self) -> Path | None:
        """Full path to reference PDB."""
        if self.path and self.reference:
            return self.path / self.reference
        return None

    @property
    def pdb_path(self) -> Path | None:
        """Full path to PDB file."""
        if self.path and self.pdb:
            return self.path / self.pdb
        return None

    @property
    def min_path(self) -> Path | None:
        """Full path to minimization folder."""
        if self.path:
            return self.path / self.min_folder
        return None

    @property
    def eq_path(self) -> Path | None:
        """Full path to equilibration folder."""
        if self.path:
            return self.path / self.eq_folder
        return None

    @property
    def md_path(self) -> Path | None:
        """Full path to production MD folder."""
        if self.path:
            return self.path / self.md_folder
        return None

    @property
    def align_path(self) -> Path | None:
        """Full path to aligned trajectory folder."""
        if self.path:
            return self.path / self.align_folder
        return None

    @property
    def density_path(self) -> Path | None:
        """Full path to density grids folder."""
        if self.path:
            return self.path / self.density_folder
        return None

    @property
    def energy_path(self) -> Path | None:
        """Full path to energy grids folder."""
        if self.path:
            return self.path / self.energy_folder
        return None

    # =========================================================================
    # Status Checks
    # =========================================================================

    def exists(self) -> bool:
        """Check if replica directory exists."""
        return self.path is not None and self.path.exists()

    def has_topology(self) -> bool:
        """Check if topology file exists."""
        return self.topology_path is not None and self.topology_path.exists()

    def has_trajectory(self) -> bool:
        """Check if trajectory file exists."""
        return self.trajectory_path is not None and self.trajectory_path.exists()

    def folders_created(self) -> bool:
        """Return True if replica directory structure is created."""
        return self._folders_created or (self.path is not None and self.path.exists())

    def md_input_written(self) -> bool:
        """Return True if MD input files have been written."""
        return self._md_input_written

    def is_created(self) -> bool:
        """Return True if replica folder and MD inputs are ready."""
        return self.folders_created() and self.md_input_written()

    def is_minimization_finished(self) -> bool:
        """Check if minimization stage is complete."""
        checker = self.get_checker(warn=False)
        if checker:
            return bool(checker.check_minimization())
        # Fallback: check for min output files
        if self.min_path and self.min_path.exists():
            out_files = list(self.min_path.glob("*.out")) + list(self.min_path.glob("*.rst7"))
            return len(out_files) > 0
        return False

    def is_equilibration_finished(self) -> bool:
        """Check if equilibration stage is complete."""
        checker = self.get_checker(warn=False)
        if checker:
            return bool(checker.check_equilibration())
        # Fallback: check for eq output files
        if self.eq_path and self.eq_path.exists():
            out_files = list(self.eq_path.glob("*.out")) + list(self.eq_path.glob("*.rst7"))
            return len(out_files) >= 2  # Typically eq1, eq2
        return False

    def is_production_finished(self, step_selection: list[int] | None = None) -> bool:
        """
        Check if MD production stage is complete.

        Parameters
        ----------
        step_selection : list[int] | None
            Specific steps to check. None checks all expected steps.

        Returns
        -------
        bool
            True if all selected production steps are complete.
        """
        if not self.settings:
            return False

        if step_selection is None:
            step_selection = list(range(1, self.settings.n_traj_files + 1))

        checker = self.get_checker(warn=False)
        if checker:
            return bool(checker.check_production(step_selection=step_selection))

        # Fallback: check trajectory files exist
        extensions = self.check_production_extension(step_selection)
        return all(ext is not None for ext in extensions.values())

    def is_aligned(self, step_selection: list[int] | None = None) -> bool:
        """
        Check if trajectory has been aligned.

        Parameters
        ----------
        step_selection : list[int] | None
            Specific steps to check.

        Returns
        -------
        bool
            True if all selected steps are aligned.
        """
        if not self.settings:
            return False

        if step_selection is None:
            step_selection = list(range(1, self.settings.n_traj_files + 1))

        extensions = self.check_align_extension(step_selection)
        return all(ext is not None for ext in extensions.values())

    def last_completed_production_step(self, start_step: int = 1) -> int:
        """
        Find the last completed production step.

        Useful for tracking progress or handling incomplete runs.

        Parameters
        ----------
        start_step : int
            Start checking from this step number.

        Returns
        -------
        int
            Last completed step number, or 0 if none complete.
        """
        if not self.settings:
            return 0

        last = 0
        for step in range(start_step, self.settings.n_traj_files + 1):
            if self.is_production_finished([step]):
                last = step
            else:
                break
        return last

    def check_production_extension(self, steps: list[int] | None = None) -> dict[int, str | None]:
        """
        Check file extensions for production trajectories.

        Parameters
        ----------
        steps : list[int] | None
            Steps to check. None checks all.

        Returns
        -------
        dict[int, str | None]
            Mapping of step number to extension (or None if not found).
        """
        if not self.md_path or not self.md_path.exists():
            return {}

        if steps is None and self.settings:
            steps = list(range(1, self.settings.n_traj_files + 1))
        elif steps is None:
            steps = []

        result = {}
        files = list(self.md_path.iterdir())

        for step in steps:
            # Build pattern from template
            base = self.md_output_template.format(step=step, extension="")
            pattern = re.compile(re.escape(base) + r"(\w+)")

            matches = []
            for f in files:
                m = pattern.match(f.name)
                if m:
                    ext = m.group(1)
                    if ext in TRAJECTORY_EXTENSIONS:
                        matches.append(ext)

            result[step] = matches[0] if matches else None

        return result

    def check_align_extension(self, steps: list[int] | None = None) -> dict[int, str | None]:
        """
        Check file extensions for aligned trajectories.

        Parameters
        ----------
        steps : list[int] | None
            Steps to check.

        Returns
        -------
        dict[int, str | None]
            Mapping of step number to extension.
        """
        if not self.align_path or not self.align_path.exists():
            if steps:
                return {s: None for s in steps}
            return {}

        if steps is None and self.settings:
            steps = list(range(1, self.settings.n_traj_files + 1))
        elif steps is None:
            steps = []

        result = {}
        files = list(self.align_path.iterdir())

        for step in steps:
            base = self.md_output_template.format(step=step, extension="")
            pattern = re.compile(re.escape(base) + r"(\w+)")

            matches = []
            for f in files:
                m = pattern.match(f.name)
                if m:
                    ext = m.group(1)
                    if ext in TRAJECTORY_EXTENSIONS:
                        matches.append(ext)

            result[step] = matches[0] if matches else None

        return result

    # =========================================================================
    # Folder Creation
    # =========================================================================

    def create_directory(self) -> None:
        """Create replica directory structure."""
        if self.path is None:
            raise ValueError("Replica path not set")

        self.path.mkdir(parents=True, exist_ok=True)

        # Create subdirectories
        (self.path / self.min_folder).mkdir(exist_ok=True)
        (self.path / self.eq_folder).mkdir(exist_ok=True)
        (self.path / self.md_folder).mkdir(exist_ok=True)

        self._folders_created = True
        self._update_modified()
        log.info(f"Created replica directory: {self.path}")

    def create_folder(
        self,
        where: Path | str | None = None,
        system: Any = None,
        fix_topology: bool = True,
    ) -> None:
        """
        Create directory tree and copy/generate system files.

        Parameters
        ----------
        where : Path | str | None
            Parent directory. Defaults to current directory.
        system : SolvatedSystem | None
            System object with topology/coordinates to copy.
        fix_topology : bool
            Remove SCEE/SCNB sections from topology (Amber compatibility).
        """
        if not self.name:
            raise ReplicaError("Unnamed replica - cannot create folder")

        # Determine base path
        if where:
            where = Path(where).resolve()
        else:
            where = Path.cwd()

        self.path = where / self.name

        # Create directory structure
        self.create_directory()

        # Copy system files if provided
        if system is not None:
            self._setup_system_files(system, fix_topology)

        # Save replica state
        self.save()
        log.info(f"Created folder structure for replica {self.name}")

    def _setup_system_files(self, system: Any, fix_topology: bool = True) -> None:
        """Copy system files to replica directory."""
        if self.path is None:
            raise ReplicaError("Replica path not set")

        basename = f"{system.name}_{self.name}"

        # Save topology and coordinates
        if hasattr(system, "save_top_crd"):
            system.save_top_crd(self.path / basename)
            self.topology = f"{basename}.prmtop"
            self.coordinates = f"{basename}.prmcrd"

        # Save PDB
        if hasattr(system, "save_pdb"):
            pdb_path = self.path / f"{basename}.pdb"
            system.save_pdb(pdb_path)
            self.pdb = f"{basename}.pdb"

        # Save reference
        if hasattr(system, "ref") and system.ref:
            ref_path = self.path / f"{basename}_ref.pdb"
            if hasattr(system.ref, "write_pdb"):
                system.ref.write_pdb(ref_path)
            self.reference = f"{basename}_ref.pdb"

        # Fix topology if requested
        if fix_topology and self.topology_path and self.topology_path.exists():
            self._fix_topology(self.topology_path)

    def _fix_topology(self, prmtop: Path) -> None:
        """
        Remove SCEE and SCNB sections from Amber topology.

        Needed for some solvent boxes that crash newer Amber versions.
        """
        backup = prmtop.with_suffix(".prmtop_back")
        shutil.copy(prmtop, backup)

        with open(backup) as f:
            lines = f.readlines()

        with open(prmtop, "w") as out:
            skip = False
            for line in lines:
                if "%FLAG SCEE_SCALE_FACTOR" in line:
                    skip = True
                elif skip and "%FLAG SOLTY" in line:
                    skip = False
                    out.write(line)
                elif not skip:
                    out.write(line)

        log.debug("Fixed topology: removed SCEE/SCNB sections")

    def create_md_input(self, **kwargs) -> bool:
        """
        Create MD input configuration files.

        Generates input files for the selected MD program (AMBER, NAMD, OpenMM).

        Returns
        -------
        bool
            True if successful.
        """
        if not self.folders_created():
            self.create_folder()

        if not self.settings:
            raise ReplicaError("MD settings not configured")

        program = self.settings.md_program.upper()
        writer: Any

        if program == "AMBER":
            from pymdmix.engines.amber import AmberWriter

            writer = AmberWriter(self)
        elif program == "NAMD":
            from pymdmix.engines.namd import NAMDWriter

            writer = NAMDWriter(self)
        elif program == "OPENMM":
            from pymdmix.engines.openmm import OpenMMWriter

            writer = OpenMMWriter(self)
        else:
            raise ReplicaError(f"Unknown MD program: {program}")

        writer.write_commands()
        writer.write_replica_input()

        # Apply H-mass repartitioning for 4fs timestep
        if self.settings.timestep == 4:
            log.info(f"4fs timestep: Applying H-mass repartitioning to {self.name}")
            self.apply_hmass_repartitioning()

        self._md_input_written = True
        self.save()
        return True

    def create_queue_input(self, queue_system: str, **kwargs) -> bool:
        """
        Create queue submission scripts.

        Parameters
        ----------
        queue_system : str
            Queue system name (slurm, pbs, sge, lsf).

        Returns
        -------
        bool
            True if successful.
        """
        if not self.folders_created():
            return False

        from pymdmix.engines.queue import QueueConfig, generate_queue_script

        log.info(f"Writing {queue_system} scripts for replica {self.name}")

        config = QueueConfig(system=queue_system, **kwargs)
        if self.path is None:
            raise ReplicaError("Replica path not set")
        script = generate_queue_script(
            config=config,
            job_name=self.name,
            commands=[f"cd {self.path}", "bash COMMANDS.sh"],
        )

        script_path = self.path / f"submit_{queue_system}.sh"
        script_path.write_text(script)

        return True

    def create_all(self, queue: str | None = None, **kwargs) -> None:
        """
        Create folders, MD input, and optionally queue scripts.

        Parameters
        ----------
        queue : str | None
            Queue system for submission scripts.
        """
        self.create_folder(**kwargs)
        self.create_md_input(**kwargs)
        if queue:
            self.create_queue_input(queue, **kwargs)

    def apply_hmass_repartitioning(self) -> None:
        """Apply hydrogen mass repartitioning to topology."""
        if not self.topology_path or not self.topology_path.exists():
            raise ReplicaError("Topology file not found")

        try:
            import parmed as pmd

            parm = pmd.load_file(str(self.topology_path))
            action = pmd.tools.actions.HMassRepartition(parm)
            action.execute()

            # Backup and save
            backup = self.topology_path.with_suffix(".prmtop_back")
            shutil.move(self.topology_path, backup)
            parm.save(str(self.topology_path))

            log.info(f"Applied H-mass repartitioning to {self.topology}")

        except ImportError:
            log.warning("parmed not available - cannot apply H-mass repartitioning")

    # =========================================================================
    # Data Access
    # =========================================================================

    def get_checker(self, warn: bool = True) -> Any:
        """
        Get MD checker for the simulation program.

        Returns
        -------
        AmberChecker | NAMDChecker | OpenMMChecker | None
            Checker instance or None if not available.
        """
        if not self.settings:
            return None

        program = self.settings.md_program.upper()

        try:
            if program == "AMBER":
                from pymdmix.engines.amber import AmberChecker

                return AmberChecker(self, warn=warn)
            elif program == "NAMD":
                from pymdmix.engines.namd import NAMDChecker

                return NAMDChecker(self, warn=warn)
            elif program == "OPENMM":
                from pymdmix.engines.openmm import OpenMMChecker

                return OpenMMChecker(self, warn=warn)
        except ImportError:
            pass

        return None

    def get_solvent(self) -> Any:
        """Get solvent instance for this replica."""
        from pymdmix.core.solvent import SolventLibrary

        library = SolventLibrary()
        return library.get(self.solvent)

    def get_probes(self) -> list[str]:
        """Get list of probe names from solvent."""
        solvent = self.get_solvent()
        probes = [p.name for p in solvent.probes]
        # Add COM probes if any
        if hasattr(solvent, "com_probes"):
            probes.extend(solvent.com_probes.keys())
        return probes

    def get_trajectory(
        self,
        step_selection: list[int] | None = None,
        use_aligned: bool = True,
        frame_step: int = 1,
    ) -> TrajectoryReader:
        """
        Get trajectory for this replica.

        Parameters
        ----------
        step_selection : list[int] | None
            Specific step numbers to include.
        use_aligned : bool
            Use aligned trajectory if available.
        frame_step : int
            Take every Nth frame.

        Returns
        -------
        TrajectoryReader
            Trajectory object for analysis.
        """
        if not self.settings:
            raise ReplicaError("MD settings not configured")

        if step_selection is None:
            step_selection = list(range(1, self.settings.n_traj_files + 1))

        # Decide which trajectory to use
        if use_aligned and self.is_aligned(step_selection):
            path = self.align_path
            extensions = self.check_align_extension(step_selection)
            log.debug("Using aligned trajectory")
        else:
            if not self.is_production_finished(step_selection):
                raise ReplicaError(f"Production not finished for steps: {step_selection}")
            path = self.md_path
            extensions = self.check_production_extension(step_selection)
            log.debug("Using production trajectory")

        if path is None:
            raise ReplicaError("Trajectory path not set")

        # Build file list
        files: list[Path] = []
        for step in sorted(step_selection):
            ext = extensions.get(step)
            if not ext:
                log.warning(f"File for step {step} not found")
                continue
            fname = self.md_output_template.format(step=step, extension=ext)
            files.append(path / fname)

        # Get topology
        topology = self.topology_path or self.pdb_path
        if topology is None:
            raise ReplicaError("Topology or PDB file not found for trajectory loading")

        return ChainedTrajectoryReader(topology, files, frame_step=frame_step)

    def get_pdb(self) -> parmed.Structure:
        """
        Get PDB structure for this replica.

        Returns
        -------
        Structure
            PDB structure with solvent information.
        """
        from pymdmix.core.structure import load_structure

        if not self.pdb_path or not self.pdb_path.exists():
            raise BadFileError(f"PDB file not found: {self.pdb_path}")

        structure = load_structure(self.pdb_path)
        return structure

    # =========================================================================
    # Grid Management
    # =========================================================================

    def fetch_grids(
        self,
        prefix: str | None = None,
        suffix: str | None = None,
        grid_type: str | None = None,
    ) -> list[Grid]:
        """
        Find and load grids from replica folders.

        Parameters
        ----------
        prefix : str | None
            Filter by filename prefix.
        suffix : str | None
            Filter by filename suffix (before extension).
        grid_type : str | None
            Filter by grid type (density, energy, etc.).

        Returns
        -------
        list[Grid]
            List of Grid objects found.
        """
        from pymdmix.core.grid import Grid

        if not self.path or not self.path.exists():
            return []

        prefix = prefix or ""
        suffix = suffix or ""

        grids = []

        # Search in replica folder and subfolders
        for dx_file in self.path.rglob("*.dx"):
            fname = dx_file.stem

            # Apply filters
            if prefix and not fname.startswith(prefix):
                continue
            if suffix and not fname.endswith(suffix):
                continue

            try:
                grid = Grid.read_dx(dx_file)
                # Try to parse probe/type from filename
                grid.metadata = {"path": dx_file, "name": fname}
                grids.append(grid)
            except Exception as e:
                log.warning(f"Failed to load grid {dx_file}: {e}")

        self._grids = grids
        return grids

    def get_grids_by_type(
        self,
        grid_type: str | None = None,
        **kwargs,
    ) -> dict[str, dict[str, Grid]]:
        """
        Get grids organized by type and probe.

        Parameters
        ----------
        grid_type : str | None
            Filter by type (density, energy). None returns all.

        Returns
        -------
        dict[str, dict[str, Grid]]
            Nested dict: {type: {probe: grid}}
        """
        if not self._grids:
            self.fetch_grids(**kwargs)

        result: dict[str, dict[str, Grid]] = {}

        for grid in self._grids or []:
            # Parse type from path or filename
            gtype = "unknown"
            if "density" in str(grid.metadata.get("path", "")):
                gtype = "density"
            elif "energy" in str(grid.metadata.get("path", "")):
                gtype = "energy"

            if grid_type and gtype != grid_type:
                continue

            if gtype not in result:
                result[gtype] = {}

            # Try to get probe from filename
            name = grid.metadata.get("name", "")
            probe = name.split("_")[0] if "_" in name else name
            result[gtype][probe] = grid

        return result

    def get_grids_by_probe(
        self,
        probe_list: list[str] | str,
        **kwargs,
    ) -> dict[str, list[Grid]]:
        """
        Get grids for specific probes.

        Parameters
        ----------
        probe_list : list[str] | str
            Probe name(s) to find.

        Returns
        -------
        dict[str, list[Grid]]
            Dict: {probe: [grids]}
        """
        if isinstance(probe_list, str):
            probe_list = [probe_list]

        if not self._grids:
            self.fetch_grids(**kwargs)

        result: dict[str, list[Grid]] = {}

        for grid in self._grids or []:
            name = grid.metadata.get("name", "")

            for probe in probe_list:
                if probe in name:
                    if probe not in result:
                        result[probe] = []
                    result[probe].append(grid)

        return result

    # =========================================================================
    # State Management
    # =========================================================================

    def set_state(self, state: ReplicaState) -> None:
        """Update replica state."""
        old_state = self.state
        self.state = state
        self._update_modified()
        log.info(f"Replica {self.name}: {old_state.name} → {state.name}")

    def _update_modified(self) -> None:
        """Update modified timestamp."""
        self.modified_at = datetime.now().isoformat()

    def to_dict(self) -> dict[str, Any]:
        """Convert to dictionary for serialization."""
        return {
            "name": self.name,
            "solvent": self.solvent,
            "path": str(self.path) if self.path else None,
            "state": self.state.name,
            "topology": self.topology,
            "coordinates": self.coordinates,
            "trajectory": self.trajectory,
            "reference": self.reference,
            "pdb": self.pdb,
            "seed": self.seed,
            "n_steps": self.n_steps,
            "settings": {
                "nanos": self.settings.nanos,
                "temperature": self.settings.temperature,
                "timestep": self.settings.timestep,
                "prod_steps": self.settings.prod_steps,
                "traj_frequency": self.settings.traj_frequency,
                "restraint_mode": self.settings.restraint_mode,
                "restraint_force": self.settings.restraint_force,
                "restraint_mask": self.settings.restraint_mask,
                "align_mask": self.settings.align_mask,
                "md_program": self.settings.md_program,
            }
            if self.settings
            else None,
            "min_folder": self.min_folder,
            "eq_folder": self.eq_folder,
            "md_folder": self.md_folder,
            "align_folder": self.align_folder,
            "density_folder": self.density_folder,
            "energy_folder": self.energy_folder,
            "md_output_template": self.md_output_template,
            "created_at": self.created_at,
            "modified_at": self.modified_at,
            "metadata": self.metadata,
            "_folders_created": self._folders_created,
            "_md_input_written": self._md_input_written,
        }

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> Replica:
        """Create Replica from dictionary."""
        # Handle state enum
        if "state" in data and isinstance(data["state"], str):
            data["state"] = ReplicaState[data["state"]]

        # Handle path
        if "path" in data and data["path"]:
            data["path"] = Path(data["path"])

        # Handle settings
        if "settings" in data and isinstance(data["settings"], dict):
            data["settings"] = MDSettings(**data["settings"])

        return cls(**data)

    def save(self, path: Path | None = None) -> None:
        """
        Save replica state to JSON file.

        Parameters
        ----------
        path : Path | None
            Output path. Defaults to {self.path}/replica.json
        """
        import json

        if path is None:
            if self.path is None:
                raise ValueError("No path specified and replica path not set")
            path = self.path / "replica.json"

        with open(path, "w") as f:
            json.dump(self.to_dict(), f, indent=2)

        log.debug(f"Saved replica to {path}")

    @classmethod
    def load(cls, path: Path) -> Replica:
        """
        Load replica from JSON file.

        Parameters
        ----------
        path : Path
            Path to replica.json file

        Returns
        -------
        Replica
            Loaded replica
        """
        import json

        with open(path) as f:
            data = json.load(f)

        # Set path from file location if not in data
        if data.get("path") is None:
            data["path"] = path.parent

        return cls.from_dict(data)

    def desc(self) -> str:
        """Return a summary description string."""
        parts = [
            f"REPLICA:{self.name}",
            f"solvent:{self.solvent}",
        ]

        if self.settings:
            parts.append(f"nanos:{self.settings.nanos}")
            if self.settings.has_restraints:
                parts.extend(
                    [
                        f"restrMode:{self.settings.restraint_mode}",
                        f"restrForce:{self.settings.restraint_force:.3f}",
                    ]
                )
            else:
                parts.append("restrMode:FREE")

        parts.extend(
            [
                f"Min:{self.is_minimization_finished()}",
                f"Eq:{self.is_equilibration_finished()}",
                f"Prod:{self.is_production_finished()}",
                f"Align:{self.is_aligned()}",
            ]
        )

        return "\t".join(parts)

    def __repr__(self) -> str:
        return f"Replica(name={self.name!r}, solvent={self.solvent!r}, state={self.state.name})"


# =========================================================================
# Helper Functions
# =========================================================================


def create_replica(
    name: str,
    solvent: str,
    base_path: Path,
    seed: int | None = None,
    settings: MDSettings | None = None,
) -> Replica:
    """
    Create a new replica with directory structure.

    Parameters
    ----------
    name : str
        Replica name
    solvent : str
        Solvent mixture name
    base_path : Path
        Base directory for replicas
    seed : int | None
        Random seed (auto-generated if None)
    settings : MDSettings | None
        MD settings

    Returns
    -------
    Replica
        New replica with created directory
    """
    import random

    if seed is None:
        seed = random.randint(1, 999999)

    replica = Replica(
        name=name,
        solvent=solvent,
        path=base_path / name,
        seed=seed,
        settings=settings or MDSettings(),
    )

    replica.create_directory()
    replica.save()

    return replica


def rename_replica_list(replicas: list[Replica]) -> None:
    """
    Set names according to solvent and count.

    E.g., 3 ETA replicas + 2 WAT → ['ETA_1', 'ETA_2', 'ETA_3', 'WAT_1', 'WAT_2']

    Parameters
    ----------
    replicas : list[Replica]
        List of replicas to rename (modified in place).
    """
    from collections import Counter

    counts: Counter[str] = Counter()

    for replica in replicas:
        counts[replica.solvent] += 1
        replica.name = f"{replica.solvent}_{counts[replica.solvent]}"


def load_replica(replica_file: Path | str | None = None) -> Replica:
    """
    Load existing replica from file.

    Parameters
    ----------
    replica_file : Path | str | None
        Path to replica.json. If None, searches current directory.

    Returns
    -------
    Replica
        Loaded replica.
    """
    if replica_file is None:
        # Search for replica files
        files = list(Path.cwd().glob("*.json")) + list(Path.cwd().glob("*.mrepl"))
        replica_files = [f for f in files if "replica" in f.name.lower()]

        if len(replica_files) > 1:
            raise ReplicaError("Multiple replica files found - specify which one")
        elif not replica_files:
            raise ReplicaError("No replica file found in current directory")

        replica_file = replica_files[0]

    return Replica.load(Path(replica_file))
