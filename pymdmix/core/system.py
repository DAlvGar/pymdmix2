"""
System Management
=================

Classes for managing molecular systems in MDMix simulations.

Defines two main classes:
    - System: Represents a macromolecule ready for solvation
    - SolvatedSystem: Represents a solvated system with topology/coordinates

Examples
--------
>>> from pymdmix.core.system import System
>>> system = System(name="protein")
>>> system.set_off("protein.off", unit_name="sys")
>>> solvated = system.solvate("ETA")
"""

from __future__ import annotations

import logging
import os
import pickle
import random
import string
import tempfile
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Any, cast

if TYPE_CHECKING:
    from pymdmix.core.solvent import Solvent
    from pymdmix.io.off_manager import OFFManager

log = logging.getLogger(__name__)


class SystemError(Exception):
    """Base exception for System errors."""

    pass


class BadFileError(SystemError):
    """Raised when a required file is not found or invalid."""

    pass


class FileLock:
    """
    File locking mechanism with context manager support.

    Cross-platform compatible file locking that doesn't rely on
    platform-specific fcntl or msvcrt.

    Parameters
    ----------
    file_name : str | Path
        File to lock
    timeout : float
        Maximum time to wait for lock
    delay : float
        Delay between lock attempts
    """

    def __init__(self, file_name: str | Path, timeout: float = 10, delay: float = 0.05):
        self.file_name = str(file_name)
        self.lockfile = f"{self.file_name}.lock"
        self.timeout = timeout
        self.delay = delay
        self.is_locked = False
        self.fd: int | None = None

    def acquire(self) -> None:
        """
        Acquire the lock.

        Raises
        ------
        SystemError
            If lock cannot be acquired within timeout
        """
        import errno

        start_time = time.time()
        while True:
            try:
                self.fd = os.open(self.lockfile, os.O_CREAT | os.O_EXCL | os.O_RDWR)
                break
            except OSError as e:
                if e.errno != errno.EEXIST:
                    raise
                if (time.time() - start_time) >= self.timeout:
                    raise SystemError("Timeout acquiring file lock")
                time.sleep(self.delay)
        self.is_locked = True

    def release(self) -> None:
        """Release the lock."""
        if self.is_locked and self.fd is not None:
            os.close(self.fd)
            try:
                os.unlink(self.lockfile)
            except OSError:
                pass
            self.is_locked = False
            self.fd = None

    def __enter__(self) -> FileLock:
        if not self.is_locked:
            self.acquire()
        return self

    def __exit__(self, exc_type: Any, exc_val: Any, exc_tb: Any) -> None:
        self.release()

    def __del__(self) -> None:
        self.release()


@dataclass
class System:
    """
    Represents a macromolecular system ready for MD preparation.

    A System can be initialized from:
    - An Amber Object File (OFF/lib) with parameterized unit
    - An Amber-ready PDB file
    - A saved System pickle file

    Attributes
    ----------
    name : str
        System identifier
    ff : list[str]
        Force field names/files used
    extra_res_list : list[str]
        Non-standard residue names to include
    off_manager : OFFManager | None
        Loaded OFF file manager
    unit_name : str | None
        Unit name in the OFF file
    ref : Any
        Reference PDB structure
    """

    name: str = ""
    ff: list[str] = field(default_factory=list)
    extra_res_list: list[str] = field(default_factory=list)
    off_manager: OFFManager | None = None
    unit_name: str | None = None
    ref: Any = None
    sys_file_path: str = ""
    _create: Any = None  # AmberCreateSystem instance (lazy loaded)

    def __post_init__(self):
        """Initialize with default name if not provided."""
        if not self.name:
            # Generate random string
            rand_str = "".join(
                random.choice(string.ascii_uppercase + string.digits) for _ in range(3)
            )
            self.name = f"mdmix_system_{rand_str}"

        if not self.sys_file_path:
            self.sys_file_path = f"{self.name}.msys"

    @classmethod
    def from_file(cls, path: str | Path) -> System:
        """
        Load a System from a saved pickle file.

        Parameters
        ----------
        path : str | Path
            Path to .msys file

        Returns
        -------
        System
            Loaded system instance
        """
        system = cls()
        system.load(path)
        return system

    @classmethod
    def from_off(
        cls,
        off_path: str | Path,
        unit_name: str | None = None,
        name: str | None = None,
        ff: list[str] | None = None,
        extra_res_list: list[str] | None = None,
    ) -> System:
        """
        Create a System from an Amber Object File.

        Parameters
        ----------
        off_path : str | Path
            Path to OFF/lib file
        unit_name : str | None
            Unit name to use (first unit if not specified)
        name : str | None
            System name
        ff : list[str] | None
            Force field names/files
        extra_res_list : list[str] | None
            Non-standard residue names

        Returns
        -------
        System
            Initialized system
        """
        system = cls(
            name=name or "",
            ff=ff or [],
            extra_res_list=extra_res_list or [],
        )
        system.set_off(off_path, unit_name)
        return system

    @classmethod
    def from_amber_pdb(
        cls,
        pdb_path: str | Path,
        name: str | None = None,
        ff: list[str] | None = None,
        extra_res_list: list[str] | None = None,
        **kwargs: Any,
    ) -> System:
        """
        Create a System from an Amber-ready PDB file.

        The PDB will be processed and converted to an OFF file.

        Parameters
        ----------
        pdb_path : str | Path
            Path to PDB file
        name : str | None
            System name
        ff : list[str] | None
            Force field names/files
        extra_res_list : list[str] | None
            Non-standard residue names
        **kwargs
            Additional arguments for PDB cleaning

        Returns
        -------
        System
            Initialized system
        """
        system = cls(
            name=name or "",
            ff=ff or [],
            extra_res_list=extra_res_list or [],
        )
        system.set_amber_pdb(pdb_path, **kwargs)
        return system

    def set_amber_pdb(
        self, pdb_path: str | Path, ff: list[str] | None = None, **kwargs: Any
    ) -> None:
        """
        Set up system from an Amber-ready PDB file.

        The PDB will be processed with cleaning steps:
        - Cap N and C terminus
        - Rename HIS residues according to protonation
        - Rename CYS to CYX if needed
        - Remove hydrogens for tLeap to add back

        Parameters
        ----------
        pdb_path : str | Path
            Path to PDB file
        ff : list[str] | None
            Additional force field files
        **kwargs
            Arguments for PDB cleaning
        """
        pdb_path = Path(pdb_path)
        if not pdb_path.exists():
            raise BadFileError(f"PDB file {pdb_path} does not exist")

        ff = ff or self.ff

        # Initialize Amber creation system
        self._init_create()

        # Create OFF from PDB
        out_off = f"{self.name}.lib"
        self._create.create_off(self.name, str(pdb_path), extra_ff=ff, **kwargs)
        self.set_off(out_off, "sys")

    def set_off(self, off_path: str | Path, unit_name: str | None = None) -> None:
        """
        Load an Amber Object File (OFF/lib).

        Parameters
        ----------
        off_path : str | Path
            Path to OFF file
        unit_name : str | None
            Unit name to use (first unit if not specified)
        """
        from pymdmix.io.off_manager import OFFManager

        off_path = Path(off_path)

        # Handle special test case
        if str(off_path).lower() == "test":
            log.info("WORKING WITH A TEST SYSTEM!")
            # Would load from test directory
            self.name = "testsystem"
            self.sys_file_path = f"{self.name}.msys"
            return

        self.off_manager = OFFManager.from_file(off_path)

        if unit_name is None:
            units = self.off_manager.get_units()
            if units:
                unit_name = units[0]

        self.unit_name = unit_name
        self._create_ref()

    def solvate(
        self,
        solvent: str | Solvent,
        suffix: str | None = None,
        tmp: bool | str = True,
    ) -> SolvatedSystem:
        """
        Solvate the system with a solvent mixture.

        Parameters
        ----------
        solvent : str | Solvent
            Solvent name or instance
        suffix : str | None
            Suffix for output files
        tmp : bool | str
            If True, work in temp directory; if str, use as path

        Returns
        -------
        SolvatedSystem
            Solvated system instance

        Raises
        ------
        SystemError
            If no OFF file is assigned
        """
        if self.off_manager is None:
            raise SystemError("Cannot solvate without an Amber Object File assigned")

        # Get solvent instance
        from pymdmix.core.solvent import SolventLibrary

        if isinstance(solvent, str):
            library = SolventLibrary()
            solvent_obj = library.get(solvent)
            if solvent_obj is None:
                raise SystemError(f"Unknown solvent: {solvent}")
            solvent_name = solvent
        else:
            solvent_obj = solvent
            solvent_name = solvent.name

        log.info(f"Solvating {self.name} with solvent mixture {solvent_name}")

        suffix = suffix or solvent_name
        name = f"{self.name}_{suffix}"
        prmtop = f"{name}.prmtop"
        prmcrd = f"{name}.prmcrd"

        # Determine working path
        if tmp:
            if isinstance(tmp, str):
                path = Path(tmp)
            else:
                path = Path(tempfile.mkdtemp())
            prmtop = str(path / prmtop)
            prmcrd = str(path / prmcrd)

        # Initialize creation system and solvate
        self._init_create()
        self._create.solvate_organic(self.unit_name, solvent=solvent_obj)
        self._create.save_amber_parm(self.unit_name, prmtop, prmcrd)
        self._clean_create()

        return SolvatedSystem(
            name=name,
            top=prmtop,
            crd=prmcrd,
            solvent=solvent_name,
            ref=self.ref,
            ff=self.ff,
            extra_res_list=self.extra_res_list,
        )

    def write_off(self, fname: str | Path) -> Path | None:
        """
        Write Object File to disk.

        Parameters
        ----------
        fname : str | Path
            Output file name

        Returns
        -------
        Path | None
            Path to written file, or None if no OFF loaded
        """
        if self.off_manager:
            return cast(Path, self.off_manager.write(fname))
        return None

    def _init_create(self) -> None:
        """Initialize AmberCreateSystem if not already loaded."""
        if self._create is None:
            # Import here to avoid circular dependency
            try:
                from pymdmix.engines.amber import AmberCreateSystem

                self._create = AmberCreateSystem(ff_list=self.ff, informative=False)

                if self.off_manager:
                    tmp_off = self.off_manager.write_tmp()
                    self._create.load_off(str(tmp_off))
                    self.off_manager.clean_tmp()
            except ImportError:
                log.warning("AmberCreateSystem not available - solvation will not work")
                self._create = None

    def _clean_create(self) -> None:
        """Clean up AmberCreateSystem resources."""
        if self._create is not None:
            if hasattr(self._create, "leap"):
                self._create.leap.close()
            self._create = None

    def _create_ref(self) -> None:
        """Create reference PDB from OFF file."""
        if not os.environ.get("AMBERHOME"):
            log.warning("AMBERHOME not set; skipping reference PDB creation")
            return

        # Initialize creation system
        self._init_create()

        if self._create is None:
            log.warning("Cannot create reference PDB - AmberCreateSystem not available")
            return

        try:
            # Save temporary topology and coordinates
            self._create.save_amber_parm(self.unit_name, "tmp.top", "tmp.crd")
            self._create.ambpdb("tmp.top", "tmp.crd", "tmp.pdb")

            # Load reference PDB
            try:
                import MDAnalysis as mda

                self.ref = mda.Universe("tmp.pdb")
            except ImportError:
                # Fallback: just store path
                self.ref = Path("tmp.pdb").read_text()

            # Clean up temp files
            for ext in ("top", "crd", "pdb"):
                Path(f"tmp.{ext}").unlink(missing_ok=True)
        except Exception as e:
            log.warning(f"Could not create reference PDB automatically: {e}")
        finally:
            self._clean_create()

    def copy(self) -> System:
        """
        Create a deep copy of the system.

        Returns
        -------
        System
            Copy of the system
        """
        import copy

        return copy.deepcopy(self)

    def write(self, outfile: str | Path | None = None) -> None:
        """
        Save system to a pickle file.

        Parameters
        ----------
        outfile : str | Path | None
            Output file path (defaults to sys_file_path)
        """
        outfile = Path(outfile or self.sys_file_path)

        with FileLock(str(outfile)) as _:
            # Prepare data for pickling
            data = {
                "name": self.name,
                "ff": self.ff,
                "extra_res_list": self.extra_res_list,
                "off_content": self.off_manager.content if self.off_manager else None,
                "unit_name": self.unit_name,
                "ref": self.ref,
                "sys_file_path": str(outfile.absolute()),
            }

            with open(outfile, "wb") as f:
                pickle.dump(data, f)

        self.sys_file_path = str(outfile.absolute())

    def load(self, sysfile: str | Path | None = None) -> None:
        """
        Load system from a pickle file.

        Parameters
        ----------
        sysfile : str | Path | None
            File to load (defaults to sys_file_path)
        """
        from pymdmix.io.off_manager import OFFManager

        sysfile = Path(sysfile or self.sys_file_path)
        if not sysfile.exists():
            raise BadFileError(f"File {sysfile} not found")

        with FileLock(str(sysfile)) as _:
            with open(sysfile, "rb") as f:
                data = pickle.load(f)

        self.name = data["name"]
        self.ff = data.get("ff", [])
        self.extra_res_list = data.get("extra_res_list", [])
        self.unit_name = data.get("unit_name")
        self.ref = data.get("ref")
        self.sys_file_path = data.get("sys_file_path", str(sysfile))

        # Restore OFF manager
        off_content = data.get("off_content")
        if off_content:
            self.off_manager = OFFManager.from_string(off_content)

    def __repr__(self) -> str:
        return f"System({self.name!r})"

    def __str__(self) -> str:
        lines = [f"SYSTEM: {self.name}"]
        if self.ff:
            lines.append(f"  FF: {self.ff}")
        if self.extra_res_list:
            lines.append(f"  ExtraRes: {self.extra_res_list}")
        return "\n".join(lines)


@dataclass
class SolvatedSystem(System):
    """
    Represents a solvated system with topology and coordinates.

    Subclass of System containing PRMTOP/PRMCRD files.

    Attributes
    ----------
    top : str | None
        PRMTOP file content
    crd : str | None
        PRMCRD file content
    pdb : str | None
        PDB file content (generated from top/crd)
    solvent : str | None
        Solvent used for solvation
    """

    top: str | None = None
    crd: str | None = None
    pdb: str | None = None
    solvent: str | None = None
    _tmp_top: Path | None = None
    _tmp_crd: Path | None = None
    _tmp_pdb: Path | None = None

    def __post_init__(self):
        """Initialize parent and set topology files if paths provided."""
        super().__post_init__()

        # If top/crd are paths, load them
        if self.top and Path(self.top).exists():
            self.set_top_crd(self.top, self.crd)

    @classmethod
    def from_top_crd(
        cls,
        name: str,
        top: str | Path,
        crd: str | Path,
        solvent: str | None = None,
        ref: Any = None,
        **kwargs: Any,
    ) -> SolvatedSystem:
        """
        Create a SolvatedSystem from PRMTOP and PRMCRD files.

        Parameters
        ----------
        name : str
            System name
        top : str | Path
            Path to PRMTOP file
        crd : str | Path
            Path to PRMCRD file
        solvent : str | None
            Solvent name
        ref : Any
            Reference structure
        **kwargs
            Additional arguments passed to constructor

        Returns
        -------
        SolvatedSystem
            Initialized solvated system
        """
        system = cls(name=name, solvent=solvent, ref=ref, **kwargs)
        system.set_top_crd(top, crd)
        return system

    def set_top_crd(self, top: str | Path, crd: str | Path) -> None:
        """
        Set PRMTOP and PRMCRD for the system.

        Parameters
        ----------
        top : str | Path
            Path to PRMTOP file
        crd : str | Path
            Path to PRMCRD file
        """
        top_path = Path(top)
        crd_path = Path(crd)

        if not top_path.exists():
            raise BadFileError(f"File {top_path} not found")
        if not crd_path.exists():
            raise BadFileError(f"File {crd_path} not found")

        self.top = top_path.read_text()
        self.crd = crd_path.read_text()
        self.set_pdb_from_top_crd()

        # Set up reference if not provided
        if self.ref is None:
            self.ref = self.get_solvated_pdb_solute()

    def save_top_crd(self, prefix: str | Path) -> bool:
        """
        Save topology and coordinates to disk files.

        Parameters
        ----------
        prefix : str | Path
            Prefix for output files (.prmtop and .prmcrd extensions added)

        Returns
        -------
        bool
            True if successful
        """
        if self.top and self.crd:
            prefix = Path(prefix)
            prefix.with_suffix(".prmtop").write_text(self.top)
            prefix.with_suffix(".prmcrd").write_text(self.crd)
            return prefix.with_suffix(".prmtop").exists() and prefix.with_suffix(".prmcrd").exists()
        return False

    def save_pdb(self, fname: str | Path) -> bool:
        """
        Save PDB to a file.

        Parameters
        ----------
        fname : str | Path
            Output file name

        Returns
        -------
        bool
            True if successful
        """
        if self.pdb:
            Path(fname).write_text(self.pdb)
            return Path(fname).exists()
        return False

    def get_solvated_pdb(self) -> Any:
        """
        Get a solvated structure object from the PDB content.

        Returns
        -------
        Any
            SolvatedPDB instance
        """
        try:
            from pymdmix.core.structure import load_structure

            pdb_path = self.get_tmp_pdb_file()
            if pdb_path is None:
                return None
            result = load_structure(pdb_path)
            self.clean_tmp()
            return result
        except ImportError:
            log.warning("Structure loader not available")
            return None

    def get_solvated_pdb_solute(self) -> Any:
        """
        Get the solute (non-solvent) portion of the solvated system.

        Returns
        -------
        Any
            Solute structure
        """
        solvated = self.get_solvated_pdb()
        if solvated and hasattr(solvated, "get_solute"):
            return solvated.get_solute()
        return None

    def get_tmp_top_crd_files(self) -> tuple[Path, Path] | None:
        """
        Save topology and coordinates to temporary files.

        Returns
        -------
        tuple[Path, Path] | None
            (prmtop_path, prmcrd_path) or None if failed
        """
        if not self.top or not self.crd:
            return None

        tmp = tempfile.mktemp()
        self._tmp_top = Path(f"{tmp}.prmtop")
        self._tmp_crd = Path(f"{tmp}.prmcrd")

        if self.save_top_crd(tmp):
            return self._tmp_top, self._tmp_crd
        return None

    def get_tmp_pdb_file(self) -> Path | None:
        """
        Save PDB to a temporary file.

        Returns
        -------
        Path | None
            Path to temporary PDB file
        """
        if self.pdb:
            fd, tmp = tempfile.mkstemp(suffix=".pdb")
            self._tmp_pdb = Path(tmp)
            with open(fd, "w") as f:
                f.write(self.pdb)
            return self._tmp_pdb
        return None

    def clean_tmp(self) -> None:
        """Remove temporary files."""
        for tmp_file in (self._tmp_top, self._tmp_crd, self._tmp_pdb):
            if tmp_file and tmp_file.exists():
                tmp_file.unlink()
        self._tmp_top = None
        self._tmp_crd = None
        self._tmp_pdb = None

    def set_pdb_from_top_crd(self) -> None:
        """Generate PDB content from topology and coordinates."""
        tmp_files = self.get_tmp_top_crd_files()
        if tmp_files is None:
            return

        tmp_top, tmp_crd = tmp_files
        fd, tmp_pdb_path = tempfile.mkstemp(suffix=".pdb")
        os.close(fd)
        self._tmp_pdb = Path(tmp_pdb_path)

        try:
            from pymdmix.engines.amber import AmberCreateSystem

            create = AmberCreateSystem(informative=False)
            create.ambpdb(str(tmp_top), str(tmp_crd), str(self._tmp_pdb))
            time.sleep(1)  # Allow ambpdb to complete

            if self._tmp_pdb.exists():
                self.pdb = self._tmp_pdb.read_text()
        except ImportError:
            log.warning("AmberCreateSystem not available - cannot generate PDB from topology")
        finally:
            self.clean_tmp()

    def solvate(
        self,
        solvent: str | Solvent,
        suffix: str | None = None,
        tmp: bool | str = True,
    ) -> SolvatedSystem:
        """Override to prevent re-solvation of already solvated system."""
        raise SystemError("Cannot solvate an already solvated system")

    def __repr__(self) -> str:
        return f"SolvatedSystem({self.name!r}, solvent={self.solvent!r})"


def load_system(systemfile: str | Path | None = None) -> System:
    """
    Load an existing system from file.

    If no path is given, searches for a .msys file in the current directory.

    Parameters
    ----------
    systemfile : str | Path | None
        Path to system file

    Returns
    -------
    System
        Loaded system

    Raises
    ------
    SystemError
        If multiple or no system files found
    """
    import glob

    if systemfile is None:
        files = glob.glob("*.msys")
        if len(files) > 1:
            raise SystemError(
                "Multiple system files in current folder. Please specify which file to load."
            )
        if not files:
            raise SystemError("No .msys file found in current folder and no path given.")
        systemfile = files[0]

    return System.from_file(systemfile)
