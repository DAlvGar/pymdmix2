"""
End-to-end tests that exercise real AmberTools binaries (tleap, cpptraj).

These tests are tagged with ``@pytest.mark.ambertools`` (and
``@pytest.mark.cpptraj`` where cpptraj is specifically needed) and are
**automatically skipped** when AmberTools is not present in the environment.

Running the tests
-----------------
Inside a Docker container (recommended — includes AmberTools)::

    ./scripts/run_ambertools_tests.sh

Directly, if AmberTools is already installed::

    pytest -m ambertools -v tests/test_ambertools_e2e.py

Coverage
--------
1. AmberTools availability smoke-tests
2. :class:`~pymdmix.engines.amber.LeapSession` interactive API
3. Solvation workflow: ``generate_leap_script`` + ``run_leap`` + ``solvate_structure``
4. ``CpptrajDensityAction`` on unsolvated and solvated systems
5. Full pipeline: solvate → cpptraj density → :class:`~pymdmix.core.grid.Grid` → hotspots
6. Full CLI project creation: project config → solvated replicas ready for Amber simulation
"""

from __future__ import annotations

import os
import shutil
import subprocess
from importlib.resources import files
from pathlib import Path

import pytest

# ---------------------------------------------------------------------------
# Test data paths (tests/data/)
# ---------------------------------------------------------------------------

_TEST_DATA_DIR = Path(__file__).parent / "data"
_PEP_DIR = _TEST_DATA_DIR / "pep"
_SOLVENTS_DIR = Path(str(files("pymdmix").joinpath("data/solvents")))

_PEP_PDB = _PEP_DIR / "pep.pdb"
_PEP_OFF = _PEP_DIR / "pep.off"
_PEP_PRMTOP = _PEP_DIR / "pep.prmtop"
_PEP_PRMCRD = _PEP_DIR / "pep.prmcrd"


# ---------------------------------------------------------------------------
# Helper: locate binaries
# ---------------------------------------------------------------------------


def _find_tleap() -> str | None:
    """Return path to tleap binary, or None if unavailable."""
    amber_home = os.environ.get("AMBERHOME")
    if amber_home:
        for sub in ("bin", "exe"):
            candidate = Path(amber_home) / sub / "tleap"
            if candidate.exists():
                return str(candidate)
    return shutil.which("tleap")


def _find_cpptraj() -> str | None:
    """Return path to cpptraj binary, or None if unavailable."""
    return os.environ.get("AMBER_PTRAJ") or shutil.which("cpptraj")


# ---------------------------------------------------------------------------
# Session-scoped solvated-pep fixture (built once, shared across groups 4 & 5)
# ---------------------------------------------------------------------------


@pytest.fixture(scope="session")
def solvated_pep(tmp_path_factory):
    """
    Solvate pep.off with ETAWAT20 (ETA solvent) once per test session.

    Returns ``(prmtop_path, inpcrd_path)`` or calls ``pytest.skip()`` if
    tleap is unavailable or solvation fails.
    """
    from pymdmix.core.solvent import SolventLibrary
    from pymdmix.setup.solvate import BoxConfig, IonConfig, generate_leap_script, run_leap

    tleap = _find_tleap()
    if not tleap:
        pytest.skip("tleap not available — skipping solvation-dependent tests")

    out_dir = tmp_path_factory.mktemp("solvated_pep")
    library = SolventLibrary()
    solvent = library.get("ETA")

    # Smaller buffer so the test runs faster
    box_config = BoxConfig(shape="box", buffer=8.0)
    ion_config = IonConfig(neutralize=True, add_ions=False)

    script = generate_leap_script(
        input_path=_PEP_OFF,
        solvent=solvent,
        output_prefix=str(out_dir / "solvated"),
        unit_name="pep",
        box_config=box_config,
        ion_config=ion_config,
    )

    success, log = run_leap(script, leap_exe=tleap, work_dir=out_dir)

    if not success:
        pytest.skip(f"tleap solvation failed (first 500 chars of log):\n{log[:500]}")

    prmtop = out_dir / "solvated.prmtop"
    inpcrd = out_dir / "solvated.inpcrd"

    if not prmtop.exists() or not inpcrd.exists():
        pytest.skip("tleap did not produce expected output files — check log for errors")

    return prmtop, inpcrd


# ===========================================================================
# Group 1 — AmberTools availability smoke-tests
# ===========================================================================


@pytest.mark.ambertools
class TestAmberToolsAvailability:
    """Basic smoke-tests that confirm the AmberTools binaries are functional."""

    def test_tleap_found(self):
        """tleap is reachable via AMBERHOME or PATH."""
        assert _find_tleap() is not None, (
            "tleap not found. Set AMBERHOME or add tleap to PATH."
        )

    def test_tleap_runs_quit(self, tmp_path):
        """tleap accepts a 'quit' command and exits successfully."""
        tleap = _find_tleap()
        script = tmp_path / "quit.in"
        script.write_text("quit\n")
        result = subprocess.run(
            [tleap, "-f", str(script)],
            capture_output=True,
            text=True,
            timeout=30,
        )
        assert result.returncode == 0, f"tleap exited with {result.returncode}: {result.stderr}"

    @pytest.mark.cpptraj
    def test_cpptraj_found(self):
        """cpptraj is reachable via AMBER_PTRAJ or PATH."""
        assert _find_cpptraj() is not None, (
            "cpptraj not found. Set AMBER_PTRAJ or add cpptraj to PATH."
        )

    @pytest.mark.cpptraj
    def test_cpptraj_version(self):
        """cpptraj --version exits without error."""
        cpptraj = _find_cpptraj()
        result = subprocess.run(
            [cpptraj, "--version"],
            capture_output=True,
            text=True,
            timeout=15,
        )
        combined = result.stdout + result.stderr
        # Either exit code 0 or "CPPTRAJ" appears in the output
        assert result.returncode == 0 or "CPPTRAJ" in combined.upper(), (
            f"Unexpected cpptraj output: {combined[:300]}"
        )

    def test_bundled_test_data_present(self):
        """pep.pdb, pep.off, pep.prmtop, and pep.prmcrd are bundled with the package."""
        for path in (_PEP_PDB, _PEP_OFF, _PEP_PRMTOP, _PEP_PRMCRD):
            assert path.exists(), f"Missing bundled test file: {path}"
            assert path.stat().st_size > 0, f"Bundled test file is empty: {path}"

    def test_eta_solvent_off_present(self):
        """ETAWAT20.off (ETA solvent box) is bundled with the package."""
        eta_off = _SOLVENTS_DIR / "ETAWAT20.off"
        assert eta_off.exists(), f"Missing ETA OFF file: {eta_off}"


# ===========================================================================
# Group 2 — LeapSession interactive API
# ===========================================================================


@pytest.mark.ambertools
class TestLeapSession:
    """Tests for :class:`~pymdmix.engines.amber.LeapSession` (interactive tleap)."""

    def test_context_manager(self, tmp_path):
        """LeapSession opens and closes cleanly as a context manager."""
        from pymdmix.engines.amber import LeapSession

        with LeapSession(cwd=tmp_path):
            pass  # just open and close — no exception means success

    def test_source_protein_ff(self, tmp_path):
        """LeapSession can source leaprc.protein.ff14SB without FATAL errors."""
        from pymdmix.engines.amber import LeapSession

        with LeapSession(cwd=tmp_path) as leap:
            output = leap.source("leaprc.protein.ff14SB")
        assert not any("FATAL" in line for line in output), (
            f"FATAL error sourcing leaprc.protein.ff14SB: {output}"
        )

    def test_load_pdb(self, tmp_path):
        """LeapSession can load pep.pdb without FATAL errors."""
        from pymdmix.engines.amber import LeapSession

        with LeapSession(extra_ff=["leaprc.protein.ff14SB"], cwd=tmp_path) as leap:
            output = leap.load_pdb(_PEP_PDB, unit_name="pep")
        assert not any("FATAL" in line for line in output), (
            f"FATAL error loading pep.pdb: {output}"
        )

    def test_load_off(self, tmp_path):
        """LeapSession can load pep.off without FATAL errors."""
        from pymdmix.engines.amber import LeapSession

        with LeapSession(extra_ff=["leaprc.protein.ff14SB"], cwd=tmp_path) as leap:
            output = leap.load_off(_PEP_OFF)
        assert not any("FATAL" in line for line in output), (
            f"FATAL error loading pep.off: {output}"
        )

    def test_charge(self, tmp_path):
        """LeapSession.charge() returns a float for the loaded peptide unit."""
        from pymdmix.engines.amber import LeapSession

        with LeapSession(extra_ff=["leaprc.protein.ff14SB"], cwd=tmp_path) as leap:
            leap.load_off(_PEP_OFF)
            charge = leap.charge("pep")
        assert isinstance(charge, float), f"charge() should return float, got {type(charge)}"

    def test_save_amber_parm_from_pdb(self, tmp_path):
        """LeapSession can save topology and coordinates from pep.pdb."""
        from pymdmix.engines.amber import LeapSession

        top = tmp_path / "out.prmtop"
        crd = tmp_path / "out.inpcrd"

        with LeapSession(extra_ff=["leaprc.protein.ff14SB"], cwd=tmp_path) as leap:
            leap.load_pdb(_PEP_PDB, unit_name="pep")
            success = leap.save_amber_parm("pep", top, crd)

        assert success, "save_amber_parm() should return True on success"
        assert top.exists() and top.stat().st_size > 0, "Topology file is missing or empty"
        assert crd.exists() and crd.stat().st_size > 0, "Coordinate file is missing or empty"

    def test_save_amber_parm_from_off(self, tmp_path):
        """LeapSession can save topology and coordinates from pep.off."""
        from pymdmix.engines.amber import LeapSession

        top = tmp_path / "pep.prmtop"
        crd = tmp_path / "pep.inpcrd"

        with LeapSession(extra_ff=["leaprc.protein.ff14SB"], cwd=tmp_path) as leap:
            leap.load_off(_PEP_OFF)
            success = leap.save_amber_parm("pep", top, crd)

        assert success, "save_amber_parm() should return True when loading from OFF"
        assert top.exists() and top.stat().st_size > 0


# ===========================================================================
# Group 3 — Solvation workflow
# ===========================================================================


@pytest.mark.ambertools
class TestSolvation:
    """Tests for the LEaP-based solvation pipeline."""

    def test_generate_script_pdb_input(self):
        """generate_leap_script() produces a valid script for PDB input."""
        from pymdmix.core.solvent import SolventLibrary
        from pymdmix.setup.solvate import generate_leap_script

        solvent = SolventLibrary().get("ETA")
        script = generate_leap_script(
            input_path=_PEP_PDB,
            solvent=solvent,
            output_prefix="solvated",
        )
        assert "leaprc.protein.ff14SB" in script
        assert "loadpdb" in script.lower()
        assert "saveamberparm" in script.lower()

    def test_generate_script_off_input(self):
        """generate_leap_script() uses loadoff + copy for OFF input."""
        from pymdmix.core.solvent import SolventLibrary
        from pymdmix.setup.solvate import generate_leap_script

        solvent = SolventLibrary().get("ETA")
        script = generate_leap_script(
            input_path=_PEP_OFF,
            solvent=solvent,
            output_prefix="solvated",
            unit_name="pep",
        )
        assert "loadoff" in script.lower()
        assert "sys = copy pep" in script

    def test_run_leap_solvation_off(self, tmp_path):
        """Full tleap solvation: pep.off + ETAWAT20 → prmtop + inpcrd."""
        from pymdmix.core.solvent import SolventLibrary
        from pymdmix.setup.solvate import BoxConfig, IonConfig, generate_leap_script, run_leap

        solvent = SolventLibrary().get("ETA")
        box_config = BoxConfig(shape="box", buffer=8.0)
        ion_config = IonConfig(neutralize=True, add_ions=False)

        script = generate_leap_script(
            input_path=_PEP_OFF,
            solvent=solvent,
            output_prefix=str(tmp_path / "solvated"),
            unit_name="pep",
            box_config=box_config,
            ion_config=ion_config,
        )

        tleap = _find_tleap()
        success, log = run_leap(script, leap_exe=tleap, work_dir=tmp_path)

        assert success, f"tleap solvation failed. Log (first 500 chars):\n{log[:500]}"
        prmtop = tmp_path / "solvated.prmtop"
        inpcrd = tmp_path / "solvated.inpcrd"
        assert prmtop.exists() and prmtop.stat().st_size > 0, "prmtop not created"
        assert inpcrd.exists() and inpcrd.stat().st_size > 0, "inpcrd not created"

    def test_solvate_structure_api(self, tmp_path):
        """High-level solvate_structure() API produces topology and coordinates."""
        from pymdmix.core.solvent import SolventLibrary
        from pymdmix.setup.solvate import SolvationOptions, solvate_structure

        solvent = SolventLibrary().get("ETA")
        options = SolvationOptions(
            box_buffer=8.0,
            box_shape="box",
            neutralize=True,
            ion_concentration=0.0,
        )

        result = solvate_structure(
            structure=_PEP_OFF,
            solvent=solvent,
            output_dir=tmp_path,
            output_prefix="solvated",
            unit_name="pep",
            options=options,
            leap_exe=_find_tleap() or "tleap",
        )

        assert result.success, f"solvate_structure() failed: {result.error}"
        assert result.topology is not None and result.topology.exists()
        assert result.coordinates is not None and result.coordinates.exists()
        assert result.topology.stat().st_size > 0
        assert result.coordinates.stat().st_size > 0


# ===========================================================================
# Group 4 — CpptrajDensityAction
# ===========================================================================


@pytest.mark.ambertools
@pytest.mark.cpptraj
class TestCpptrajDensity:
    """Tests for :class:`~pymdmix.analysis.density.CpptrajDensityAction`."""

    def test_density_backbone_atoms(self, tmp_path):
        """
        CpptrajDensityAction runs on pep.prmtop + pep.prmcrd using the
        backbone CA atoms as a probe (no solvent required).
        """
        from pymdmix.analysis.density import CpptrajDensityAction

        action = CpptrajDensityAction()
        result = action.run(
            topology=_PEP_PRMTOP,
            trajectory_pattern=[str(_PEP_PRMCRD)],
            probe_masks={"CA": "@CA"},
            grid_dimensions=(20, 20, 20),
            grid_origin=(0.0, 0.0, 0.0),
            grid_spacing=1.0,
            output_dir=tmp_path,
        )

        assert result.success, f"CpptrajDensityAction failed: {result.error}"
        dx_file = tmp_path / "CA.dx"
        assert dx_file.exists(), "Density DX file was not produced"
        assert dx_file.stat().st_size > 0, "Density DX file is empty"

    def test_density_grid_can_be_loaded(self, tmp_path):
        """DX output from cpptraj can be loaded as a pyMDMix Grid object."""
        from pymdmix.analysis.density import CpptrajDensityAction
        from pymdmix.core.grid import Grid

        action = CpptrajDensityAction()
        result = action.run(
            topology=_PEP_PRMTOP,
            trajectory_pattern=[str(_PEP_PRMCRD)],
            probe_masks={"CA": "@CA"},
            grid_dimensions=(20, 20, 20),
            grid_origin=(0.0, 0.0, 0.0),
            grid_spacing=1.0,
            output_dir=tmp_path,
        )

        assert result.success, f"CpptrajDensityAction failed: {result.error}"
        grid = Grid.read_dx(tmp_path / "CA.dx")
        assert grid is not None
        assert grid.data.shape == (20, 20, 20)

    def test_density_solvated_system(self, tmp_path, solvated_pep):
        """
        CpptrajDensityAction computes density for ETA probes on the solvated
        peptide system produced by the session-scoped solvation fixture.
        """
        from pymdmix.analysis.density import CpptrajDensityAction
        from pymdmix.core.solvent import SolventLibrary

        prmtop, inpcrd = solvated_pep
        solvent = SolventLibrary().get("ETA")

        # Build Amber masks from the first probe in the ETA definition
        probe = solvent.probes[0]
        atoms_joined = ",".join(probe.atoms)
        probe_masks = {probe.name: f":{probe.residue}@{atoms_joined}"}

        action = CpptrajDensityAction()
        result = action.run(
            topology=prmtop,
            trajectory_pattern=[str(inpcrd)],
            probe_masks=probe_masks,
            grid_dimensions=(20, 20, 20),
            grid_origin=(-10.0, -10.0, -10.0),
            grid_spacing=1.0,
            output_dir=tmp_path,
        )

        assert result.success, f"CpptrajDensityAction on solvated system failed: {result.error}"
        dx_file = tmp_path / f"{probe.name}.dx"
        assert dx_file.exists(), f"Density file {dx_file.name} not produced"


# ===========================================================================
# Group 5 — Full pipeline: solvate → cpptraj density → Grid → hotspots
# ===========================================================================


@pytest.mark.ambertools
@pytest.mark.cpptraj
class TestFullSolvationDensityWorkflow:
    """
    Mimics a realistic pyMDMix workflow:
    pep.off → solvate with ETA → prmtop/inpcrd → cpptraj density → Grid → hotspots
    """

    def test_full_workflow(self, tmp_path, solvated_pep):
        """End-to-end: solvated prmtop → density → Grid → hotspot detection."""
        from pymdmix.analysis.density import CpptrajDensityAction
        from pymdmix.analysis.hotspots import HotspotAction
        from pymdmix.core.grid import Grid
        from pymdmix.core.solvent import SolventLibrary

        prmtop, inpcrd = solvated_pep
        solvent = SolventLibrary().get("ETA")

        # Use the OH probe (H-bond donor/acceptor) — most scientifically relevant
        oh_probe = next((p for p in solvent.probes if p.name == "OH"), solvent.probes[0])
        atoms_joined = ",".join(oh_probe.atoms)
        probe_masks = {oh_probe.name: f":{oh_probe.residue}@{atoms_joined}"}

        # ---- Step 1: compute density with cpptraj ----
        density_dir = tmp_path / "density"
        density_action = CpptrajDensityAction()
        density_result = density_action.run(
            topology=prmtop,
            trajectory_pattern=[str(inpcrd)],
            probe_masks=probe_masks,
            grid_dimensions=(30, 30, 30),
            grid_origin=(-15.0, -15.0, -15.0),
            grid_spacing=1.0,
            output_dir=density_dir,
        )

        assert density_result.success, f"Density step failed: {density_result.error}"
        dx_file = density_dir / f"{oh_probe.name}.dx"
        assert dx_file.exists(), "DX density file not produced"

        # ---- Step 2: load Grid ----
        grid = Grid.read_dx(dx_file)
        assert grid is not None
        assert grid.data.shape == (30, 30, 30)

        # ---- Step 3: hotspot detection ----
        hotspot_dir = tmp_path / "hotspots"
        hotspot_action = HotspotAction()
        hotspot_result = hotspot_action.run(
            grids={oh_probe.name: grid},
            output_dir=hotspot_dir,
        )

        assert hotspot_result.success, f"Hotspot step failed: {hotspot_result.error}"


# ===========================================================================
# Group 6 — Full CLI project creation workflow
# ===========================================================================


@pytest.mark.ambertools
class TestProjectCreationWorkflow:
    """
    Full end-to-end CLI test: project config → solvated replicas ready for Amber simulation.

    Starting from a ``.cfg`` project config file the test drives the
    ``pymdmix create project`` command and verifies the complete output
    directory tree including:

    * Solvated topology / coordinates (tleap)
    * Amber MD input files (min.in, eq1.in, eq2.in, prod.in) with the
      correct restraint settings for heavy-atom (HA) mode
    * COMMANDS.sh submission script
    * Replicas serialised as ``replica.json`` with state ``READY``
    * Solvated PDB contains more atoms than the dry input system

    The project is created once per class via the ``_project_dir`` fixture so
    tleap only runs once even though multiple test methods inspect the output.
    """

    @pytest.fixture(scope="class")
    def _project_dir(self, tmp_path_factory):
        """
        Create the full pymdmix project once and share it across all tests
        in this class.  tleap is only invoked once per test class.
        """
        import textwrap

        from click.testing import CliRunner

        from pymdmix.cli import cli as pymdmix_cli

        tmp_path = tmp_path_factory.mktemp("project_creation")
        config_content = textwrap.dedent(f"""\
            [SYSTEM]
            NAME = pep
            OFF = {_PEP_OFF}
            UNAME = pep

            [MDSETTINGS]
            SOLVENTS = ETA
            NREPL = 2
            NANOS = 100
            TEMP = 300
            RESTR = HA
            FORCE = 10.0
        """)
        config_file = tmp_path / "project.cfg"
        config_file.write_text(config_content)

        project_dir = tmp_path / "myproject"

        runner = CliRunner()
        result = runner.invoke(
            pymdmix_cli,
            [
                "create", "project",
                "-n", "myproject",
                "-f", str(config_file),
                "-d", str(project_dir),
            ],
        )

        # Surface any exception immediately so the fixture fails clearly
        if result.exception:
            import traceback

            tb = "".join(
                traceback.format_exception(
                    type(result.exception),
                    result.exception,
                    result.exception.__traceback__,
                )
            )
            pytest.fail(
                f"CLI raised an exception (exit code {result.exit_code}):\n{tb}\n"
                f"--- stdout ---\n{result.output}"
            )

        if result.exit_code != 0:
            pytest.fail(
                f"pymdmix create project exited with code {result.exit_code}.\n"
                f"--- output ---\n{result.output}"
            )

        return project_dir

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _count_pdb_atoms(pdb_path: Path) -> int:
        """Return the number of ATOM / HETATM records in *pdb_path*."""
        return sum(
            1
            for line in pdb_path.read_text().splitlines()
            if line.startswith(("ATOM  ", "HETATM"))
        )

    # ------------------------------------------------------------------
    # Test methods
    # ------------------------------------------------------------------

    def test_create_project_eta_2replicas_100ns_ha_restraints(self, _project_dir):
        """
        Validate the complete output of ``pymdmix create project``.

        Configuration:
        - Input system : pep.off (bundled test peptide, pre-parameterised)
        - Solvent      : ETA (20% ethanol)
        - Replicas     : 2
        - Length       : 100 ns
        - Restraints   : HA (soft restraints on all non-hydrogen atoms,
                         10 kcal/mol·Å²)

        Assertions cover the complete folder tree, solvated topologies,
        all Amber input files, restraint parameters, and replica state.
        """
        import json

        project_dir = _project_dir

        # ---- 1. Project directory structure ----------------------------------
        assert (project_dir / "project.json").exists(), (
            "project.json not created"
        )
        for subdir in ("replicas", "systems", "input"):
            assert (project_dir / subdir).is_dir(), f"Missing project subdir: {subdir}"

        # Input config should have been copied into the project
        assert (project_dir / "input" / "project.cfg").exists(), (
            "Input config not copied to project/input/"
        )

        # ---- 2. Solvated system files ----------------------------------------
        eta_systems_dir = project_dir / "systems" / "ETA"
        prmtop = eta_systems_dir / "pep_ETA.prmtop"
        inpcrd = eta_systems_dir / "pep_ETA.inpcrd"
        assert prmtop.exists() and prmtop.stat().st_size > 0, (
            f"Solvated topology missing or empty: {prmtop}"
        )
        assert inpcrd.exists() and inpcrd.stat().st_size > 0, (
            f"Solvated coordinates missing or empty: {inpcrd}"
        )

        # ---- 3. Two replica directories, each fully populated ----------------
        for i in (1, 2):
            rep_dir = project_dir / "replicas" / f"pep_ETA_{i}"
            assert rep_dir.is_dir(), f"Replica directory missing: {rep_dir}"

            # Topology / coordinates copied from systems/
            for fname in ("pep_ETA.prmtop", "pep_ETA.inpcrd"):
                fpath = rep_dir / fname
                assert fpath.exists() and fpath.stat().st_size > 0, (
                    f"{fname} missing or empty in replica {i}"
                )

            # All Amber input files present
            expected_inputs = {
                "min/min.in": "minimization input",
                "eq/eq1.in": "NVT heating input",
                "eq/eq2.in": "NPT equilibration input",
                "md/prod.in": "production MD input",
                "COMMANDS.sh": "submission script",
            }
            for rel, label in expected_inputs.items():
                fpath = rep_dir / rel
                assert fpath.exists() and fpath.stat().st_size > 0, (
                    f"Replica {i}: {label} ({rel}) missing or empty"
                )

            # ---- 4. Restraint parameters in input files ----------------------
            # Minimization: HA restraints enabled (ntr=1) with correct mask
            min_in = (rep_dir / "min" / "min.in").read_text()
            assert "ntr=1," in min_in.replace(" ", ""), (
                f"Replica {i}: min.in should have ntr=1 for HA restraints"
            )
            assert "!@H=" in min_in, (
                f"Replica {i}: min.in should contain HA mask '!@H='"
            )
            assert "restraint_wt=10.0," in min_in.replace(" ", ""), (
                f"Replica {i}: min.in should have restraint_wt=10.0"
            )

            # NVT heating (eq1.in): restraints carried over from minimization
            eq1_in = (rep_dir / "eq" / "eq1.in").read_text()
            assert "ntr=1," in eq1_in.replace(" ", ""), (
                f"Replica {i}: eq1.in should have ntr=1 for HA restraints"
            )
            assert "!@H=" in eq1_in, (
                f"Replica {i}: eq1.in should contain HA mask '!@H='"
            )

            # Production: no positional restraints (free production run)
            prod_in = (rep_dir / "md" / "prod.in").read_text()
            assert "ntr=0," in prod_in.replace(" ", ""), (
                f"Replica {i}: prod.in should have ntr=0 (no restraints in production)"
            )

            # ---- 5. Replica metadata: state READY, 100 ns -------------------
            replica_data = json.loads((rep_dir / "replica.json").read_text())
            assert replica_data.get("state") == "READY", (
                f"Replica {i}: expected state READY, got {replica_data.get('state')}"
            )
            # nanos is stored inside the nested settings dict
            settings_data = replica_data.get("settings") or {}
            assert settings_data.get("nanos") == 100, (
                f"Replica {i}: expected nanos=100 in settings, got {settings_data.get('nanos')}"
            )
            assert settings_data.get("restraint_mode", "").upper() == "HA", (
                f"Replica {i}: expected restraint_mode=HA, got {settings_data.get('restraint_mode')}"
            )

    def test_solvated_pdb_has_more_atoms_than_dry_system(self, _project_dir):
        """
        Solvated PDB produced by tleap contains more atoms than the dry input.

        tleap embeds all atoms (protein + solvent + counter-ions) into
        ``systems/ETA/pep_ETA.pdb``.  The dry reference is ``pep.pdb`` which
        holds only the 122 peptide atoms.  After solvation the file must
        contain strictly more ATOM / HETATM records.
        """
        solvated_pdb = _project_dir / "systems" / "ETA" / "pep_ETA.pdb"
        assert solvated_pdb.exists(), (
            f"Solvated PDB not found: {solvated_pdb}"
        )

        dry_atoms = self._count_pdb_atoms(_PEP_PDB)
        solvated_atoms = self._count_pdb_atoms(solvated_pdb)

        assert solvated_atoms > dry_atoms, (
            f"Solvated PDB should have more atoms than the dry system "
            f"({solvated_atoms} <= {dry_atoms}).  "
            f"Solvation may have failed or the PDB was not written correctly."
        )
