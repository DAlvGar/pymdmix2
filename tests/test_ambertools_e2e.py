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
"""

from __future__ import annotations

import os
import shutil
import subprocess
from importlib.resources import files
from pathlib import Path

import pytest

# ---------------------------------------------------------------------------
# Bundled test data paths (accessed via importlib.resources)
# ---------------------------------------------------------------------------

_DATA_DIR = Path(str(files("pymdmix").joinpath("data")))
_PEP_DIR = _DATA_DIR / "test" / "pep"
_SOLVENTS_DIR = _DATA_DIR / "solvents"

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
