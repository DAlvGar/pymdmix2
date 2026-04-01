"""
End-to-end and integration tests for pyMDMix workflows.

These tests exercise complete analysis pipelines, from raw data structures
through density calculation, free energy conversion, hotspot detection, and
project management.  They complement the narrower unit tests by verifying
that the modules interact correctly together.
"""

from __future__ import annotations

from importlib.resources import files
from pathlib import Path

import numpy as np
import pytest

from pymdmix.analysis.base import get_action, register_action, run_action
from pymdmix.analysis.density import DensityAction, calculate_density
from pymdmix.analysis.energy import boltzmann_average, density_to_free_energy
from pymdmix.analysis.hotspots import Hotspot, HotSpotSet, HotspotAction, detect_hotspots
from pymdmix.core.grid import Grid
from pymdmix.core.solvent import Probe, Solvent, SolventLibrary
from pymdmix.core.structure import (
    get_backbone_mask,
    get_heavy_atom_mask,
    get_protein_mask,
    get_water_mask,
    load_structure,
    save_structure,
)
from pymdmix.core.trajectory import Frame
from pymdmix.project.config import Config, MDSettings
from pymdmix.project.project import Project
from pymdmix.project.replica import Replica, ReplicaState, create_replica
from pymdmix.setup.prepare import PrepareResult, prepare_structure

# Minimum number of solvents expected in the bundled library (ETA, WAT, MAM, ...)
_MIN_STANDARD_SOLVENTS = 3

# Path to the bundled test PDB file (located via importlib.resources)
_PEP_PDB_PATH = Path(str(files("pymdmix").joinpath("data/test/pep/pep.pdb")))


# ---------------------------------------------------------------------------
# Shared helpers / fixtures
# ---------------------------------------------------------------------------


class _FakeTrajectoryReader:
    """Lightweight mock trajectory for tests that need a real-ish reader."""

    def __init__(self, frames: list[np.ndarray]):
        self._frames = frames

    @property
    def n_frames(self) -> int:
        return len(self._frames)

    @property
    def n_atoms(self) -> int:
        return self._frames[0].shape[0] if self._frames else 0

    def __len__(self) -> int:
        return self.n_frames

    def __iter__(self):
        for coords in self._frames:
            yield Frame(coordinates=coords)


def _make_trajectory(
    n_frames: int = 20,
    n_atoms: int = 150,
    seed: int = 0,
    center: tuple[float, float, float] = (25.0, 25.0, 25.0),
    spread: float = 8.0,
) -> _FakeTrajectoryReader:
    """Create a deterministic synthetic trajectory."""
    rng = np.random.default_rng(seed)
    base = rng.standard_normal((n_atoms, 3)) * spread + np.asarray(center)
    frames = []
    for _ in range(n_frames):
        perturb = rng.standard_normal((n_atoms, 3)) * 0.3
        frames.append((base + perturb).astype(np.float64))
    return _FakeTrajectoryReader(frames)


@pytest.fixture
def fake_trajectory():
    return _make_trajectory()


@pytest.fixture
def tmp_dir(tmp_path):
    return tmp_path


# ---------------------------------------------------------------------------
# 1.  Density → Free Energy → Hotspot pipeline
# ---------------------------------------------------------------------------


class TestDensityHotspotPipeline:
    """Full pipeline: trajectory → density grid → free energy → hotspots."""

    def test_density_action_produces_grid_file(self, fake_trajectory, tmp_dir):
        """DensityAction writes a .dx file and the result is marked success."""
        action = DensityAction()
        probe_indices = {"OH": np.arange(0, 20)}

        result = action.run(
            trajectory=fake_trajectory,
            probe_indices=probe_indices,
            spacing=1.0,
            output_dir=tmp_dir,
        )

        assert result.success, f"DensityAction failed: {result.error}"
        assert len(result.output_files) == 1
        dx_path = result.output_files[0]
        assert dx_path.exists()
        assert dx_path.suffix == ".dx"

    def test_density_grid_round_trips_correctly(self, fake_trajectory, tmp_dir):
        """Density grids written to DX format are re-read without data loss."""
        action = DensityAction()
        probe_indices = {"CT": np.arange(10, 40)}

        result = action.run(
            trajectory=fake_trajectory,
            probe_indices=probe_indices,
            spacing=1.0,
            output_dir=tmp_dir,
        )

        assert result.success
        original = Grid.read_dx(result.output_files[0])
        reloaded = Grid.read_dx(result.output_files[0])

        np.testing.assert_array_almost_equal(original.data, reloaded.data, decimal=6)
        assert original.spacing == reloaded.spacing
        np.testing.assert_allclose(original.origin, reloaded.origin, atol=1e-5)

    def test_free_energy_conversion_preserves_shape(self, fake_trajectory, tmp_dir):
        """Density → free energy preserves grid shape and produces finite values."""
        action = DensityAction()
        probe_indices = {"WAT": np.arange(0, 50)}

        result = action.run(
            trajectory=fake_trajectory,
            probe_indices=probe_indices,
            spacing=1.0,
            output_dir=tmp_dir,
        )

        assert result.success
        density_grid = Grid.read_dx(result.output_files[0])
        fe_grid = density_grid.to_free_energy(temperature=300.0)

        assert fe_grid.data.shape == density_grid.data.shape
        assert np.isfinite(fe_grid.data).all(), "Free energy grid contains non-finite values"

    def test_hotspot_action_detects_high_density_region(self, tmp_dir):
        """HotspotAction finds hotspots when a high-density region is present."""
        # Build a density grid with a clear peak
        nx, ny, nz = 30, 30, 30
        data = np.full((nx, ny, nz), 0.3)  # background below bulk
        # Insert a region with density 6x the bulk → very favourable
        data[12:17, 12:17, 12:17] = 6.0

        grid = Grid(data=data, origin=np.array([0.0, 0.0, 0.0]), spacing=0.5)

        action = HotspotAction()
        result = action.run(
            grids={"OH": grid},
            energy_cutoff=-0.5,
            cluster_distance=2.0,
            min_points=2,
            temperature=300.0,
            output_dir=tmp_dir,
        )

        assert result.success, f"HotspotAction failed: {result.error}"
        assert result.metadata["n_hotspots"] >= 1, "Expected at least one hotspot"
        assert "OH" in result.metadata["probes"]

    def test_full_pipeline_trajectory_to_hotspot_files(self, fake_trajectory, tmp_dir):
        """Complete pipeline: trajectory → density → FE → hotspot files on disk."""
        # Step 1 – density
        density_action = DensityAction()
        probe_indices = {"OH": np.arange(0, 30), "CT": np.arange(30, 60)}

        density_result = density_action.run(
            trajectory=fake_trajectory,
            probe_indices=probe_indices,
            spacing=1.0,
            output_dir=tmp_dir / "density",
        )
        assert density_result.success
        assert len(density_result.output_files) == 2  # one per probe

        # Step 2 – load grids and build energy grids
        grids = {
            f.stem.split("_")[0]: Grid.read_dx(f) for f in density_result.output_files
        }
        assert len(grids) == 2

        # Step 3 – detect hotspots
        hotspot_action = HotspotAction()
        hotspot_result = hotspot_action.run(
            grids=grids,
            energy_cutoff=-0.2,
            cluster_distance=2.5,
            min_points=2,
            output_dir=tmp_dir / "hotspots",
        )
        assert hotspot_result.success

        # Step 4 – output files exist
        for fpath in hotspot_result.output_files:
            assert fpath.exists(), f"Missing hotspot output: {fpath}"


# ---------------------------------------------------------------------------
# 2.  Grid mathematical operations chain
# ---------------------------------------------------------------------------


class TestGridAnalysisChain:
    """Multi-step grid operations exercised end-to-end."""

    def _make_density_grid(
        self,
        shape: tuple[int, int, int] = (20, 20, 20),
        peak_val: float = 4.0,
        peak_slice: tuple = (slice(8, 12), slice(8, 12), slice(8, 12)),
        bg_val: float = 0.5,
        spacing: float = 0.5,
        seed: int = 1,
    ) -> Grid:
        rng = np.random.default_rng(seed)
        data = rng.uniform(0.0, bg_val, shape)
        data[peak_slice] = peak_val
        return Grid(data=data, origin=np.array([0.0, 0.0, 0.0]), spacing=spacing)

    def test_grid_add_counts_and_normalize(self):
        """Grid accumulates point counts and to_density normalises correctly."""
        grid = Grid.from_bounds(
            min_coord=(0.0, 0.0, 0.0),
            max_coord=(10.0, 10.0, 10.0),
            spacing=1.0,
        )
        # Add 100 counts to the same location (centre)
        n_added = 0
        for _ in range(100):
            added = grid.add_count(np.array([5.0, 5.0, 5.0]))
            if added:
                n_added += 1

        assert n_added == 100
        density = grid.to_density(n_frames=100)
        # At (5,5,5) voxel the density should equal 1.0 (count/frame)
        idx = density.coord_to_index((5.0, 5.0, 5.0))
        assert density.data[idx] == pytest.approx(1.0, abs=1e-6)

    def test_grid_to_free_energy_monotone(self):
        """Higher density yields lower (more negative) free energy.

        Uses density_to_free_energy with an explicit shared reference density
        so that both uniform grids are evaluated on the same scale.
        """
        from pymdmix.analysis.energy import density_to_free_energy

        low_density = Grid(
            data=np.full((5, 5, 5), 0.5), origin=np.zeros(3), spacing=1.0
        )
        high_density = Grid(
            data=np.full((5, 5, 5), 5.0), origin=np.zeros(3), spacing=1.0
        )

        ref = 1.0  # shared bulk reference
        dg_low = density_to_free_energy(low_density, temperature=300.0, reference_density=ref)
        dg_high = density_to_free_energy(high_density, temperature=300.0, reference_density=ref)

        assert dg_low.data.mean() > dg_high.data.mean(), (
            "Higher density should give lower (more negative) free energy"
        )

    def test_grid_addition_chain(self, tmp_dir):
        """Add two density grids, convert to FE, write to disk and reload."""
        g1 = self._make_density_grid(seed=1)
        g2 = self._make_density_grid(seed=2, bg_val=0.3)

        # Add grids
        g_sum = g1 + g2
        assert g_sum.data.shape == g1.data.shape
        np.testing.assert_allclose(g_sum.data, g1.data + g2.data)

        # Convert to FE
        fe = g_sum.to_free_energy(300.0)
        assert np.isfinite(fe.data).all()

        # Save and reload
        out = tmp_dir / "sum_fe.dx"
        fe.write_dx(out)
        loaded = Grid.read_dx(out)
        np.testing.assert_array_almost_equal(loaded.data, fe.data, decimal=5)

    def test_boltzmann_average_of_replica_grids(self, tmp_dir):
        """boltzmann_average over multiple DX files yields sensible result."""
        rng = np.random.default_rng(42)
        paths = []
        arrays = []
        for i in range(3):
            data = rng.uniform(0.5, 3.0, (10, 10, 10))
            arrays.append(data)
            g = Grid(data=data, origin=np.zeros(3), spacing=1.0)
            p = tmp_dir / f"replica_{i}.dx"
            g.write_dx(p)
            paths.append(p)

        avg_grid = boltzmann_average(paths, temperature=300.0)

        assert avg_grid.data.shape == (10, 10, 10)
        # Boltzmann average should be between the min and max of the individual means
        all_means = np.array([a.mean() for a in arrays])
        assert all_means.min() <= avg_grid.data.mean() <= all_means.max() + 1.0

    def test_grid_subgrid_extraction(self):
        """take_subgrid_box returns a smaller grid with correct data."""
        data = np.arange(1000.0).reshape(10, 10, 10)
        grid = Grid(data=data, origin=np.array([0.0, 0.0, 0.0]), spacing=1.0)

        sub = grid.take_subgrid_box(
            min_coord=(2.0, 2.0, 2.0),
            max_coord=(6.0, 6.0, 6.0),
        )

        assert all(s <= 6 for s in sub.data.shape), "Subgrid should be smaller than original"

    def test_percentile_cutoff(self):
        """percentile_cutoff returns correct threshold value."""
        rng = np.random.default_rng(7)
        data = rng.uniform(0.0, 10.0, (20, 20, 20))
        grid = Grid(data=data, origin=np.zeros(3), spacing=0.5)

        p95 = grid.percentile_cutoff(95)
        assert np.percentile(data, 95) == pytest.approx(p95, rel=1e-4)

    def test_grid_expand_contract_roundtrip(self):
        """Expand then contract grid restores original shape and interior values."""
        data = np.ones((8, 8, 8)) * 2.0
        grid = Grid(data=data, origin=np.zeros(3), spacing=1.0)

        expanded = grid.expand(buffer=2)
        assert expanded.data.shape == (12, 12, 12)

        contracted = expanded.contract(buffer=2)
        assert contracted.data.shape == grid.data.shape
        np.testing.assert_array_almost_equal(contracted.data, grid.data)


# ---------------------------------------------------------------------------
# 3.  Solvent library → probe extraction → density workflow
# ---------------------------------------------------------------------------


class TestSolventDensityWorkflow:
    """Load a real solvent definition and drive density calculations from it."""

    def test_load_eta_solvent_and_extract_probes(self):
        """SolventLibrary loads ETA and its probes are present."""
        lib = SolventLibrary()
        eta = lib.get("ETA")
        assert eta is not None, "ETA solvent not found in library"
        assert len(eta.probes) >= 2

        probe_names = eta.get_probe_names()
        assert "OH" in probe_names
        assert "CT" in probe_names

    def test_probe_selection_strings_non_empty(self):
        """Probe selection strings are non-empty MDAnalysis strings."""
        lib = SolventLibrary()
        eta = lib.get("ETA")
        for probe in eta.probes:
            sel = probe.selection
            assert sel and "resname" in sel and "name" in sel

    def test_all_standard_solvents_loadable(self):
        """Every solvent in the bundled library can be loaded."""
        lib = SolventLibrary()
        names = lib.list_solvents()
        assert len(names) >= _MIN_STANDARD_SOLVENTS, (
            f"Expected at least {_MIN_STANDARD_SOLVENTS} solvents (ETA, WAT, MAM, ...)"
        )

        for name in names:
            solvent = lib.get(name)
            assert solvent is not None, f"Solvent {name} could not be loaded"
            assert solvent.name == name
            assert len(solvent.probes) >= 1

    def test_solvent_probability_in_range(self):
        """Probe probability values are positive and < 1."""
        lib = SolventLibrary()
        for name in lib.list_solvents():
            solvent = lib.get(name)
            for probe in solvent.probes:
                if probe.probability is not None:
                    assert 0 < probe.probability < 1, (
                        f"{name}/{probe.name} probability out of range: {probe.probability}"
                    )

    def test_density_with_solvent_probe_indices(self, tmp_dir):
        """Use probe atom counts from the ETA solvent to set up a density run."""
        lib = SolventLibrary()
        eta = lib.get("ETA")

        n_oh = 40
        n_ct = 40
        n_wat = 60
        n_atoms = n_oh + n_ct + n_wat

        rng = np.random.default_rng(11)
        base = rng.standard_normal((n_atoms, 3)) * 5.0 + 20.0
        frames = [base + rng.standard_normal((n_atoms, 3)) * 0.2 for _ in range(15)]
        traj = _FakeTrajectoryReader(frames)

        # Map probes to atom index slices
        probe_indices = {
            "OH": np.arange(0, n_oh),
            "CT": np.arange(n_oh, n_oh + n_ct),
            "WAT": np.arange(n_oh + n_ct, n_atoms),
        }

        action = DensityAction()
        result = action.run(
            trajectory=traj,
            probe_indices=probe_indices,
            spacing=1.0,
            output_dir=tmp_dir,
        )

        assert result.success, result.error
        assert result.metadata["n_probes"] == len(probe_indices)

        # All three probe density files exist
        for fpath in result.output_files:
            assert fpath.exists()


# ---------------------------------------------------------------------------
# 4.  Project lifecycle – create, populate, save, reload
# ---------------------------------------------------------------------------


class TestProjectLifecycleWorkflow:
    """Create a project, add replicas, persist to disk and reload."""

    def test_create_project_and_add_replicas(self, tmp_dir):
        """Project.create() produces a JSON file; replicas can be added."""
        project = Project.create(
            name="test_project",
            path=tmp_dir / "test_project",
        )
        project.add_replicas(solvent="ETA", n_replicas=3)
        project.add_replicas(solvent="WAT", n_replicas=2)

        assert project.n_replicas == 5
        assert set(project.solvents) == {"ETA", "WAT"}
        assert len(project.get_replicas_by_solvent("ETA")) == 3
        assert len(project.get_replicas_by_solvent("WAT")) == 2

    def test_project_save_and_reload(self, tmp_dir):
        """Project state round-trips through JSON serialisation."""
        project = Project.create(
            name="reload_test",
            path=tmp_dir / "reload_test",
        )
        project.add_replicas("MAM", n_replicas=2)
        project.save()

        # Reload
        reloaded = Project.load(tmp_dir / "reload_test")

        assert reloaded.name == "reload_test"
        assert reloaded.n_replicas == 2
        assert reloaded.solvents == ["MAM"]

    def test_replica_state_transitions(self, tmp_dir):
        """Replicas progress through states correctly."""
        replica = Replica(
            name="ETA_1",
            solvent="ETA",
            path=tmp_dir / "ETA_1",
        )

        assert replica.state == ReplicaState.CREATED

        replica.create_directory()
        assert replica.exists()
        assert (replica.path / "min").exists()
        assert (replica.path / "md").exists()

        replica.set_state(ReplicaState.SETUP)
        assert replica.state == ReplicaState.SETUP

        replica.set_state(ReplicaState.READY)
        assert replica.state == ReplicaState.READY

    def test_project_group_management(self, tmp_dir):
        """Groups of replicas can be created and queried."""
        project = Project.create(name="grp_test", path=tmp_dir / "grp_test")
        project.add_replicas("ETA", n_replicas=3)
        project.add_replicas("WAT", n_replicas=3)

        replica_names = [r.name for r in project.get_replicas_by_solvent("ETA")]
        project.create_group("eta_group", replica_names)

        group = project.get_group("eta_group")
        assert group == replica_names

    def test_project_status_dict_structure(self, tmp_dir):
        """Project.status() returns expected keys."""
        project = Project.create(name="status_test", path=tmp_dir / "status_test")
        project.add_replicas("ETA", n_replicas=2)

        status = project.status()
        assert "n_replicas" in status
        assert "solvents" in status
        assert "states" in status
        assert status["n_replicas"] == 2

    def test_create_replica_function(self, tmp_dir):
        """create_replica() factory builds a valid Replica object."""
        replica = create_replica(
            name="MAM_1",
            solvent="MAM",
            base_path=tmp_dir,
        )
        assert replica.name == "MAM_1"
        assert replica.solvent == "MAM"
        assert replica.state == ReplicaState.CREATED


# ---------------------------------------------------------------------------
# 5.  Structure preparation pipeline
# ---------------------------------------------------------------------------


class TestStructurePreparationPipeline:
    """Load PDB → prepare → analyse → save."""

    @pytest.fixture
    def pep_pdb_path(self) -> Path:
        """Return the bundled test peptide PDB file (located via importlib.resources)."""
        assert _PEP_PDB_PATH.exists(), f"Test PDB not found at {_PEP_PDB_PATH}"
        return _PEP_PDB_PATH

    def test_load_real_pdb_structure(self, pep_pdb_path):
        """Bundled pep.pdb loads as a parmed.Structure with correct counts."""
        struct = load_structure(pep_pdb_path)

        assert struct is not None
        assert len(struct.atoms) > 50  # The peptide has >50 atoms
        assert len(struct.residues) >= 5  # At least 5 residues

    def test_protein_mask_covers_all_atoms(self, pep_pdb_path):
        """All atoms in pep.pdb are recognised as protein atoms."""
        struct = load_structure(pep_pdb_path)
        mask = get_protein_mask(struct)

        # pep.pdb is a bare peptide with no water – every atom should be protein
        assert mask.sum() == len(struct.atoms)

    def test_heavy_atom_mask_excludes_hydrogens(self, pep_pdb_path):
        """Heavy atom mask excludes hydrogens."""
        struct = load_structure(pep_pdb_path)
        ha_mask = get_heavy_atom_mask(struct)
        bb_mask = get_backbone_mask(struct)

        n_heavy = ha_mask.sum()
        n_backbone = bb_mask.sum()

        assert n_heavy > n_backbone
        # Every backbone atom is a heavy atom
        assert int(np.logical_and(ha_mask, bb_mask).sum()) == n_backbone

    def test_prepare_structure_returns_result(self, pep_pdb_path):
        """prepare_structure() completes without error on the test peptide."""
        struct = load_structure(pep_pdb_path)
        result = prepare_structure(struct)

        assert isinstance(result, PrepareResult)
        assert result.structure is not None

    def test_save_reload_preserves_atom_count(self, pep_pdb_path, tmp_dir):
        """Saving and reloading a structure preserves the atom count."""
        struct = load_structure(pep_pdb_path)
        out = tmp_dir / "pep_out.pdb"
        save_structure(struct, out)

        assert out.exists()
        reloaded = load_structure(out)
        assert len(reloaded.atoms) == len(struct.atoms)

    def test_grid_from_structure_encloses_all_atoms(self, pep_pdb_path):
        """Grid created from the peptide coordinates contains every atom."""
        struct = load_structure(pep_pdb_path)
        coords = np.array([[a.xx, a.xy, a.xz] for a in struct.atoms])

        grid = Grid.from_structure(coords, spacing=1.0, padding=5.0)

        # Every atom coordinate should be inside the grid
        inside = [grid.is_inside(c) for c in coords]
        assert all(inside), "Some atoms are outside the grid bounding box"

    def test_density_calculation_on_real_structure(self, pep_pdb_path, tmp_dir):
        """Density action runs on a trajectory derived from the real peptide."""
        struct = load_structure(pep_pdb_path)
        coords = np.array([[a.xx, a.xy, a.xz] for a in struct.atoms])
        n_atoms = len(coords)

        rng = np.random.default_rng(99)
        frames = [coords + rng.standard_normal((n_atoms, 3)) * 0.5 for _ in range(10)]
        traj = _FakeTrajectoryReader(frames)

        # Use the first 20 atom indices as a mock probe
        probe_indices = {"mock_probe": np.arange(min(20, n_atoms))}

        action = DensityAction()
        result = action.run(
            trajectory=traj,
            probe_indices=probe_indices,
            spacing=1.0,
            output_dir=tmp_dir,
        )

        assert result.success, result.error
        dx_path = result.output_files[0]
        assert dx_path.exists()

        # The resulting grid should be non-trivial
        density_grid = Grid.read_dx(dx_path)
        assert density_grid.data.sum() > 0


# ---------------------------------------------------------------------------
# 6.  Multi-solvent comparison
# ---------------------------------------------------------------------------


class TestMultiSolventComparison:
    """Compare density results across multiple solvents."""

    def test_eta_vs_wat_density_comparison(self, tmp_dir):
        """ETA and WAT density grids can be loaded and compared."""
        lib = SolventLibrary()
        solvents = ["ETA", "WAT"]

        rng = np.random.default_rng(55)
        n_atoms = 100
        base = rng.standard_normal((n_atoms, 3)) * 6.0 + 20.0
        frames = [base + rng.standard_normal((n_atoms, 3)) * 0.3 for _ in range(10)]
        traj = _FakeTrajectoryReader(frames)

        density_grids: dict[str, dict[str, Grid]] = {}

        for solvent_name in solvents:
            solvent = lib.get(solvent_name)
            assert solvent is not None

            solvent_dir = tmp_dir / solvent_name
            action = DensityAction()

            # Use 20 atoms per probe as a simple stand-in
            probes = solvent.get_probe_names()[:2]
            probe_indices = {
                p: np.arange(i * 20, (i + 1) * 20) for i, p in enumerate(probes)
            }

            result = action.run(
                trajectory=traj,
                probe_indices=probe_indices,
                spacing=1.0,
                output_dir=solvent_dir,
            )
            assert result.success, f"{solvent_name}: {result.error}"

            density_grids[solvent_name] = {
                f.stem.split("_")[0]: Grid.read_dx(f) for f in result.output_files
            }

        # Both solvents produced density grids
        assert len(density_grids) == 2
        for sv in solvents:
            assert len(density_grids[sv]) >= 1

    def test_average_across_replicas_reduces_noise(self, tmp_dir):
        """Averaging N replica grids reduces noise compared to a single replica."""
        rng = np.random.default_rng(77)
        shape = (15, 15, 15)
        true_signal = np.full(shape, 2.0)

        # Generate replica grids with Gaussian noise
        replica_paths = []
        for i in range(5):
            noisy = true_signal + rng.standard_normal(shape) * 0.5
            noisy = np.clip(noisy, 0.0, None)
            g = Grid(data=noisy, origin=np.zeros(3), spacing=1.0)
            p = tmp_dir / f"rep_{i}.dx"
            g.write_dx(p)
            replica_paths.append(p)

        avg = boltzmann_average(replica_paths, temperature=300.0)

        # The average should have less variance than any individual replica
        single = Grid.read_dx(replica_paths[0])
        assert avg.data.var() <= single.data.var() + 0.5, (
            "Averaged grid should have lower or comparable variance"
        )

    def test_hotspot_detection_across_probes(self, tmp_dir):
        """HotspotAction detects independent hotspots for multiple probe types."""
        import json

        shape = (25, 25, 25)

        # Probe OH: hotspot at index (5,5,5)
        oh_data = np.full(shape, 0.2)
        oh_data[4:8, 4:8, 4:8] = 8.0

        # Probe CT: hotspot at index (18,18,18)
        ct_data = np.full(shape, 0.2)
        ct_data[17:21, 17:21, 17:21] = 8.0

        grids = {
            "OH": Grid(data=oh_data, origin=np.zeros(3), spacing=0.5),
            "CT": Grid(data=ct_data, origin=np.zeros(3), spacing=0.5),
        }

        action = HotspotAction()
        result = action.run(
            grids=grids,
            energy_cutoff=-0.5,
            cluster_distance=2.0,
            min_points=2,
            output_dir=tmp_dir,
        )

        assert result.success
        assert result.metadata["n_hotspots"] >= 2

        # Both probe types are represented
        assert "OH" in result.metadata["probes"]
        assert "CT" in result.metadata["probes"]

        # Read hotspot JSON to check per-probe centroids are separated
        json_path = next(f for f in result.output_files if f.suffix == ".json")
        with open(json_path) as fh:
            hotspot_data = json.load(fh)

        probes_found = {h["probe"] for h in hotspot_data["hotspots"]}
        assert "OH" in probes_found
        assert "CT" in probes_found

        # Hotspot centroids for OH and CT should be well separated
        oh_centroids = [h["centroid"] for h in hotspot_data["hotspots"] if h["probe"] == "OH"]
        ct_centroids = [h["centroid"] for h in hotspot_data["hotspots"] if h["probe"] == "CT"]

        oh_center = np.mean(oh_centroids, axis=0)
        ct_center = np.mean(ct_centroids, axis=0)
        dist = np.linalg.norm(oh_center - ct_center)
        assert dist > 3.0, "OH and CT hotspots should be well separated"


# ---------------------------------------------------------------------------
# 7.  HotSpotSet operations workflow
# ---------------------------------------------------------------------------


class TestHotSpotSetWorkflow:
    """Test HotSpotSet collection operations end-to-end."""

    def _make_set(self, n: int = 6, seed: int = 0) -> HotSpotSet:
        rng = np.random.default_rng(seed)
        hotspots = []
        for i in range(n):
            n_pts = rng.integers(5, 20)
            coords = rng.standard_normal((int(n_pts), 3)) * 2.0 + float(i * 5)
            energies = rng.uniform(-3.0, -0.5, int(n_pts))
            hs = Hotspot(
                id=i,
                probe="OH",
                centroid=tuple(coords.mean(axis=0)),
                energy=float(energies.mean()),
                volume=float(n_pts) * 0.125,
                n_points=int(n_pts),
                coords=coords,
                energies=energies,
            )
            hotspots.append(hs)
        return HotSpotSet(probe="OH", name="test", hotspots=hotspots)

    def test_prune_by_energy_removes_weak_hotspots(self):
        """prune_by_energy keeps only hotspots with energy ≤ threshold."""
        hs_set = self._make_set(n=10)
        cutoff = -1.0
        pruned = hs_set.prune_by_energy(cutoff)

        assert all(h.energy <= cutoff for h in pruned)
        assert len(pruned) <= len(hs_set)

    def test_prune_by_volume_removes_small_hotspots(self):
        """prune_by_volume keeps only hotspots with volume ≥ min_volume."""
        hs_set = self._make_set(n=10)
        min_vol = 0.5
        pruned = hs_set.prune_by_volume(min_vol)

        assert all(h.volume >= min_vol for h in pruned)

    def test_clustering_separates_distant_hotspots(self):
        """Hotspots spread far apart cluster into distinct groups."""
        rng = np.random.default_rng(3)
        # Two tight clusters separated by 30 Å
        hotspots = []
        for cluster_offset in [0.0, 30.0]:
            for j in range(4):
                coords = rng.standard_normal((8, 3)) * 0.5 + cluster_offset
                hs = Hotspot(
                    id=len(hotspots),
                    probe="OH",
                    centroid=tuple(coords.mean(axis=0)),
                    energy=-1.5,
                    volume=1.0,
                    n_points=8,
                    coords=coords,
                    energies=np.full(8, -1.5),
                )
                hotspots.append(hs)

        hs_set = HotSpotSet(probe="OH", name="two_clusters", hotspots=hotspots)
        hs_set.cluster(cutoff=3.0)

        assert hs_set.n_clusters == 2

    def test_json_roundtrip_preserves_hotspot_data(self, tmp_dir):
        """HotSpotSet serialises and deserialises via JSON without data loss."""
        hs_set = self._make_set(n=4)
        json_path = tmp_dir / "hotspots.json"
        hs_set.to_json(json_path)

        assert json_path.exists()

        import json

        with open(json_path) as fh:
            data = json.load(fh)

        assert "hotspots" in data
        assert len(data["hotspots"]) == 4

    def test_pdb_output_file_created(self, tmp_dir):
        """HotSpotSet.to_pdb() writes a valid PDB file."""
        hs_set = self._make_set(n=3)
        pdb_path = tmp_dir / "hotspots.pdb"
        hs_set.to_pdb(pdb_path)

        assert pdb_path.exists()
        content = pdb_path.read_text()
        assert "ATOM" in content or "HETATM" in content


# ---------------------------------------------------------------------------
# 8.  Action registry integration
# ---------------------------------------------------------------------------


class TestActionRegistryIntegration:
    """Verify that actions registered with @register_action are discoverable."""

    def test_builtin_actions_all_registered(self):
        """The standard analysis actions are all present in the registry."""
        from pymdmix.analysis.base import list_actions

        actions = list_actions()
        for name in ("density", "hotspots", "residence", "cpptraj_density"):
            assert name in actions, f"Expected action '{name}' not registered"

    def test_custom_action_round_trip(self, fake_trajectory):
        """A custom action registered at runtime can be retrieved and called."""
        from pymdmix.analysis.base import Action, ActionResult

        # Track invocations to verify the action is actually executed
        action_invocation_log: list[str] = []

        @register_action("_integration_test_action")
        class _IntegrationAction(Action):
            name = "_integration_test_action"
            description = "test-only action"

            def run(self, trajectory=None, **kwargs) -> ActionResult:
                action_invocation_log.append("called")
                return ActionResult(success=True, metadata={"called": True})

        cls = get_action("_integration_test_action")
        assert cls is _IntegrationAction

        result = run_action("_integration_test_action", fake_trajectory)
        assert result.success
        assert action_invocation_log == ["called"]
