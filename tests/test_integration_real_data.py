"""
Integration tests using bundled real test data.

These tests load the actual data files shipped in ``pymdmix/data/test/pep/``
and the solvent JSON definitions to verify that the I/O layer, parsers, and
analysis modules work correctly end-to-end with representative inputs.

Data used
---------
* ``pymdmix/data/test/pep/pep.pdb``       – peptide structure (8 residues)
* ``pymdmix/data/test/pep/pep.prmtop``    – Amber topology
* ``pymdmix/data/test/pep/pep.off``       – LEaP object file
* ``pymdmix/data/test/pep/pep_amber_mdmix.cfg`` – legacy Amber MDMix config
* ``pymdmix/data/test/pep/multisettings.cfg``  – multi-section settings file
* ``pymdmix/data/solvents/*.json``         – solvent definitions
"""

from __future__ import annotations

import json
import os
from importlib.resources import files
from pathlib import Path

import numpy as np
import pytest

# ---------------------------------------------------------------------------
# Data path helpers
# ---------------------------------------------------------------------------

# Use importlib.resources so the paths work whether pymdmix is installed as a
# regular package or as an editable install.
_DATA_DIR = Path(str(files("pymdmix").joinpath("data")))
_PEP_DIR = _DATA_DIR / "test" / "pep"
_SOLVENTS_DIR = _DATA_DIR / "solvents"

# Minimum number of solvent JSON files expected in the bundled library
# (ETA, WAT, MAM, ISO, ANT, MOH, ION, PYR, ISO5 are all present)
_MIN_BUNDLED_SOLVENTS = 4


def _pep(name: str) -> Path:
    """Return full path of a pep test data file."""
    return _PEP_DIR / name


def _solvent(name: str) -> Path:
    """Return full path of a solvent data file."""
    return _SOLVENTS_DIR / name


# ---------------------------------------------------------------------------
# 1.  Real PDB structure analysis
# ---------------------------------------------------------------------------


class TestRealPDBStructure:
    """Tests using the bundled pep.pdb test peptide."""

    def test_pdb_exists(self):
        """The bundled pep.pdb test file is present."""
        assert _pep("pep.pdb").exists(), f"Missing test PDB: {_pep('pep.pdb')}"

    def test_load_peptide_pdb(self):
        """pep.pdb loads as a parmed.Structure with the expected residue count."""
        from pymdmix.core.structure import load_structure

        struct = load_structure(_pep("pep.pdb"))
        # pep.pdb contains 8 residues (GLN-PHE-GLY-TRP-SER-PHE-TYR-ALA)
        assert len(struct.residues) == 8
        assert len(struct.atoms) > 100  # 122 heavy + H atoms

    def test_peptide_has_no_water(self):
        """The bundled peptide PDB contains no water residues."""
        from pymdmix.core.structure import get_water_mask, load_structure

        struct = load_structure(_pep("pep.pdb"))
        water_mask = get_water_mask(struct)
        assert water_mask.sum() == 0, "pep.pdb should contain no water molecules"

    def test_all_atoms_are_protein(self):
        """Every atom in pep.pdb is identified as a protein atom."""
        from pymdmix.core.structure import get_protein_mask, load_structure

        struct = load_structure(_pep("pep.pdb"))
        mask = get_protein_mask(struct)
        assert mask.sum() == len(struct.atoms)

    def test_heavy_atom_mask_excludes_hydrogens(self):
        """Heavy-atom mask gives fewer atoms than the full atom count."""
        from pymdmix.core.structure import get_heavy_atom_mask, load_structure

        struct = load_structure(_pep("pep.pdb"))
        ha_mask = get_heavy_atom_mask(struct)

        assert ha_mask.sum() > 0
        assert ha_mask.sum() < len(struct.atoms), "Heavy atom mask should exclude H atoms"

    def test_backbone_mask_cardinality(self):
        """Backbone (N, CA, C, O) has 4 atoms per residue in the peptide."""
        from pymdmix.core.structure import get_backbone_mask, load_structure

        struct = load_structure(_pep("pep.pdb"))
        bb_mask = get_backbone_mask(struct)
        n_residues = len(struct.residues)

        # 4 backbone atoms per residue (last residue has OXT instead, still 4)
        # Allow some tolerance for protonated termini
        assert n_residues * 3 <= bb_mask.sum() <= n_residues * 5

    def test_structure_coordinates_are_finite(self):
        """All atomic coordinates in pep.pdb are finite."""
        from pymdmix.core.structure import load_structure

        struct = load_structure(_pep("pep.pdb"))
        coords = np.array([[a.xx, a.xy, a.xz] for a in struct.atoms])
        assert np.isfinite(coords).all()

    def test_prepare_peptide_structure(self):
        """prepare_structure() completes on pep.pdb and returns PrepareResult."""
        from pymdmix.core.structure import load_structure
        from pymdmix.setup.prepare import PrepareResult, prepare_structure

        struct = load_structure(_pep("pep.pdb"))
        result = prepare_structure(struct)

        assert isinstance(result, PrepareResult)
        assert result.structure is not None
        assert len(result.structure.atoms) > 0

    def test_grid_creation_from_peptide_coordinates(self):
        """Grid.from_structure() creates a grid that contains all peptide atoms."""
        from pymdmix.core.grid import Grid
        from pymdmix.core.structure import load_structure

        struct = load_structure(_pep("pep.pdb"))
        coords = np.array([[a.xx, a.xy, a.xz] for a in struct.atoms])

        grid = Grid.from_structure(coords, spacing=0.5, padding=5.0)

        # All atoms should be inside the grid
        outside = [not grid.is_inside(c) for c in coords]
        assert sum(outside) == 0, f"{sum(outside)} atoms lie outside the grid"

    def test_density_simulation_on_peptide(self, tmp_path):
        """Synthetic density run using peptide atom positions as mock probes."""
        from pymdmix.analysis.density import DensityAction
        from pymdmix.core.grid import Grid
        from pymdmix.core.structure import load_structure
        from pymdmix.core.trajectory import Frame

        struct = load_structure(_pep("pep.pdb"))
        coords = np.array([[a.xx, a.xy, a.xz] for a in struct.atoms])
        n_atoms = len(coords)

        # Build a 10-frame synthetic trajectory with small thermal noise
        rng = np.random.default_rng(42)
        frames = []
        for _ in range(10):
            noise = rng.standard_normal((n_atoms, 3)) * 0.3
            frames.append(coords + noise)

        class _TrajReader:
            def __init__(self, frames):
                self._frames = frames

            @property
            def n_frames(self):
                return len(self._frames)

            @property
            def n_atoms(self):
                return self._frames[0].shape[0]

            def __len__(self):
                return self.n_frames

            def __iter__(self):
                for f in self._frames:
                    yield Frame(coordinates=f)

        traj = _TrajReader(frames)

        # Use backbone-like atoms (first 32) as mock probes
        action = DensityAction()
        result = action.run(
            trajectory=traj,
            probe_indices={"backbone": np.arange(min(32, n_atoms))},
            spacing=0.5,
            output_dir=tmp_path,
        )

        assert result.success, result.error
        dx_path = result.output_files[0]
        density = Grid.read_dx(dx_path)

        # The density around the peptide structure should be non-zero
        assert density.data.sum() > 0


# ---------------------------------------------------------------------------
# 2.  Amber prmtop / topology loading
# ---------------------------------------------------------------------------


class TestAmberTopology:
    """Tests that load the bundled Amber topology file."""

    def test_prmtop_file_exists(self):
        """pep.prmtop is present in the test data directory."""
        assert _pep("pep.prmtop").exists()

    def test_prmtop_has_content(self):
        """pep.prmtop is non-empty and contains expected Amber header."""
        content = _pep("pep.prmtop").read_text()
        assert len(content) > 0
        assert "%FLAG" in content or "TITLE" in content or "ATOM_NAME" in content

    def test_load_structure_from_prmtop(self):
        """parmed can load pep.prmtop without errors."""
        import parmed

        struct = parmed.load_file(str(_pep("pep.prmtop")))
        assert struct is not None
        assert len(struct.atoms) > 100

    def test_load_structure_with_coordinates(self):
        """parmed can load prmtop + prmcrd coordinate file."""
        import parmed

        # prmcrd is the coordinate file
        struct = parmed.load_file(
            str(_pep("pep.prmtop")),
            xyz=str(_pep("pep.prmcrd")),
        )
        assert struct is not None
        coords = struct.coordinates
        assert coords is not None
        assert np.isfinite(coords).all()
        assert coords.shape[1] == 3


# ---------------------------------------------------------------------------
# 3.  OFF (object file format) loading
# ---------------------------------------------------------------------------


class TestOFFFile:
    """Tests using pep.off and the solvent OFF boxes."""

    def test_pep_off_exists(self):
        """pep.off is present in the test data directory."""
        assert _pep("pep.off").exists()

    def test_load_pep_off(self):
        """OFFManager parses pep.off without raising."""
        from pymdmix.io.off_manager import OFFManager

        manager = OFFManager.from_file(_pep("pep.off"))
        assert manager is not None

    def test_pep_off_contains_residue(self):
        """pep.off can be parsed and the 'pep' unit has coordinates."""
        from pymdmix.io.off_manager import OFFManager

        manager = OFFManager.from_file(_pep("pep.off"))
        # get_coords returns the raw coordinates for the unit
        coords = manager.get_coords("pep")
        assert coords is not None
        assert coords.shape[0] > 0, "Expected at least one atom in pep unit"

    def test_pep_off_coordinates(self):
        """pep.off coordinates can be extracted and are finite."""
        from pymdmix.io.off_manager import OFFManager

        manager = OFFManager.from_file(_pep("pep.off"))
        coords = manager.get_coords("pep")
        assert coords is not None
        assert coords.shape[1] == 3
        assert np.isfinite(coords).all()

    def test_solvent_off_files_loadable(self):
        """All bundled solvent OFF files parse without error."""
        from pymdmix.io.off_manager import OFFManager

        off_files = list(_SOLVENTS_DIR.glob("*.off"))
        assert len(off_files) >= 3, "Expected at least 3 bundled solvent OFF files (ETAWAT20, MAMWAT20, etc.)"

        for off_path in off_files:
            manager = OFFManager.from_file(off_path)
            assert manager is not None, f"Failed to load {off_path.name}"


# ---------------------------------------------------------------------------
# 4.  Legacy configuration file parsing
# ---------------------------------------------------------------------------


class TestLegacyConfigParsing:
    """Tests for the legacy MDMix config file parser."""

    def test_amber_config_file_exists(self):
        """pep_amber_mdmix.cfg is present."""
        assert _pep("pep_amber_mdmix.cfg").exists()

    def test_parse_amber_config(self):
        """SystemConfigFileParser reads pep_amber_mdmix.cfg."""
        from pymdmix.io.parsers import SystemConfigFileParser

        parser = SystemConfigFileParser()
        params = parser.parse(_pep("pep_amber_mdmix.cfg"))

        assert params is not None
        # Should have a name and off file
        assert "name" in params or "off" in params or "NAME" in str(params).upper()

    def test_mdsettings_parser_reads_amber_config(self):
        """MDSettingsConfigFileParser extracts settings from pep_amber_mdmix.cfg."""
        from pymdmix.io.parsers import MDSettingsConfigFileParser

        parser = MDSettingsConfigFileParser()
        settings_list = parser.parse(_pep("pep_amber_mdmix.cfg"))

        assert len(settings_list) >= 1

        first = settings_list[0]
        # The parser returns dicts with solvent, temperature, etc.
        assert isinstance(first, dict), f"Expected dict, got {type(first)}"
        assert "solvent" in first or "temp" in first, (
            f"Expected temperature/solvent info in settings, got keys: {list(first.keys())}"
        )

    def test_multisettings_produces_multiple_blocks(self):
        """multisettings.cfg yields more than one settings block."""
        from pymdmix.io.parsers import MDSettingsConfigFileParser

        parser = MDSettingsConfigFileParser()
        settings_list = parser.parse(_pep("multisettings.cfg"))

        # The file has [MDSETTINGS1] and [MDSETTINGS2] sections
        assert len(settings_list) >= 2


# ---------------------------------------------------------------------------
# 5.  Solvent library – JSON round-trip and probe consistency
# ---------------------------------------------------------------------------


class TestSolventJSONData:
    """Tests for the bundled solvent JSON definition files."""

    def test_eta_json_exists(self):
        """eta.json is present in the solvents data directory."""
        assert _solvent("eta.json").exists()

    def test_eta_json_valid_structure(self):
        """eta.json contains all required fields."""
        with open(_solvent("eta.json")) as fh:
            data = json.load(fh)

        required = {"name", "probes", "residues"}
        assert required.issubset(data.keys()), f"Missing keys: {required - set(data)}"
        assert len(data["probes"]) >= 2
        assert data["name"] == "ETA"

    def test_load_eta_from_json(self):
        """Solvent.from_json() loads eta.json and exposes probes correctly."""
        from pymdmix.core.solvent import Solvent

        solvent = Solvent.from_json(_solvent("eta.json"))
        assert solvent.name == "ETA"
        assert len(solvent.probes) >= 2

        probe_names = {p.name for p in solvent.probes}
        assert "OH" in probe_names
        assert "CT" in probe_names

    def test_all_solvent_json_parseable(self):
        """Every *.json file in solvents/ parses without error."""
        from pymdmix.core.solvent import Solvent

        json_files = list(_SOLVENTS_DIR.glob("*.json"))
        assert len(json_files) >= _MIN_BUNDLED_SOLVENTS, (
            f"Expected at least {_MIN_BUNDLED_SOLVENTS} solvent JSON files"
        )

        for json_path in json_files:
            solvent = Solvent.from_json(json_path)
            assert solvent.name, f"{json_path.name}: solvent has no name"
            assert len(solvent.probes) >= 1, f"{json_path.name}: no probes defined"

    def test_solvent_json_roundtrip(self, tmp_path):
        """Solvent.to_json() → Solvent.from_json() preserves all probe data."""
        from pymdmix.core.solvent import Solvent

        original = Solvent.from_json(_solvent("eta.json"))
        out_path = tmp_path / "eta_copy.json"
        original.to_json(out_path)

        reloaded = Solvent.from_json(out_path)
        assert reloaded.name == original.name
        assert len(reloaded.probes) == len(original.probes)

        for orig_p, rel_p in zip(original.probes, reloaded.probes):
            assert orig_p.name == rel_p.name
            assert orig_p.residue == rel_p.residue
            assert orig_p.atoms == rel_p.atoms

    def test_solvent_library_includes_bundled_solvents(self):
        """SolventLibrary lists at least ETA, WAT, and MAM."""
        from pymdmix.core.solvent import SolventLibrary

        lib = SolventLibrary()
        names = set(lib.list_solvents())
        for expected in ("ETA", "WAT", "MAM"):
            assert expected in names, f"Solvent '{expected}' not in library"

    def test_probe_probability_derived_from_volume(self):
        """Probe probability in eta.json is consistent with the solvent volume."""
        from pymdmix.core.solvent import Solvent

        solvent = Solvent.from_json(_solvent("eta.json"))
        assert solvent.volume is not None and solvent.volume > 0

        for probe in solvent.probes:
            # calculate_probability should return a value close to the stored one
            calc = solvent.calculate_probability(probe.name)
            if probe.probability is not None and calc is not None:
                assert abs(calc - probe.probability) < 1e-6, (
                    f"{probe.name}: stored={probe.probability}, calculated={calc}"
                )


# ---------------------------------------------------------------------------
# 6.  Grid I/O with DX files created from real structure data
# ---------------------------------------------------------------------------


class TestGridWithRealStructureData:
    """Create and verify DX grids derived from the real pep.pdb geometry."""

    def test_dx_written_from_real_structure_reloads(self, tmp_path):
        """DX grid created from the peptide bounding box round-trips correctly."""
        from pymdmix.core.grid import Grid
        from pymdmix.core.structure import load_structure

        struct = load_structure(_pep("pep.pdb"))
        coords = np.array([[a.xx, a.xy, a.xz] for a in struct.atoms])

        grid = Grid.from_structure(coords, spacing=1.0, padding=4.0)

        # Write and reload
        dx_path = tmp_path / "pep_grid.dx"
        grid.write_dx(dx_path)

        reloaded = Grid.read_dx(dx_path)
        assert reloaded.data.shape == grid.data.shape
        np.testing.assert_allclose(reloaded.origin, grid.origin, atol=1e-4)
        np.testing.assert_allclose(reloaded.spacing, grid.spacing, rtol=1e-5)

    def test_free_energy_grid_from_real_structure(self, tmp_path):
        """FE grid from the peptide structure has finite values and correct shape."""
        from pymdmix.analysis.energy import density_to_free_energy
        from pymdmix.core.grid import Grid
        from pymdmix.core.structure import load_structure

        struct = load_structure(_pep("pep.pdb"))
        coords = np.array([[a.xx, a.xy, a.xz] for a in struct.atoms])

        grid = Grid.from_structure(coords, spacing=1.0, padding=3.0)

        # Populate with synthetic density
        rng = np.random.default_rng(10)
        grid.data[:] = rng.uniform(0.5, 2.0, grid.data.shape)

        fe_grid = density_to_free_energy(grid, temperature=300.0, reference_density=1.0)

        assert fe_grid.data.shape == grid.data.shape
        assert np.isfinite(fe_grid.data).all()

        fe_path = tmp_path / "pep_fe.dx"
        fe_grid.write_dx(fe_path)
        assert fe_path.exists()

    def test_grid_add_bulk_counts_from_peptide_frames(self, tmp_path):
        """add_counts_bulk accumulates counts from a synthetic trajectory."""
        from pymdmix.core.grid import Grid
        from pymdmix.core.structure import load_structure

        struct = load_structure(_pep("pep.pdb"))
        coords = np.array([[a.xx, a.xy, a.xz] for a in struct.atoms])

        grid = Grid.from_structure(coords, spacing=1.0, padding=4.0)

        rng = np.random.default_rng(5)
        n_frames = 20
        for _ in range(n_frames):
            frame_coords = coords + rng.standard_normal(coords.shape) * 0.5
            grid.add_counts_bulk(frame_coords)

        density_grid = grid.to_density(n_frames)

        # Density should be concentrated around the peptide core
        max_density = density_grid.data.max()
        mean_density = density_grid.data[density_grid.data > 0].mean()
        assert max_density > mean_density, "Peak density should exceed mean density"


# ---------------------------------------------------------------------------
# 7.  CLI integration – using real data files
# ---------------------------------------------------------------------------


class TestCLIWithRealData:
    """CLI commands exercised against the bundled pep data."""

    @pytest.fixture
    def runner(self):
        from click.testing import CliRunner

        return CliRunner()

    def test_cli_info_solvents_shows_eta(self, runner):
        """``pymdmix info --solvents`` lists ETA from the bundled JSON files."""
        from pymdmix.cli import cli

        result = runner.invoke(cli, ["info", "--solvents"])
        assert result.exit_code == 0
        assert "ETA" in result.output

    def test_cli_tools_energy_on_synthetic_dx(self, runner, tmp_path):
        """``pymdmix tools energy`` converts a DX density grid to free energy."""
        from pymdmix.cli import cli
        from pymdmix.core.grid import Grid

        # Create a simple density grid
        density = Grid(
            data=np.ones((8, 8, 8)) * 2.5,
            origin=np.zeros(3),
            spacing=1.0,
        )
        in_dx = tmp_path / "density.dx"
        out_dx = tmp_path / "energy.dx"
        density.write_dx(in_dx)

        result = runner.invoke(
            cli,
            ["tools", "energy", "-i", str(in_dx), "-o", str(out_dx), "-T", "300"],
        )

        assert result.exit_code == 0
        assert out_dx.exists()

        fe_grid = Grid.read_dx(out_dx)
        assert fe_grid.data.shape == density.data.shape

    def test_cli_tools_sumgrids_with_real_like_grids(self, runner, tmp_path):
        """``pymdmix tools sumgrids`` adds two grids created from pep.pdb geometry."""
        from pymdmix.cli import cli
        from pymdmix.core.grid import Grid
        from pymdmix.core.structure import load_structure

        struct = load_structure(_pep("pep.pdb"))
        coords = np.array([[a.xx, a.xy, a.xz] for a in struct.atoms])

        g1 = Grid.from_structure(coords, spacing=1.0, padding=3.0)
        g1.data[:] = np.ones_like(g1.data)
        g2 = Grid.from_structure(coords, spacing=1.0, padding=3.0)
        g2.data[:] = np.ones_like(g2.data) * 2.0

        p1 = tmp_path / "g1.dx"
        p2 = tmp_path / "g2.dx"
        out = tmp_path / "sum.dx"
        g1.write_dx(p1)
        g2.write_dx(p2)

        result = runner.invoke(
            cli,
            ["tools", "sumgrids", "-g1", str(p1), "-g2", str(p2), "-o", str(out)],
        )

        assert result.exit_code == 0, result.output
        assert out.exists()

        total = Grid.read_dx(out)
        np.testing.assert_allclose(total.data, 3.0, atol=1e-5)

    def test_cli_create_project_and_check_json(self, runner, tmp_path):
        """``pymdmix create project`` creates a project.json with correct fields."""
        from pymdmix.cli import cli

        with runner.isolated_filesystem(temp_dir=tmp_path):
            result = runner.invoke(cli, ["create", "project", "-n", "pep_study"])
            assert result.exit_code == 0, result.output
            assert Path("pep_study").exists()

            proj_json = Path("pep_study") / "project.json"
            assert proj_json.exists()

            with open(proj_json) as fh:
                data = json.load(fh)

            assert data["name"] == "pep_study"
            assert "replicas" in data
