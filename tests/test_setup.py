"""
Tests for pymdmix.setup module.
"""

from pymdmix.core.solvent import Solvent
from pymdmix.core.structure import load_structure
from pymdmix.setup.prepare import (
    PrepareResult,
    find_and_fix_disulfides,
    prepare_structure,
    renumber_residues,
)
from pymdmix.setup.solvate import (
    BoxConfig,
    IonConfig,
    SolvateResult,
    SolvationOptions,
    generate_leap_script,
)


class TestPrepareResult:
    """Test PrepareResult dataclass."""

    def test_default_result(self):
        """Test default prepare result."""
        import parmed

        struct = parmed.Structure()
        result = PrepareResult(structure=struct)

        assert result.n_caps_added == 0
        assert result.n_disulfides == 0
        assert result.warnings == []


class TestPrepareStructure:
    """Test structure preparation."""

    def test_prepare_simple_structure(self, sample_pdb_file):
        """Test preparing a simple structure."""
        struct = load_structure(sample_pdb_file)
        result = prepare_structure(struct)

        assert isinstance(result, PrepareResult)
        assert result.structure is not None

    def test_prepare_remove_waters(self, sample_pdb_file):
        """Test removing waters during preparation."""
        struct = load_structure(sample_pdb_file)

        # Sample has 3 water atoms
        result = prepare_structure(struct, remove_waters=True)

        assert result.n_waters_removed == 3

    def test_prepare_with_options(self, sample_pdb_file):
        """Test preparation with various options."""
        struct = load_structure(sample_pdb_file)

        result = prepare_structure(
            struct,
            cap_termini=True,
            fix_disulfides=True,
            remove_waters=False,
        )

        assert isinstance(result, PrepareResult)


class TestDisulfides:
    """Test disulfide handling."""

    def test_find_no_disulfides(self, sample_pdb_file):
        """Test finding disulfides when none present."""
        struct = load_structure(sample_pdb_file)
        disulfides = find_and_fix_disulfides(struct)

        # Sample PDB has no CYS
        assert disulfides == []


class TestRenumberResidues:
    """Test residue renumbering."""

    def test_renumber_from_one(self, sample_pdb_file):
        """Test renumbering residues from 1."""
        struct = load_structure(sample_pdb_file)
        renumber_residues(struct, start=1)

        assert struct.residues[0].number == 1
        assert struct.residues[1].number == 2

    def test_renumber_from_custom(self, sample_pdb_file):
        """Test renumbering from custom start."""
        struct = load_structure(sample_pdb_file)
        renumber_residues(struct, start=100)

        assert struct.residues[0].number == 100


class TestBoxConfig:
    """Test BoxConfig dataclass."""

    def test_default_config(self):
        """Test default box configuration."""
        config = BoxConfig()

        assert config.shape == "truncated_octahedron"
        assert config.buffer == 12.0

    def test_leap_command_oct(self):
        """Test LEaP command for truncated octahedron."""
        config = BoxConfig(shape="truncated_octahedron", buffer=10.0)
        cmd = config.to_leap_command("sys", "ETABOX")

        assert "solvateoct" in cmd
        assert "sys" in cmd
        assert "ETABOX" in cmd
        assert "10.0" in cmd

    def test_leap_command_box(self):
        """Test LEaP command for rectangular box."""
        config = BoxConfig(shape="box", buffer=15.0)
        cmd = config.to_leap_command("sys", "WATBOX")

        assert "solvatebox" in cmd


class TestIonConfig:
    """Test IonConfig dataclass."""

    def test_default_config(self):
        """Test default ion configuration."""
        config = IonConfig()

        assert config.neutralize is True
        assert config.cation == "Na+"
        assert config.anion == "Cl-"

    def test_leap_commands(self):
        """Test LEaP command generation."""
        config = IonConfig(neutralize=True)
        commands = config.to_leap_commands("sys")

        assert len(commands) >= 2
        assert any("Na+" in cmd for cmd in commands)


class TestLeapScript:
    """Test LEaP script generation."""

    def test_generate_basic_script(self, tmp_output_dir):
        """Test generating a basic LEaP script."""
        pdb_path = tmp_output_dir / "test.pdb"
        pdb_path.write_text("ATOM      1  CA  ALA A   1       0.0   0.0   0.0  1.0  0.0\nEND\n")

        solvent = Solvent(
            name="ETA",
            description="Ethanol",
        )

        script = generate_leap_script(
            input_path=pdb_path,
            solvent=solvent,
            output_prefix="system",
        )

        assert "source leaprc" in script
        assert "loadpdb" in script
        assert "saveamberparm" in script
        assert "quit" in script

    def test_script_with_custom_ff(self, tmp_output_dir):
        """Test script with custom force fields."""
        pdb_path = tmp_output_dir / "test.pdb"
        pdb_path.write_text("END\n")

        solvent = Solvent(name="WAT")

        script = generate_leap_script(
            input_path=pdb_path,
            solvent=solvent,
            output_prefix="system",
            force_fields=["leaprc.protein.ff14SB", "leaprc.water.tip3p"],
        )

        assert "ff14SB" in script
        assert "tip3p" in script


class TestSolvationOptions:
    """Test SolvationOptions dataclass."""

    def test_defaults(self):
        """Test default values."""
        opts = SolvationOptions()

        assert opts.box_buffer == 12.0
        assert opts.box_shape == "truncated_octahedron"
        assert opts.neutralize is True
        assert opts.ion_concentration == 0.15
        assert opts.cation == "Na+"
        assert opts.anion == "Cl-"
        assert opts.extra_forcefields == []

    def test_custom_values(self):
        """Test creating with custom values."""
        opts = SolvationOptions(
            box_buffer=15.0,
            box_shape="box",
            neutralize=False,
            ion_concentration=0.0,
        )

        assert opts.box_buffer == 15.0
        assert opts.box_shape == "box"
        assert opts.neutralize is False

    def test_extra_forcefields(self):
        """Test extra_forcefields list."""
        opts = SolvationOptions(extra_forcefields=["leaprc.gaff2"])

        assert "leaprc.gaff2" in opts.extra_forcefields


class TestSolvateResult:
    """Test SolvateResult dataclass."""

    def test_default_result(self):
        """Test default solvate result."""
        result = SolvateResult()

        assert result.success is False
        assert result.topology is None
        assert result.error is None

    def test_successful_result(self, tmp_output_dir):
        """Test successful solvation result."""
        result = SolvateResult(
            success=True,
            topology=tmp_output_dir / "system.prmtop",
            coordinates=tmp_output_dir / "system.inpcrd",
            n_solvent_residues=1000,
        )

        assert result.success is True
        assert result.n_solvent_residues == 1000

    def test_save_coordinates(self, tmp_output_dir):
        """Test save_coordinates copies the file."""
        src = tmp_output_dir / "system.inpcrd"
        src.write_text("dummy crd content")
        dest = tmp_output_dir / "output.inpcrd"

        result = SolvateResult(success=True, coordinates=src)
        result.save_coordinates(dest)

        assert dest.exists()
        assert dest.read_text() == "dummy crd content"

    def test_save_coordinates_no_file_raises(self):
        """save_coordinates raises ValueError when no coordinates."""
        import pytest

        result = SolvateResult(success=False)
        with pytest.raises(ValueError, match="No coordinates"):
            result.save_coordinates("/tmp/nowhere.inpcrd")

    def test_save_topology(self, tmp_output_dir):
        """Test save_topology copies the topology file."""
        src = tmp_output_dir / "system.prmtop"
        src.write_text("dummy prmtop content")
        dest = tmp_output_dir / "output.prmtop"

        result = SolvateResult(success=True, topology=src)
        result.save_topology(dest)

        assert dest.exists()
        assert dest.read_text() == "dummy prmtop content"

    def test_save_topology_no_file_raises(self):
        """save_topology raises ValueError when no topology."""
        import pytest

        result = SolvateResult(success=False)
        with pytest.raises(ValueError, match="No topology"):
            result.save_topology("/tmp/nowhere.prmtop")


# =============================================================================
# Solvate - fixed behaviours matching the original pyMDMix workflow
# =============================================================================


class TestLeapScriptFixed:
    """Tests covering fixes to generate_leap_script to match original pyMDMix."""

    @staticmethod
    def _dummy_off(path, stem="prot"):
        """Create a minimal placeholder OFF file and return its Path."""
        off = path / f"{stem}.off"
        off.write_text(f"!entry.{stem}.unit.name single str\n \"{stem}\"\n")
        return off

    def test_solvent_box_uses_box_unit_not_nameBox(self, tmp_output_dir):
        """The tleap solvation command uses solvent.box_unit, not '{name}BOX'.

        The original code solvated with the unit name from the OFF file
        (e.g. 'ETAWAT20'), not a made-up 'ETABOX' name.
        """
        pdb = tmp_output_dir / "prot.pdb"
        pdb.write_text("END\n")

        solvent = Solvent(name="ETA", box_unit="ETAWAT20")
        script = generate_leap_script(input_path=pdb, solvent=solvent, output_prefix="out")

        # Must contain the actual box unit name
        assert "ETAWAT20" in script
        # Must NOT contain the wrong made-up name
        assert "ETABOX" not in script

    def test_solvent_box_fallback_when_box_unit_empty(self, tmp_output_dir):
        """Falls back to '{name}BOX' when box_unit is not set."""
        pdb = tmp_output_dir / "prot.pdb"
        pdb.write_text("END\n")

        solvent = Solvent(name="ETA", box_unit="")
        script = generate_leap_script(input_path=pdb, solvent=solvent, output_prefix="out")

        assert "ETABOX" in script

    def test_off_input_uses_loadoff_not_loadpdb(self, tmp_output_dir):
        """When input is an OFF file, script uses loadoff + copy, not loadpdb.

        The original workflow's primary input was an Amber Object File
        (pre-parameterised protein unit), not a PDB.
        """
        off = self._dummy_off(tmp_output_dir)

        solvent = Solvent(name="ETA", box_unit="ETAWAT20")
        script = generate_leap_script(
            input_path=off,
            solvent=solvent,
            output_prefix="out",
            unit_name="prot",
        )

        assert "loadoff" in script
        assert f"loadoff {off}" in script
        assert "sys = copy prot" in script
        # Must NOT use loadpdb for the protein
        # (it may use loadoff for the solvent, but loadpdb must not appear)
        assert "loadpdb" not in script

    def test_off_input_uses_stem_as_unit_when_unit_name_not_given(self, tmp_output_dir):
        """When unit_name is None, the file stem is used as the unit name."""
        off = self._dummy_off(tmp_output_dir, stem="myprotein")

        solvent = Solvent(name="ETA", box_unit="ETAWAT20")
        script = generate_leap_script(
            input_path=off,
            solvent=solvent,
            output_prefix="out",
        )

        assert "sys = copy myprotein" in script

    def test_pdb_input_still_uses_loadpdb(self, tmp_output_dir):
        """PDB input still generates the loadpdb command (unchanged behaviour)."""
        pdb = tmp_output_dir / "prot.pdb"
        pdb.write_text("END\n")

        solvent = Solvent(name="ETA", box_unit="ETAWAT20")
        script = generate_leap_script(input_path=pdb, solvent=solvent, output_prefix="out")

        assert f"sys = loadpdb {pdb}" in script

    def test_extra_forcefields_threaded_into_leap_script(self, tmp_output_dir):
        """extra_forcefields from SolvationOptions reach the leap script."""
        from pymdmix.setup.solvate import SolvateResult, SolvationOptions, solvate_structure
        import unittest.mock as mock

        pdb = tmp_output_dir / "prot.pdb"
        pdb.write_text("END\n")

        solvent = Solvent(name="ETA", box_unit="ETAWAT20")

        # Capture the script passed to run_leap
        captured = {}

        def fake_run_leap(script, leap_exe, work_dir):
            captured["script"] = script
            return False, "mocked"  # Simulate failure so we just inspect script

        with mock.patch("pymdmix.setup.solvate.run_leap", side_effect=fake_run_leap):
            solvate_structure(
                pdb,
                solvent,
                output_dir=tmp_output_dir / "out",
                output_prefix="out",
                options=SolvationOptions(extra_forcefields=["leaprc.gaff2"]),
            )

        assert "leaprc.gaff2" in captured.get("script", "")

    def test_default_forcefields_use_ff14sb_and_gaff(self, tmp_output_dir):
        """Default force fields match original pyMDMix defaults (ff14SB + gaff)."""
        pdb = tmp_output_dir / "prot.pdb"
        pdb.write_text("END\n")

        solvent = Solvent(name="ETA", box_unit="ETAWAT20")
        script = generate_leap_script(input_path=pdb, solvent=solvent, output_prefix="out")

        assert "ff14SB" in script
        assert "gaff" in script
