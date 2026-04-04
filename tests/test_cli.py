"""Tests for pymdmix CLI."""

from pathlib import Path

import pytest
from click.testing import CliRunner

from pymdmix.cli import cli


@pytest.fixture
def runner():
    """Create a CLI test runner."""
    return CliRunner()


@pytest.fixture
def temp_dir(tmp_path):
    """Provide a temporary directory."""
    return tmp_path


# =============================================================================
# Basic CLI Tests
# =============================================================================


class TestBasicCLI:
    def test_version(self, runner):
        """Test --version flag."""
        result = runner.invoke(cli, ["--version"])
        assert result.exit_code == 0
        assert "pyMDMix" in result.output

    def test_help(self, runner):
        """Test --help flag."""
        result = runner.invoke(cli, ["--help"])
        assert result.exit_code == 0
        assert "pyMDMix" in result.output
        assert "create" in result.output
        assert "analyze" in result.output


# =============================================================================
# Create Command Tests
# =============================================================================


class TestCreateCommands:
    def test_create_help(self, runner):
        """Test create subcommand help."""
        result = runner.invoke(cli, ["create", "--help"])
        assert result.exit_code == 0
        assert "project" in result.output
        assert "replica" in result.output

    def test_create_project(self, runner, temp_dir):
        """Test creating a project."""
        with runner.isolated_filesystem(temp_dir=temp_dir):
            result = runner.invoke(cli, ["create", "project", "-n", "test_project"])

            assert result.exit_code == 0
            assert "created successfully" in result.output

            # Check project directory was created
            project_dir = Path("test_project")
            assert project_dir.exists()

    def test_create_project_custom_dir(self, runner, temp_dir):
        """Test creating a project in custom directory."""
        with runner.isolated_filesystem(temp_dir=temp_dir):
            result = runner.invoke(
                cli, ["create", "project", "-n", "myproject", "-d", "custom_location"]
            )

            assert result.exit_code == 0
            assert Path("custom_location").exists()


# =============================================================================
# Setup Command Tests
# =============================================================================


class TestSetupCommands:
    def test_setup_help(self, runner):
        """Test setup subcommand help."""
        result = runner.invoke(cli, ["setup", "--help"])
        assert result.exit_code == 0
        assert "prepare" in result.output
        assert "solvate" in result.output

    def test_prepare_help(self, runner):
        """Test prepare command help."""
        result = runner.invoke(cli, ["setup", "prepare", "--help"])
        assert result.exit_code == 0
        assert "cap" in result.output
        assert "disulfide" in result.output


# =============================================================================
# Analyze Command Tests
# =============================================================================


class TestAnalyzeCommands:
    def test_analyze_help(self, runner):
        """Test analyze subcommand help."""
        result = runner.invoke(cli, ["analyze", "--help"])
        assert result.exit_code == 0
        assert "density" in result.output
        assert "hotspots" in result.output

    def test_density_help(self, runner):
        """Test density command help."""
        result = runner.invoke(cli, ["analyze", "density", "--help"])
        assert result.exit_code == 0
        assert "trajectory" in result.output
        assert "probe" in result.output

    def test_hotspots_help(self, runner):
        """Test hotspots command help."""
        result = runner.invoke(cli, ["analyze", "hotspots", "--help"])
        assert result.exit_code == 0
        assert "threshold" in result.output or "cutoff" in result.output
        assert "min-size" in result.output


# =============================================================================
# Info Command Tests
# =============================================================================


class TestInfoCommand:
    def test_info_solvents(self, runner):
        """Test listing solvents."""
        result = runner.invoke(cli, ["info", "--solvents"])

        assert result.exit_code == 0
        assert "Available solvents" in result.output

    def test_info_help(self, runner):
        """Test info command help."""
        result = runner.invoke(cli, ["info", "--help"])
        assert result.exit_code == 0
        assert "project" in result.output
        assert "solvents" in result.output


# =============================================================================
# Queue Command Tests
# =============================================================================


class TestQueueCommand:
    def test_queue_help(self, runner):
        """Test queue command help."""
        result = runner.invoke(cli, ["queue", "--help"])
        assert result.exit_code == 0
        assert "generate" in result.output.lower()
        assert "submit" in result.output.lower()


# =============================================================================
# Verbose Mode Tests
# =============================================================================


class TestVerboseMode:
    def test_verbose_flag(self, runner, temp_dir):
        """Test verbose flag passes through."""
        with runner.isolated_filesystem(temp_dir=temp_dir):
            result = runner.invoke(cli, ["-v", "create", "project", "-n", "test"])
            assert result.exit_code == 0


# =============================================================================
# Error Handling Tests
# =============================================================================


class TestErrorHandling:
    def test_missing_required_arg(self, runner):
        """Test error when required argument missing."""
        result = runner.invoke(cli, ["analyze", "density"])
        # Should fail or show help since neither --project nor --trajectory given
        assert result.exit_code != 0 or "Specify either" in result.output

    def test_file_not_found(self, runner):
        """Test error when file doesn't exist."""
        result = runner.invoke(cli, ["setup", "prepare", "/nonexistent/file.pdb"])
        assert result.exit_code != 0

    def test_unknown_command(self, runner):
        """Test error for unknown command."""
        result = runner.invoke(cli, ["foobar"])
        assert result.exit_code != 0


# =============================================================================
# Full Project Config Tests (create project -f project.cfg)
# =============================================================================


class TestCreateProjectWithConfig:
    """Test `create project -f project.cfg` for both PDB and OFF inputs."""

    def _write_project_cfg(self, path: Path, input_path: Path, unit_name: str | None = None) -> Path:
        """Write a minimal project config referencing an input file."""
        cfg = path / "project.cfg"
        suffix = input_path.suffix.lower()
        if suffix in (".off", ".lib"):
            input_key = "OFF"
        else:
            input_key = "PDB"
        uname_line = f"UNAME = {unit_name}\n" if unit_name else ""
        cfg.write_text(
            f"[SYSTEM]\n"
            f"NAME = prot\n"
            f"{input_key} = {input_path}\n"
            f"{uname_line}"
            f"\n"
            f"[MDSETTINGS]\n"
            f"SOLVENT = ETA\n"
            f"NREPL = 1\n"
        )
        return cfg

    def _write_dummy_pdb(self, path: Path) -> Path:
        """Write a minimal PDB file."""
        pdb = path / "prot.pdb"
        pdb.write_text(
            "ATOM      1  CA  ALA A   1       1.000   1.000   1.000  1.00  0.00           C  \n"
            "END\n"
        )
        return pdb

    def _write_dummy_off(self, path: Path) -> Path:
        """Write a minimal placeholder OFF file."""
        off = path / "prot.off"
        off.write_text("!entry.prot.unit.name single str\n \"prot\"\n")
        return off

    def test_create_project_cfg_no_tleap(self, runner, temp_dir):
        """create project with a .cfg file warns when tleap is absent."""
        import unittest.mock as mock

        with runner.isolated_filesystem(temp_dir=temp_dir):
            pdb = self._write_dummy_pdb(Path("."))
            cfg = self._write_project_cfg(Path("."), pdb)

            # Pretend tleap is not on PATH
            with mock.patch("pymdmix.cli.shutil.which", return_value=None):
                result = runner.invoke(
                    cli,
                    ["create", "project", "-n", "myproj", "-f", str(cfg)],
                )

            assert result.exit_code == 0, result.output
            assert "myproj" in result.output or "prot" in result.output

            # Project directory must exist
            assert Path("myproj").exists()

            # Replica directories must exist
            replicas_dir = Path("myproj") / "replicas"
            assert replicas_dir.exists()
            assert any(replicas_dir.iterdir())

            # Systems directory must exist (created even without solvation)
            assert (Path("myproj") / "systems").exists()

    def test_create_project_cfg_off_input_no_tleap(self, runner, temp_dir):
        """create project with an OFF input file warns when tleap is absent."""
        import unittest.mock as mock

        with runner.isolated_filesystem(temp_dir=temp_dir):
            off = self._write_dummy_off(Path("."))
            cfg = self._write_project_cfg(Path("."), off, unit_name="prot")

            with mock.patch("pymdmix.cli.shutil.which", return_value=None):
                result = runner.invoke(
                    cli,
                    ["create", "project", "-n", "myproj", "-f", str(cfg)],
                )

            assert result.exit_code == 0, result.output
            # Replicas must still be created
            assert any((Path("myproj") / "replicas").iterdir())

    def test_create_project_cfg_replica_has_input_files_when_tleap_available(
        self, runner, temp_dir
    ):
        """When tleap succeeds (PDB input), replicas get topology, input files, and COMMANDS.sh."""
        import unittest.mock as mock
        from pymdmix.setup.solvate import SolvateResult

        with runner.isolated_filesystem(temp_dir=temp_dir):
            pdb = self._write_dummy_pdb(Path("."))
            cfg = self._write_project_cfg(Path("."), pdb)

            # Mock tleap as available and solvation as successful
            fake_top = Path("fake.prmtop")
            fake_crd = Path("fake.inpcrd")
            fake_top.write_text("FAKE_PRMTOP")
            fake_crd.write_text("FAKE_INPCRD")

            fake_result = SolvateResult(
                topology=fake_top.resolve(),
                coordinates=fake_crd.resolve(),
                success=True,
            )

            with (
                mock.patch("pymdmix.cli.shutil.which", return_value="/usr/bin/tleap"),
                mock.patch(
                    "pymdmix.setup.solvate.solvate_structure",
                    return_value=fake_result,
                ),
            ):
                result = runner.invoke(
                    cli,
                    ["create", "project", "-n", "myproj", "-f", str(cfg)],
                )

            assert result.exit_code == 0, result.output

            # Each replica should have COMMANDS.sh and input files
            for rep_dir in (Path("myproj") / "replicas").iterdir():
                if rep_dir.is_dir():
                    assert (rep_dir / "COMMANDS.sh").exists(), f"COMMANDS.sh missing in {rep_dir}"
                    assert (rep_dir / "min" / "min.in").exists(), f"min.in missing in {rep_dir}"
                    assert (rep_dir / "eq" / "eq1.in").exists(), f"eq1.in missing in {rep_dir}"
                    assert (rep_dir / "md" / "prod.in").exists(), f"prod.in missing in {rep_dir}"

    def test_create_project_cfg_off_input_replica_files(self, runner, temp_dir):
        """When tleap succeeds with OFF input, replicas get input files and COMMANDS.sh."""
        import unittest.mock as mock
        from pymdmix.setup.solvate import SolvateResult

        with runner.isolated_filesystem(temp_dir=temp_dir):
            off = self._write_dummy_off(Path("."))
            cfg = self._write_project_cfg(Path("."), off, unit_name="prot")

            fake_top = Path("fake.prmtop")
            fake_crd = Path("fake.inpcrd")
            fake_top.write_text("FAKE_PRMTOP")
            fake_crd.write_text("FAKE_INPCRD")

            fake_result = SolvateResult(
                topology=fake_top.resolve(),
                coordinates=fake_crd.resolve(),
                success=True,
            )

            with (
                mock.patch("pymdmix.cli.shutil.which", return_value="/usr/bin/tleap"),
                mock.patch(
                    "pymdmix.setup.solvate.solvate_structure",
                    return_value=fake_result,
                ),
            ):
                result = runner.invoke(
                    cli,
                    ["create", "project", "-n", "myproj", "-f", str(cfg)],
                )

            assert result.exit_code == 0, result.output

            for rep_dir in (Path("myproj") / "replicas").iterdir():
                if rep_dir.is_dir():
                    assert (rep_dir / "COMMANDS.sh").exists(), f"COMMANDS.sh missing in {rep_dir}"
                    assert (rep_dir / "min" / "min.in").exists(), f"min.in missing in {rep_dir}"
