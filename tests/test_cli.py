"""Tests for pymdmix CLI."""
import json
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
            result = runner.invoke(cli, [
                "create", "project", 
                "-n", "myproject",
                "-d", "custom_location"
            ])
            
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
        assert "percentile" in result.output
        assert "cutoff" in result.output


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
        assert "slurm" in result.output
        assert "sge" in result.output
        assert "pbs" in result.output


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
        result = runner.invoke(cli, [
            "setup", "prepare", "/nonexistent/file.pdb"
        ])
        assert result.exit_code != 0
    
    def test_unknown_command(self, runner):
        """Test error for unknown command."""
        result = runner.invoke(cli, ["foobar"])
        assert result.exit_code != 0
