"""End-to-end CLI tests."""
import tempfile
from pathlib import Path

import numpy as np
import pytest
from click.testing import CliRunner

from pymdmix.cli import cli
from pymdmix.core.grid import Grid


@pytest.fixture
def runner():
    """Create CLI test runner."""
    return CliRunner()


@pytest.fixture
def temp_dir():
    """Create a temporary directory."""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield Path(tmpdir)


class TestCLIBasic:
    """Basic CLI tests."""
    
    def test_version(self, runner):
        """Test version command."""
        result = runner.invoke(cli, ["--version"])
        assert result.exit_code == 0
        assert "pyMDMix" in result.output or "0." in result.output
    
    def test_help(self, runner):
        """Test help command."""
        result = runner.invoke(cli, ["--help"])
        assert result.exit_code == 0
        assert "pyMDMix" in result.output
        assert "Commands:" in result.output


class TestInfoCommand:
    """Tests for info command."""
    
    def test_info_solvents(self, runner):
        """Test listing solvents."""
        result = runner.invoke(cli, ["info", "--solvents"])
        assert result.exit_code == 0
        assert "ETA" in result.output
        assert "WAT" in result.output


class TestCreateCommand:
    """Tests for create command."""
    
    def test_create_project(self, runner, temp_dir):
        """Test creating a project."""
        with runner.isolated_filesystem(temp_dir=temp_dir):
            result = runner.invoke(cli, ["create", "project", "-n", "test_proj"])
            assert result.exit_code == 0
            assert "created successfully" in result.output
            assert Path("test_proj").exists()
            assert Path("test_proj/project.json").exists()


class TestToolsCommands:
    """Tests for tools commands."""
    
    @pytest.fixture
    def sample_grids(self, temp_dir):
        """Create sample grids for testing."""
        grid1 = Grid(
            data=np.ones((10, 10, 10)) * 2.0,
            origin=(0.0, 0.0, 0.0),
            spacing=1.0,
        )
        grid2 = Grid(
            data=np.ones((10, 10, 10)) * 1.0,
            origin=(0.0, 0.0, 0.0),
            spacing=1.0,
        )
        
        path1 = temp_dir / "grid1.dx"
        path2 = temp_dir / "grid2.dx"
        grid1.write_dx(path1)
        grid2.write_dx(path2)
        
        return path1, path2
    
    def test_tools_diffgrids(self, runner, sample_grids, temp_dir):
        """Test grid subtraction."""
        g1, g2 = sample_grids
        output = temp_dir / "diff.dx"
        
        result = runner.invoke(cli, [
            "tools", "diffgrids",
            "-g1", str(g1),
            "-g2", str(g2),
            "-o", str(output),
        ])
        
        assert result.exit_code == 0
        assert output.exists()
        
        # Verify result: 2.0 - 1.0 = 1.0
        diff = Grid.read_dx(output)
        assert np.allclose(diff.data, 1.0)
    
    def test_tools_sumgrids(self, runner, sample_grids, temp_dir):
        """Test grid addition."""
        g1, g2 = sample_grids
        output = temp_dir / "sum.dx"
        
        result = runner.invoke(cli, [
            "tools", "sumgrids",
            "-g1", str(g1),
            "-g2", str(g2),
            "-o", str(output),
        ])
        
        assert result.exit_code == 0
        assert output.exists()
        
        # Verify result: 2.0 + 1.0 = 3.0
        total = Grid.read_dx(output)
        assert np.allclose(total.data, 3.0)
    
    def test_tools_avggrids(self, runner, sample_grids, temp_dir):
        """Test grid averaging."""
        g1, g2 = sample_grids
        output = temp_dir / "avg.dx"
        
        result = runner.invoke(cli, [
            "tools", "avggrids",
            "-i", str(g1),
            "-i", str(g2),
            "-o", str(output),
        ])
        
        assert result.exit_code == 0
        assert output.exists()
        
        # Verify result: (2.0 + 1.0) / 2 = 1.5
        avg = Grid.read_dx(output)
        assert np.allclose(avg.data, 1.5)
    
    def test_tools_energy(self, runner, temp_dir):
        """Test density to energy conversion."""
        # Create density grid
        density = Grid(
            data=np.ones((5, 5, 5)) * 2.0,  # 2x bulk density
            origin=(0.0, 0.0, 0.0),
            spacing=1.0,
        )
        input_path = temp_dir / "density.dx"
        output_path = temp_dir / "energy.dx"
        density.write_dx(input_path)
        
        result = runner.invoke(cli, [
            "tools", "energy",
            "-i", str(input_path),
            "-o", str(output_path),
            "-T", "300",
        ])
        
        assert result.exit_code == 0
        assert output_path.exists()
        assert "Converted to free energy" in result.output


class TestSetupCommands:
    """Tests for setup commands."""
    
    def test_setup_prepare_help(self, runner):
        """Test prepare help."""
        result = runner.invoke(cli, ["setup", "prepare", "--help"])
        assert result.exit_code == 0
        assert "cap" in result.output.lower() or "prepare" in result.output.lower()
    
    def test_setup_solvate_help(self, runner):
        """Test solvate help."""
        result = runner.invoke(cli, ["setup", "solvate", "--help"])
        assert result.exit_code == 0
        assert "--solvent" in result.output


class TestAnalyzeCommands:
    """Tests for analyze commands."""
    
    def test_analyze_density_help(self, runner):
        """Test density analysis help."""
        result = runner.invoke(cli, ["analyze", "density", "--help"])
        assert result.exit_code == 0
        assert "probe density" in result.output.lower()
    
    def test_analyze_hotspots_help(self, runner):
        """Test hotspots analysis help."""
        result = runner.invoke(cli, ["analyze", "hotspots", "--help"])
        assert result.exit_code == 0
        assert "hotspots" in result.output.lower()


class TestQueueCommand:
    """Tests for queue command."""
    
    def test_queue_help(self, runner):
        """Test queue help."""
        result = runner.invoke(cli, ["queue", "--help"])
        assert result.exit_code == 0
        assert "slurm" in result.output.lower()


class TestGridGetValue:
    """Tests for Grid.get_value method."""
    
    def test_get_value_inside(self):
        """Test getting value inside grid."""
        grid = Grid(
            data=np.arange(27).reshape(3, 3, 3).astype(float),
            origin=(0.0, 0.0, 0.0),
            spacing=1.0,
        )
        
        # Center of first voxel
        val = grid.get_value(0.5, 0.5, 0.5)
        assert val == 0.0
        
        # Second voxel along x
        val = grid.get_value(1.5, 0.5, 0.5)
        assert val == 9.0  # data[1,0,0]
    
    def test_get_value_outside(self):
        """Test getting value outside grid returns None."""
        grid = Grid(
            data=np.ones((3, 3, 3)),
            origin=(0.0, 0.0, 0.0),
            spacing=1.0,
        )
        
        assert grid.get_value(-1.0, 0.0, 0.0) is None
        assert grid.get_value(10.0, 0.0, 0.0) is None
