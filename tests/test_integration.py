"""
Integration tests for pyMDMix.

Tests full workflows from structure loading through analysis.
"""
import tempfile
from pathlib import Path

import numpy as np
import pytest

from pymdmix.core.grid import Grid
from pymdmix.core.solvent import Solvent, SolventLibrary, Probe
from pymdmix.core.trajectory import Frame
from pymdmix.project.config import Config
from pymdmix.project.replica import Replica, ReplicaState
from pymdmix.project.project import Project
from pymdmix.analysis.base import get_action, list_actions
from pymdmix.analysis.density import DensityAction


# =============================================================================
# Fixtures
# =============================================================================

@pytest.fixture
def temp_dir():
    """Create a temporary directory for test files."""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield Path(tmpdir)


@pytest.fixture
def sample_solvent():
    """Create a sample solvent for testing."""
    return Solvent(
        name="TEST",
        description="Test solvent mixture",
        probes=[
            Probe("WAT", "WAT", ["O"], "Water oxygen"),
            Probe("CT", "TST", ["C1"], "Test carbon"),
        ],
    )


@pytest.fixture
def sample_grid():
    """Create a sample grid for testing."""
    return Grid(
        data=np.zeros((20, 20, 20)),
        origin=(0.0, 0.0, 0.0),
        spacing=1.0,
    )


@pytest.fixture
def mock_trajectory():
    """Create a mock trajectory with predictable coordinates."""
    class MockTrajectory:
        def __init__(self):
            self.n_frames = 10
            self.n_atoms = 100
            
        def __len__(self):
            return self.n_frames
        
        def __iter__(self):
            rng = np.random.RandomState(42)
            for i in range(self.n_frames):
                # Generate coordinates in a 20x20x20 box
                coords = rng.uniform(0, 20, (self.n_atoms, 3))
                yield Frame(
                    coordinates=coords,
                    time=i * 1.0,
                    box=np.array([20.0, 20.0, 20.0, 90.0, 90.0, 90.0]),
                )
    
    return MockTrajectory()


# =============================================================================
# Workflow Tests
# =============================================================================

class TestProjectWorkflow:
    """Test complete project creation and management workflow."""
    
    def test_create_project(self, temp_dir):
        """Test creating a project."""
        project = Project(
            name="test_project",
            path=temp_dir / "test_project",
        )
        
        assert project.name == "test_project"
        assert project.n_replicas == 0
    
    def test_add_replica_to_project(self, temp_dir, sample_solvent):
        """Test adding replicas to a project."""
        project = Project(
            name="test_project",
            path=temp_dir / "test_project",
        )
        
        replica = Replica(
            name="rep1",
            solvent=sample_solvent.name,
            path=temp_dir / "test_project" / "rep1",
        )
        
        project.add_replica(replica)
        
        assert project.n_replicas == 1
        assert project.get_replica("rep1") == replica
    
    def test_project_save_load_roundtrip(self, temp_dir):
        """Test saving and loading a project."""
        # Create project
        project = Project(
            name="roundtrip_test",
            path=temp_dir / "roundtrip_test",
        )
        project.path.mkdir(parents=True, exist_ok=True)
        
        replica = Replica(
            name="rep1",
            solvent="MAM",
            path=project.path / "rep1",
        )
        replica.path.mkdir(parents=True, exist_ok=True)
        project.add_replica(replica)
        
        # Save
        project.save()
        
        # Load
        loaded = Project.load(project.path)
        
        assert loaded.name == project.name
        assert loaded.n_replicas == 1


class TestAnalysisWorkflow:
    """Test complete analysis workflow."""
    
    def test_density_analysis_workflow(self, mock_trajectory, sample_grid):
        """Test running density analysis on a trajectory."""
        grid = sample_grid
        
        # Create atom mask (first 10 atoms are "probe atoms")
        probe_mask = np.zeros(mock_trajectory.n_atoms, dtype=bool)
        probe_mask[:10] = True
        
        # Run density calculation
        for frame in mock_trajectory:
            coords = frame.coordinates[probe_mask]
            grid.add_counts_bulk(coords)
        
        # Check results
        assert grid.data.sum() > 0
        assert grid.data.max() > 0
    
    def test_action_registry(self):
        """Test that all expected actions are registered."""
        actions = list_actions()
        
        assert "density" in actions
        assert "residence" in actions
        assert "hotspots" in actions
    
    def test_get_action_by_name(self):
        """Test retrieving actions by name."""
        density_action = get_action("density")
        
        assert density_action is not None
        assert density_action.name == "density"


class TestSolventLibraryWorkflow:
    """Test solvent library operations."""
    
    def test_library_loads_all_solvents(self):
        """Test that the library loads all migrated solvents."""
        library = SolventLibrary()
        solvents = library.list_solvents()
        
        # Should have at least the migrated solvents
        expected = {"ANT", "ETA", "ION", "ISO", "ISO5", "MAM", "MOH", "PYR1", "WAT"}
        assert expected.issubset(set(solvents))
    
    def test_solvent_has_probes(self):
        """Test that loaded solvents have probe definitions."""
        library = SolventLibrary()
        eta = library.get("ETA")
        
        assert eta is not None
        assert len(eta.probes) > 0
        
        # Check probe has selection string
        probe = eta.probes[0]
        assert probe.selection  # MDAnalysis selection string
    
    def test_solvent_json_roundtrip(self, temp_dir):
        """Test saving and loading a solvent as JSON."""
        original = Solvent(
            name="CUSTOM",
            description="Custom test solvent",
            probes=[
                Probe("P1", "RES", ["A1", "A2"], "Test probe"),
            ],
        )
        
        # Save
        json_path = temp_dir / "custom.json"
        original.to_json(json_path)
        
        # Load
        loaded = Solvent.from_json(json_path)
        
        assert loaded.name == original.name
        assert len(loaded.probes) == 1
        assert loaded.probes[0].atoms == ["A1", "A2"]


class TestGridIOWorkflow:
    """Test grid I/O operations."""
    
    def test_grid_dx_roundtrip(self, temp_dir):
        """Test saving and loading a grid in DX format."""
        # Create grid with data
        grid = Grid(
            data=np.random.rand(10, 10, 10),
            origin=(0.0, 0.0, 0.0),
            spacing=0.5,
        )
        
        # Save
        dx_path = temp_dir / "test.dx"
        grid.write_dx(dx_path)
        
        # Load
        loaded = Grid.read_dx(dx_path)
        
        np.testing.assert_array_almost_equal(loaded.data, grid.data)
        assert loaded.origin == grid.origin
        assert loaded.spacing == grid.spacing
    
    def test_grid_free_energy_conversion(self):
        """Test converting density to free energy."""
        # Create grid with varying density (high spot in center)
        data = np.ones((5, 5, 5)) * 0.5  # Background
        data[2, 2, 2] = 5.0  # High density spot
        grid = Grid(
            data=data,
            origin=(0.0, 0.0, 0.0),
            spacing=1.0,
        )
        
        # Convert to free energy
        fe_grid = grid.to_free_energy(temperature=300.0)
        
        # High density point should have more negative free energy
        # Low density points should have positive free energy
        assert fe_grid.data[2, 2, 2] < fe_grid.data[0, 0, 0]


class TestCLIIntegration:
    """Test CLI commands work correctly."""
    
    def test_cli_imports(self):
        """Test that CLI module imports without errors."""
        from pymdmix.cli import cli, create, analyze, info
        
        assert cli is not None
        assert create is not None
    
    def test_cli_help(self):
        """Test CLI help runs without errors."""
        from click.testing import CliRunner
        from pymdmix.cli import cli
        
        runner = CliRunner()
        result = runner.invoke(cli, ["--help"])
        
        assert result.exit_code == 0
        assert "pyMDMix" in result.output


# =============================================================================
# End-to-End Test
# =============================================================================

class TestEndToEnd:
    """Full end-to-end workflow test."""
    
    def test_full_workflow(self, temp_dir):
        """
        Test complete workflow:
        1. Create project
        2. Add replica
        3. Create mock density grid
        4. Save outputs
        """
        # 1. Create project
        project = Project(
            name="e2e_test",
            path=temp_dir / "e2e_test",
        )
        project.path.mkdir(parents=True, exist_ok=True)
        
        # 2. Add replica
        replica = Replica(
            name="rep1",
            solvent="ETA",
            path=project.path / "rep1",
        )
        replica.path.mkdir(parents=True, exist_ok=True)
        project.add_replica(replica)
        
        # 3. Create density grid (simulating analysis result)
        rng = np.random.RandomState(42)
        grid = Grid(
            data=rng.exponential(scale=1.0, size=(20, 20, 20)),
            origin=(0.0, 0.0, 0.0),
            spacing=0.5,
        )
        
        # 4. Save outputs
        # Save grid
        grid_path = replica.path / "density_WAT.dx"
        grid.write_dx(grid_path)
        
        # Save project
        project.save()
        
        # Update replica state
        replica.state = ReplicaState.ANALYZED
        
        # Verify outputs
        assert grid_path.exists()
        assert (project.path / "project.json").exists()
        
        # Load and verify
        loaded_grid = Grid.read_dx(grid_path)
        assert loaded_grid.shape == (20, 20, 20)
        
        loaded_project = Project.load(project.path)
        assert loaded_project.name == "e2e_test"
        assert loaded_project.n_replicas == 1
