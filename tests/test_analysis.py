"""
Tests for pymdmix.analysis module.
"""

import pytest
import numpy as np
from pathlib import Path

from pymdmix.analysis.base import (
    Action,
    ActionResult,
    register_action,
    get_action,
    list_actions,
    run_action,
)
from pymdmix.analysis.density import DensityAction, calculate_density
from pymdmix.analysis.residence import ResidenceAction
from pymdmix.analysis.hotspots import HotspotAction, Hotspot
from pymdmix.core.grid import Grid


class TestActionResult:
    """Test ActionResult dataclass."""
    
    def test_success_result(self):
        """Test successful action result."""
        result = ActionResult(
            success=True,
            output_files=[Path("test.dx")],
            metadata={"n_frames": 100},
        )
        
        assert result.success is True
        assert len(result.output_files) == 1
        assert result.error is None
    
    def test_failure_result(self):
        """Test failed action result."""
        result = ActionResult(
            success=False,
            error="Something went wrong",
        )
        
        assert result.success is False
        assert result.error == "Something went wrong"
    
    def test_result_repr(self):
        """Test result string representation."""
        result = ActionResult(success=True, elapsed_time=1.5)
        
        repr_str = repr(result)
        assert "✓" in repr_str
        assert "1.50s" in repr_str


class TestActionRegistry:
    """Test action registration and discovery."""
    
    def test_register_action(self):
        """Test registering a custom action."""
        @register_action("test_custom")
        class CustomAction(Action):
            def run(self, trajectory, **kwargs):
                return ActionResult(success=True)
        
        assert "test_custom" in list_actions()
        assert get_action("test_custom") is CustomAction
    
    def test_get_unknown_action(self):
        """Test getting unknown action returns None."""
        assert get_action("nonexistent_action") is None
    
    def test_list_actions(self):
        """Test listing actions."""
        actions = list_actions()
        
        assert isinstance(actions, list)
        assert "density" in actions  # Should be registered
    
    def test_run_action_unknown(self, mock_trajectory):
        """Test running unknown action raises error."""
        with pytest.raises(ValueError, match="Unknown action"):
            run_action("nonexistent", mock_trajectory)


class TestDensityAction:
    """Test density calculation action."""
    
    def test_density_action_exists(self):
        """Test density action is registered."""
        assert get_action("density") is DensityAction
    
    def test_density_action_validation(self, mock_trajectory):
        """Test density action validation."""
        action = DensityAction()
        
        # Should fail without probes
        errors = action.validate(mock_trajectory)
        assert len(errors) > 0
        assert any("probe" in e.lower() for e in errors)
        
        # Should pass with probes
        errors = action.validate(
            mock_trajectory,
            probe_indices={"test": np.array([0, 1, 2])},
        )
        assert len(errors) == 0
    
    def test_density_calculation(self, mock_trajectory, sample_coordinates, tmp_output_dir):
        """Test full density calculation."""
        action = DensityAction()
        
        # Use first 10 atoms as probe
        probe_indices = {"test_probe": np.arange(10)}
        
        result = action.run(
            trajectory=mock_trajectory,
            reference=None,  # Will use first frame
            probe_indices=probe_indices,
            spacing=1.0,
            output_dir=tmp_output_dir,
        )
        
        assert result.success is True
        assert len(result.output_files) == 1
        assert result.metadata["n_frames"] == 10
        assert result.metadata["n_probes"] == 1
        
        # Check output file exists
        assert result.output_files[0].exists()
        assert result.output_files[0].suffix == ".dx"
    
    def test_density_with_free_energy(self, mock_trajectory, tmp_output_dir):
        """Test density calculation with free energy."""
        action = DensityAction()
        
        probe_indices = {"test": np.arange(10)}
        
        result = action.run(
            trajectory=mock_trajectory,
            probe_indices=probe_indices,
            spacing=1.0,
            output_dir=tmp_output_dir,
            compute_free_energy=True,
        )
        
        assert result.success is True
        assert len(result.output_files) == 2  # density + dg
        
        # Check both files exist
        density_file = [f for f in result.output_files if "density" in f.name][0]
        dg_file = [f for f in result.output_files if "dg" in f.name][0]
        
        assert density_file.exists()
        assert dg_file.exists()
    
    def test_density_multiple_probes(self, mock_trajectory, tmp_output_dir):
        """Test density calculation with multiple probes."""
        action = DensityAction()
        
        probe_indices = {
            "probe1": np.arange(0, 10),
            "probe2": np.arange(10, 20),
            "probe3": np.arange(20, 30),
        }
        
        result = action.run(
            trajectory=mock_trajectory,
            probe_indices=probe_indices,
            spacing=1.0,
            output_dir=tmp_output_dir,
        )
        
        assert result.success is True
        assert len(result.output_files) == 3
        assert result.metadata["n_probes"] == 3
    
    def test_density_no_probes(self, mock_trajectory, tmp_output_dir):
        """Test density calculation fails without probes."""
        action = DensityAction()
        
        result = action(
            trajectory=mock_trajectory,
            output_dir=tmp_output_dir,
        )
        
        assert result.success is False
        assert "probe" in result.error.lower()


class TestCalculateDensityFunction:
    """Test convenience function for density calculation."""
    
    def test_calculate_density(self, mock_trajectory, sample_coordinates):
        """Test calculate_density function."""
        # Get reference from first frame
        ref_coords = sample_coordinates
        
        probe_indices = {
            "probe1": np.arange(0, 10),
            "probe2": np.arange(10, 20),
        }
        
        grids = calculate_density(
            trajectory=mock_trajectory,
            probe_indices=probe_indices,
            reference_coords=ref_coords,
            spacing=1.0,
        )
        
        assert len(grids) == 2
        assert "probe1" in grids
        assert "probe2" in grids
        
        # Each grid should have counts
        for name, grid in grids.items():
            assert grid.data.sum() > 0  # Should have some density


class TestActionTiming:
    """Test action timing functionality."""
    
    def test_action_measures_time(self, mock_trajectory, tmp_output_dir):
        """Test that action measures elapsed time."""
        action = DensityAction()
        
        result = action(
            trajectory=mock_trajectory,
            probe_indices={"test": np.arange(10)},
            output_dir=tmp_output_dir,
        )
        
        assert result.elapsed_time > 0


class TestResidenceAction:
    """Test residence time analysis."""
    
    def test_residence_action_exists(self):
        """Test residence action is registered."""
        assert get_action("residence") is ResidenceAction
    
    def test_residence_validation(self, mock_trajectory):
        """Test residence action validation."""
        action = ResidenceAction()
        
        # Should fail without hotspot coords
        errors = action.validate(mock_trajectory)
        assert len(errors) > 0
        
        # Should pass with hotspot coords
        errors = action.validate(
            mock_trajectory,
            hotspot_coords=[(25.0, 25.0, 25.0)],
        )
        assert len(errors) == 0
    
    def test_residence_calculation(self, mock_trajectory, tmp_output_dir):
        """Test residence calculation."""
        action = ResidenceAction()
        
        # Define hotspot near center of coordinates
        hotspot_coords = [(25.0, 25.0, 25.0)]
        
        # Use all atoms, create mock residue IDs (one per atom for simplicity)
        n_atoms = mock_trajectory.n_atoms
        residue_ids = np.arange(n_atoms)
        
        result = action.run(
            trajectory=mock_trajectory,
            hotspot_coords=hotspot_coords,
            tolerance=15.0,  # Large tolerance to ensure hits
            residue_ids=residue_ids,
            output_dir=tmp_output_dir,
        )
        
        assert result.success is True
        assert len(result.output_files) == 2  # JSON + summary
        assert result.metadata["n_hotspots"] == 1
        assert result.metadata["n_frames"] == 10


class TestHotspotAction:
    """Test hotspot detection."""
    
    def test_hotspot_action_exists(self):
        """Test hotspot action is registered."""
        assert get_action("hotspots") is HotspotAction
    
    def test_hotspot_validation(self):
        """Test hotspot action validation."""
        action = HotspotAction()
        
        # Should fail without grids
        errors = action.validate()
        assert len(errors) > 0
        
        # Should pass with grids
        grid = Grid.from_bounds((0, 0, 0), (10, 10, 10), spacing=1.0)
        errors = action.validate(grids={"test": grid})
        assert len(errors) == 0
    
    def test_hotspot_detection(self, tmp_output_dir):
        """Test hotspot detection from grid."""
        # Create a density grid with a hotspot
        grid = Grid.from_bounds((0, 0, 0), (20, 20, 20), spacing=1.0)
        
        # Add high density at center (simulating a hotspot)
        grid.data[8:12, 8:12, 8:12] = 100.0  # High density region
        grid.data[grid.data == 0] = 1.0  # Background
        
        action = HotspotAction()
        result = action.run(
            grids={"test_probe": grid},
            energy_cutoff=-1.0,
            cluster_distance=3.0,
            min_points=2,
            output_dir=tmp_output_dir,
        )
        
        assert result.success is True
        assert len(result.output_files) == 3  # JSON + PDB + summary
        assert result.metadata["n_hotspots"] >= 1
    
    def test_hotspot_dataclass(self):
        """Test Hotspot dataclass."""
        coords = np.array([[10.0, 10.0, 10.0], [11.0, 10.0, 10.0]])
        energies = np.array([-1.0, -0.8])
        
        hotspot = Hotspot(
            id=0,
            probe="OH",
            centroid=(10.0, 10.0, 10.0),
            energy=-0.9,
            volume=2.0,
            n_points=2,
            coords=coords,
            energies=energies,
        )
        
        assert hotspot.min_energy == -1.0
        assert hotspot.mean_energy == -0.9
        assert hotspot.extent == (1.0, 0.0, 0.0)
        
        # Test serialization
        data = hotspot.to_dict()
        assert data["id"] == 0
        assert data["probe"] == "OH"
