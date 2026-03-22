"""
Tests for pymdmix.analysis module.
"""

from pathlib import Path

import numpy as np
import pytest

from pymdmix.analysis.base import (
    Action,
    ActionResult,
    get_action,
    list_actions,
    register_action,
    run_action,
)
from pymdmix.analysis.density import (
    CpptrajDensityAction,
    DensityAction,
    DensityAllHAAction,
    DensityProteinAction,
    calculate_density,
)
from pymdmix.analysis.hotspots import Hotspot, HotspotAction
from pymdmix.analysis.residence import ResidenceAction, ResidenceResult
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


class TestDensityActionCOM:
    """Test center of mass (COM) probe features."""

    def test_density_with_include_com(self, mock_trajectory, tmp_output_dir):
        """Test density calculation with include_com=True."""
        action = DensityAction()

        # Regular probe indices
        probe_indices = {"test_probe": np.arange(10)}

        # COM residue indices: 2 residues, each with 5 atoms
        com_residue_indices = {
            "solvent": [np.arange(0, 5), np.arange(5, 10)],
        }

        result = action.run(
            trajectory=mock_trajectory,
            probe_indices=probe_indices,
            spacing=1.0,
            output_dir=tmp_output_dir,
            include_com=True,
            com_residue_indices=com_residue_indices,
        )

        assert result.success is True
        assert result.metadata["n_probes"] == 2  # test_probe + solvent_COM
        assert any("COM" in str(f) for f in result.output_files)

    def test_density_with_only_com(self, mock_trajectory, tmp_output_dir):
        """Test density calculation with only_com=True."""
        action = DensityAction()

        # COM residue indices
        com_residue_indices = {
            "solvent": [np.arange(0, 5), np.arange(5, 10)],
        }

        result = action.run(
            trajectory=mock_trajectory,
            spacing=1.0,
            output_dir=tmp_output_dir,
            only_com=True,
            com_residue_indices=com_residue_indices,
        )

        assert result.success is True
        assert result.metadata["n_probes"] == 1  # Only solvent COM
        assert "solvent" in result.metadata["probe_names"][0]

    def test_density_only_com_without_indices_fails(self, mock_trajectory, tmp_output_dir):
        """Test that only_com without com_residue_indices fails validation."""
        action = DensityAction()

        errors = action.validate(
            mock_trajectory,
            only_com=True,
        )

        assert len(errors) > 0
        assert any("com_residue_indices" in e for e in errors)


class TestDensitySubregion:
    """Test subregion functionality for focused density calculation."""

    def test_density_with_subregion(self, mock_trajectory, tmp_output_dir):
        """Test density calculation with subregion."""
        action = DensityAction()

        probe_indices = {"test": np.arange(10)}

        # Define a subregion in the center of the coordinate space
        subregion = ((20.0, 20.0, 20.0), (30.0, 30.0, 30.0))

        result = action.run(
            trajectory=mock_trajectory,
            probe_indices=probe_indices,
            spacing=1.0,
            output_dir=tmp_output_dir,
            subregion=subregion,
        )

        assert result.success is True
        assert result.metadata["subregion"] == subregion

    def test_density_subregion_validation(self, mock_trajectory):
        """Test subregion validation."""
        action = DensityAction()

        # Invalid subregion (wrong format)
        errors = action.validate(
            mock_trajectory,
            probe_indices={"test": np.arange(10)},
            subregion=((1, 2, 3),),  # Only one point
        )
        assert len(errors) > 0

        # Valid subregion
        errors = action.validate(
            mock_trajectory,
            probe_indices={"test": np.arange(10)},
            subregion=((0, 0, 0), (10, 10, 10)),
        )
        assert len(errors) == 0


class TestDensityParallel:
    """Test parallel processing for density calculation."""

    def test_density_parallel_processing(self, mock_trajectory, tmp_output_dir):
        """Test density calculation with multiple workers."""
        action = DensityAction()

        probe_indices = {"test": np.arange(10)}

        result = action.run(
            trajectory=mock_trajectory,
            probe_indices=probe_indices,
            spacing=1.0,
            output_dir=tmp_output_dir,
            n_workers=2,
        )

        assert result.success is True
        assert result.metadata["parallel_workers"] == 2

    def test_density_sequential_vs_parallel_same_result(self, mock_trajectory, tmp_output_dir):
        """Test that sequential and parallel give same results."""
        action = DensityAction()
        probe_indices = {"test": np.arange(20)}

        # Sequential
        result_seq = action.run(
            trajectory=mock_trajectory,
            probe_indices=probe_indices,
            spacing=1.0,
            output_dir=tmp_output_dir,
            output_prefix="seq_",
            n_workers=1,
        )

        # Parallel (need to reset trajectory)
        mock_trajectory._current_frame = 0  # Reset for re-iteration
        result_par = action.run(
            trajectory=mock_trajectory,
            probe_indices=probe_indices,
            spacing=1.0,
            output_dir=tmp_output_dir,
            output_prefix="par_",
            n_workers=2,
        )

        assert result_seq.success is True
        assert result_par.success is True
        # Both should have processed same number of frames
        assert result_seq.metadata["n_frames"] == result_par.metadata["n_frames"]


class TestDensityProteinAction:
    """Test protein/solute density calculation."""

    def test_density_protein_action_exists(self):
        """Test density_protein action is registered."""
        assert get_action("density_protein") is DensityProteinAction

    def test_density_protein_calculation(self, mock_trajectory, tmp_output_dir):
        """Test protein density calculation."""

        action = DensityProteinAction()

        # Use first 50 atoms as "solute"
        solute_mask = np.zeros(mock_trajectory.n_atoms, dtype=bool)
        solute_mask[:50] = True

        result = action.run(
            trajectory=mock_trajectory,
            solute_mask=solute_mask,
            spacing=1.0,
            output_dir=tmp_output_dir,
        )

        assert result.success is True
        assert len(result.output_files) == 1
        assert "protein" in result.output_files[0].name

    def test_density_protein_with_indices(self, mock_trajectory, tmp_output_dir):
        """Test protein density with solute_indices instead of mask."""

        action = DensityProteinAction()

        # Use first 50 atoms as "solute"
        solute_indices = np.arange(50)

        result = action.run(
            trajectory=mock_trajectory,
            solute_indices=solute_indices,
            spacing=1.0,
            output_dir=tmp_output_dir,
        )

        assert result.success is True

    def test_density_protein_requires_mask_or_indices(self, mock_trajectory, tmp_output_dir):
        """Test that protein density requires solute mask or indices."""

        action = DensityProteinAction()

        result = action.run(
            trajectory=mock_trajectory,
            spacing=1.0,
            output_dir=tmp_output_dir,
        )

        assert result.success is False
        assert "solute" in result.error.lower()


class TestDensityAllHAAction:
    """Test all heavy atoms density calculation."""

    def test_density_all_ha_action_exists(self):
        """Test density_all_ha action is registered."""
        assert get_action("density_all_ha") is DensityAllHAAction

    def test_density_all_ha_calculation(self, mock_trajectory, tmp_output_dir):
        """Test all heavy atoms density calculation."""

        action = DensityAllHAAction()

        # Define heavy atom info for mock residues
        heavy_atom_info = {
            "ETA": {
                "C1": [0, 5, 10],  # Indices for C1 atoms across residues
                "O1": [1, 6, 11],  # Indices for O1 atoms
            },
            "ACN": {
                "N1": [20, 25],
            },
        }

        result = action.run(
            trajectory=mock_trajectory,
            heavy_atom_info=heavy_atom_info,
            spacing=1.0,
            output_dir=tmp_output_dir,
        )

        assert result.success is True
        # Should have grids for: ETA_C1, ETA_O1, ETA_COM, ACN_N1, ACN_COM
        assert result.metadata["n_probes"] == 5

    def test_density_all_ha_excludes_water(self, mock_trajectory, tmp_output_dir):
        """Test that WAT is excluded by default."""

        action = DensityAllHAAction()

        heavy_atom_info = {
            "ETA": {"C1": [0, 5]},
            "WAT": {"O": [10, 15]},  # Should be excluded
            "HOH": {"O": [20, 25]},  # Should also be excluded
        }

        result = action.run(
            trajectory=mock_trajectory,
            heavy_atom_info=heavy_atom_info,
            spacing=1.0,
            output_dir=tmp_output_dir,
        )

        assert result.success is True
        # Only ETA should be processed: ETA_C1 + ETA_COM = 2
        assert result.metadata["n_probes"] == 2


class TestCpptrajDensityAction:
    """Test cpptraj-based density calculation."""

    def test_cpptraj_density_action_exists(self):
        """Test cpptraj_density action is registered."""
        assert get_action("cpptraj_density") is CpptrajDensityAction

    def test_cpptraj_density_without_cpptraj(self, tmp_output_dir):
        """Test cpptraj action fails gracefully without cpptraj."""
        import os
        import shutil

        action = CpptrajDensityAction()

        # Create mock topology and trajectory files
        topo_path = tmp_output_dir / "system.prmtop"
        traj_path = tmp_output_dir / "aligned.nc"
        topo_path.touch()
        traj_path.touch()

        # Clear any cpptraj environment variable
        old_ptraj = os.environ.pop("AMBER_PTRAJ", None)

        # Only run test if cpptraj is not installed
        if shutil.which("cpptraj") is None:
            result = action.run(
                topology=topo_path,
                trajectory_pattern=[str(traj_path)],
                probe_masks={"OH": ":ETA@O1"},
                grid_dimensions=(50, 50, 50),
                grid_origin=(0.0, 0.0, 0.0),
                grid_spacing=0.5,
                output_dir=tmp_output_dir,
            )

            # Should fail with informative message
            assert result.success is False
            assert "cpptraj" in result.error.lower()

        # Restore environment
        if old_ptraj:
            os.environ["AMBER_PTRAJ"] = old_ptraj


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
        assert len(result.output_files) == 3  # JSON + summary + data
        assert result.metadata["n_hotspots"] == 1
        assert result.metadata["n_frames"] == 10


class TestResidenceParallel:
    """Test parallel processing for residence calculation."""

    def test_residence_parallel_processing(self, mock_trajectory, tmp_output_dir):
        """Test residence calculation with multiple workers."""
        action = ResidenceAction()

        hotspot_coords = [(25.0, 25.0, 25.0)]
        n_atoms = mock_trajectory.n_atoms
        residue_ids = np.arange(n_atoms)
        residue_names = np.array(["ETA"] * n_atoms)

        result = action.run(
            trajectory=mock_trajectory,
            hotspot_coords=hotspot_coords,
            tolerance=15.0,
            residue_ids=residue_ids,
            residue_names=residue_names,
            output_dir=tmp_output_dir,
            n_workers=2,
        )

        assert result.success is True
        assert result.metadata["parallel_workers"] == 2

    def test_residence_with_track_residues(self, mock_trajectory, tmp_output_dir):
        """Test residence calculation with residue type filtering."""
        action = ResidenceAction()

        hotspot_coords = [(25.0, 25.0, 25.0)]
        n_atoms = mock_trajectory.n_atoms
        residue_ids = np.arange(n_atoms)
        # Mix of residue types
        residue_names = np.array(["ETA" if i % 2 == 0 else "WAT" for i in range(n_atoms)])

        result = action.run(
            trajectory=mock_trajectory,
            hotspot_coords=hotspot_coords,
            tolerance=15.0,
            residue_ids=residue_ids,
            residue_names=residue_names,
            track_residues=["ETA"],  # Only track ETA
            output_dir=tmp_output_dir,
        )

        assert result.success is True
        # Check that output mentions tracked residues
        import json

        with open(result.output_files[0]) as f:
            data = json.load(f)
        assert data["track_residues"] == ["ETA"]

    def test_residence_with_non_hydrogen_mask(self, mock_trajectory, tmp_output_dir):
        """Test residence calculation with non-hydrogen mask."""
        action = ResidenceAction()

        hotspot_coords = [(25.0, 25.0, 25.0)]
        n_atoms = mock_trajectory.n_atoms
        residue_ids = np.arange(n_atoms)

        # Mask out every other atom (simulating hydrogen exclusion)
        non_h_mask = np.array([i % 2 == 0 for i in range(n_atoms)])

        result = action.run(
            trajectory=mock_trajectory,
            hotspot_coords=hotspot_coords,
            tolerance=15.0,
            residue_ids=residue_ids,
            non_hydrogen_mask=non_h_mask,
            output_dir=tmp_output_dir,
        )

        assert result.success is True

    def test_residence_multiple_hotspots(self, mock_trajectory, tmp_output_dir):
        """Test residence calculation with multiple hotspots."""
        action = ResidenceAction()

        hotspot_coords = [
            (20.0, 20.0, 20.0),
            (30.0, 30.0, 30.0),
            (40.0, 40.0, 40.0),
        ]
        n_atoms = mock_trajectory.n_atoms
        residue_ids = np.arange(n_atoms)

        result = action.run(
            trajectory=mock_trajectory,
            hotspot_coords=hotspot_coords,
            tolerance=10.0,
            residue_ids=residue_ids,
            output_dir=tmp_output_dir,
        )

        assert result.success is True
        assert result.metadata["n_hotspots"] == 3
        assert len(result.metadata["hotspot_occupancies"]) == 3

    def test_residence_data_file_format(self, mock_trajectory, tmp_output_dir):
        """Test that residence data file has correct format."""
        action = ResidenceAction()

        hotspot_coords = [(25.0, 25.0, 25.0)]
        n_atoms = mock_trajectory.n_atoms
        residue_ids = np.arange(n_atoms)
        residue_names = np.array(["ETA"] * n_atoms)

        result = action.run(
            trajectory=mock_trajectory,
            hotspot_coords=hotspot_coords,
            tolerance=15.0,
            residue_ids=residue_ids,
            residue_names=residue_names,
            output_dir=tmp_output_dir,
        )

        assert result.success is True

        # Check data file format
        data_file = [f for f in result.output_files if "data" in f.name][0]
        with open(data_file) as f:
            content = f.read()

        # Should have header comments
        assert "# Residence results" in content
        # Should have RESNAME-RESID MAP
        assert "# RESNAME-RESID MAP" in content


class TestResidenceResult:
    """Test ResidenceResult dataclass."""

    def test_residence_result_occupancy(self):
        """Test occupancy calculation."""

        result = ResidenceResult(
            hotspot_id=0,
            hotspot_coord=(0.0, 0.0, 0.0),
            frame_residues={
                1: [1, 2],  # Occupied
                2: [0],  # Not occupied (0 marker)
                3: [3],  # Occupied
                4: [0],  # Not occupied
            },
            total_frames=4,
        )

        # 2 out of 4 frames occupied
        assert result.occupancy == 0.5

    def test_residence_result_top_residues(self):
        """Test top residues calculation."""

        result = ResidenceResult(
            hotspot_id=0,
            hotspot_coord=(0.0, 0.0, 0.0),
            residue_counts={1: 10, 2: 5, 3: 8, 4: 3},
            total_frames=10,
        )

        top = result.top_residues(3)
        assert len(top) == 3
        assert top[0] == (1, 10)  # Most frequent
        assert top[1] == (3, 8)
        assert top[2] == (2, 5)


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
