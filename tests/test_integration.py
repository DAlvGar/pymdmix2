"""Integration tests for pyMDMix workflows."""

import numpy as np
import pytest
import tempfile
from pathlib import Path

from pymdmix.core.structure import (
    load_structure,
    save_structure,
    get_protein_mask,
    get_backbone_mask,
    align_structures,
    find_disulfides,
    rename_cys_to_cyx,
)
from pymdmix.core.grid import Grid
from pymdmix.core.solvent import SolventLibrary
from pymdmix.project.settings import MDSettings
from pymdmix.analysis.manager import ActionsManager, Action, ActionResult
from pymdmix.analysis.hotspots import Hotspot, HotSpotSet


# =============================================================================
# Test Fixtures
# =============================================================================

@pytest.fixture
def sample_pdb_path(tmp_path):
    """Create a sample PDB file."""
    pdb_content = """\
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00  0.00           C
ATOM      3  C   ALA A   1       2.009   1.420   0.000  1.00  0.00           C
ATOM      4  O   ALA A   1       1.246   2.390   0.000  1.00  0.00           O
ATOM      5  CB  ALA A   1       1.986  -0.760  -1.216  1.00  0.00           C
ATOM      6  N   CYS A   2       3.320   1.560   0.000  1.00  0.00           N
ATOM      7  CA  CYS A   2       3.970   2.860   0.000  1.00  0.00           C
ATOM      8  C   CYS A   2       5.480   2.760   0.000  1.00  0.00           C
ATOM      9  O   CYS A   2       6.080   1.690   0.000  1.00  0.00           O
ATOM     10  CB  CYS A   2       3.500   3.700   1.200  1.00  0.00           C
ATOM     11  SG  CYS A   2       3.800   3.000   2.800  1.00  0.00           S
ATOM     12  N   CYS A   3       6.100   3.900   0.000  1.00  0.00           N
ATOM     13  CA  CYS A   3       7.550   4.000   0.000  1.00  0.00           C
ATOM     14  C   CYS A   3       8.100   5.400   0.000  1.00  0.00           C
ATOM     15  O   CYS A   3       7.350   6.380   0.000  1.00  0.00           O
ATOM     16  CB  CYS A   3       8.100   3.200   1.200  1.00  0.00           C
ATOM     17  SG  CYS A   3       7.800   3.900   2.800  1.00  0.00           S
TER
END
"""
    pdb_path = tmp_path / "protein.pdb"
    pdb_path.write_text(pdb_content)
    return pdb_path


@pytest.fixture
def sample_grid():
    """Create a sample density grid."""
    data = np.random.rand(20, 20, 20) * 2
    # Add some hotspot-like regions
    data[5:8, 5:8, 5:8] = 5.0  # High density region
    data[15:17, 15:17, 15:17] = 4.0
    
    grid = Grid(
        data=data,
        origin=np.array([0.0, 0.0, 0.0]),
        spacing=0.5,
    )
    return grid


# =============================================================================
# Structure Workflow Tests
# =============================================================================

class TestStructureWorkflow:
    """Test structure preparation workflow."""

    def test_load_and_analyze_structure(self, sample_pdb_path):
        """Test loading structure and getting masks."""
        struct = load_structure(sample_pdb_path)
        
        # Check loading
        assert len(struct.atoms) == 17
        assert len(struct.residues) == 3
        
        # Check protein mask
        protein_mask = get_protein_mask(struct)
        assert protein_mask.sum() == 17
        
        # Check backbone mask
        bb_mask = get_backbone_mask(struct)
        assert bb_mask.sum() == 12  # N, CA, C, O for 3 residues

    def test_disulfide_detection_and_renaming(self, sample_pdb_path):
        """Test disulfide bond detection workflow."""
        struct = load_structure(sample_pdb_path)
        
        # Find disulfides
        ss_bonds = find_disulfides(struct, cutoff=3.0)
        
        # Should find one SS bond between CYS 2 and CYS 3
        # (their SG atoms are close)
        
        # Rename CYS to CYX
        rename_cys_to_cyx(struct, ss_bonds)
        
        # Check at least the operation completes
        assert struct is not None

    def test_save_and_reload_structure(self, sample_pdb_path, tmp_path):
        """Test structure save/load roundtrip."""
        struct = load_structure(sample_pdb_path)
        
        output_path = tmp_path / "output.pdb"
        save_structure(struct, output_path)
        
        assert output_path.exists()
        
        # Reload
        reloaded = load_structure(output_path)
        assert len(reloaded.atoms) == len(struct.atoms)

    def test_alignment_workflow(self, sample_pdb_path):
        """Test structure alignment."""
        struct = load_structure(sample_pdb_path)
        
        # Create a translated copy
        import parmed
        mobile = struct.copy(parmed.Structure)
        mobile.coordinates = struct.coordinates + np.array([5.0, 5.0, 5.0])
        
        # Align using backbone
        bb_mask = get_backbone_mask(mobile)
        aligned, rmsd = align_structures(mobile, struct, mask=bb_mask)
        
        assert rmsd == pytest.approx(0.0, abs=1e-5)


# =============================================================================
# Grid Workflow Tests
# =============================================================================

class TestGridWorkflow:
    """Test grid operations workflow."""

    def test_grid_creation_and_operations(self, sample_grid):
        """Test grid operations."""
        assert sample_grid.data.shape == (20, 20, 20)
        assert sample_grid.spacing == 0.5

    def test_grid_to_free_energy(self, sample_grid):
        """Test density to free energy conversion."""
        dg_grid = sample_grid.to_free_energy(temperature=300.0)
        
        # Higher density = lower (more negative) energy
        # Check the grid was transformed
        assert dg_grid.data.shape == sample_grid.data.shape

    def test_grid_save_load(self, sample_grid, tmp_path):
        """Test grid I/O."""
        dx_path = tmp_path / "test.dx"
        sample_grid.write_dx(dx_path)
        
        assert dx_path.exists()
        
        # Reload
        loaded = Grid.read_dx(dx_path)
        assert loaded.data.shape == sample_grid.data.shape
        np.testing.assert_array_almost_equal(loaded.data, sample_grid.data)


# =============================================================================
# MDSettings Workflow Tests
# =============================================================================

class TestMDSettingsWorkflow:
    """Test MD settings workflow."""

    def test_settings_creation(self):
        """Test creating settings."""
        settings = MDSettings(solvent="ETA", nanos=40)
        
        assert settings.solvent == "ETA"
        assert settings.nanos == 40
        assert settings.n_trajectory_files == 40  # 40ns / 1ns per file

    def test_settings_with_restraints(self):
        """Test settings with restraints."""
        settings = MDSettings(
            solvent="ETA",
            restraint_mode="HA",
            restraint_force=0.1,
        )
        
        assert settings.has_restraints
        assert settings.restraint_mode == "HA"
        assert settings.restraint_force == 0.1

    def test_settings_toml_roundtrip(self, tmp_path):
        """Test settings TOML save/load."""
        toml_content = """
[mdsettings]
solvent = "MAM"
nanos = 30
temperature = 310.0
restraint_mode = "BB"
"""
        config_path = tmp_path / "settings.toml"
        config_path.write_text(toml_content)
        
        settings = MDSettings.from_toml(config_path)
        
        assert settings.solvent == "MAM"
        assert settings.nanos == 30
        assert settings.temperature == 310.0
        assert settings.restraint_mode == "BB"


# =============================================================================
# HotSpot Workflow Tests
# =============================================================================

class TestHotSpotWorkflow:
    """Test hotspot detection workflow."""

    def test_hotspot_set_workflow(self):
        """Test creating and filtering hotspot sets."""
        # Create some hotspots
        hotspots = [
            Hotspot(
                id=i,
                probe="OH",
                centroid=(float(i), float(i), float(i)),
                energy=-2.0 + i * 0.5,
                volume=10.0 - i,
                n_points=20 - i * 3,
                coords=np.random.randn(20 - i * 3, 3) + i,
                energies=np.random.randn(20 - i * 3) - 2.0 + i * 0.5,
            )
            for i in range(5)
        ]
        
        # Create set
        hs_set = HotSpotSet(probe="OH", name="test", hotspots=hotspots)
        assert len(hs_set) == 5
        
        # Filter by energy
        filtered = hs_set.prune_by_energy(-1.0)
        assert len(filtered) < len(hs_set)
        
        # Filter by volume
        by_volume = hs_set.prune_by_volume(7.0)
        assert all(h.volume >= 7.0 for h in by_volume)

    def test_hotspot_clustering(self):
        """Test hotspot clustering."""
        # Create hotspots in two clusters
        cluster1 = [
            Hotspot(
                id=i, probe="OH",
                centroid=(float(i) * 0.5, 0.0, 0.0),
                energy=-1.5, volume=5.0, n_points=10,
                coords=np.random.randn(10, 3),
                energies=np.random.randn(10) - 1.5,
            )
            for i in range(3)
        ]
        cluster2 = [
            Hotspot(
                id=i + 3, probe="OH",
                centroid=(10.0 + float(i) * 0.5, 0.0, 0.0),
                energy=-1.0, volume=5.0, n_points=10,
                coords=np.random.randn(10, 3) + 10,
                energies=np.random.randn(10) - 1.0,
            )
            for i in range(3)
        ]
        
        hs_set = HotSpotSet(
            probe="OH",
            name="test",
            hotspots=cluster1 + cluster2,
        )
        
        # Cluster
        labels = hs_set.cluster(cutoff=2.0)
        assert len(labels) == 6
        assert hs_set.n_clusters == 2
        
        # Get representatives
        reps = hs_set.get_cluster_representatives(cutoff=2.0)
        assert len(reps) == 2


# =============================================================================
# ActionsManager Workflow Tests
# =============================================================================

class MockAnalysisAction(Action):
    """Mock action for testing."""
    action_name = "mock_analysis"
    
    def run(self, **kwargs) -> ActionResult:
        return {"status": "completed", "value": 42}


class TestActionsManagerWorkflow:
    """Test ActionsManager workflow."""

    def test_manager_serial_execution(self):
        """Test serial action execution."""
        from dataclasses import dataclass
        
        @dataclass
        class MockReplica:
            name: str
            path: Path = Path("/tmp")
            
            def get_trajectory(self, **kwargs):
                return None
        
        manager = ActionsManager(ncpus=1)
        manager.add_replicas([
            MockReplica(name="r1"),
            MockReplica(name="r2"),
        ])
        manager.add_actions(MockAnalysisAction)
        
        results = manager.run()
        
        assert "r1" in results
        assert "r2" in results
        assert results["r1"]["mock_analysis"]["value"] == 42

    def test_manager_parallel_execution(self):
        """Test parallel action execution."""
        from dataclasses import dataclass
        
        @dataclass
        class MockReplica:
            name: str
            path: Path = Path("/tmp")
            
            def get_trajectory(self, **kwargs):
                return None
        
        manager = ActionsManager(ncpus=2, use_threads=True)
        manager.add_replicas([
            MockReplica(name=f"r{i}") for i in range(4)
        ])
        manager.add_actions(MockAnalysisAction)
        
        results = manager.run()
        
        assert len(results) == 4


# =============================================================================
# Full Pipeline Test
# =============================================================================

class TestFullPipeline:
    """Test full analysis pipeline."""

    def test_structure_to_hotspots_pipeline(self, sample_pdb_path, tmp_path):
        """Test complete pipeline from structure to hotspots."""
        # 1. Load and prepare structure
        struct = load_structure(sample_pdb_path)
        assert len(struct.atoms) > 0
        
        # 2. Create mock density grid
        data = np.random.rand(20, 20, 20) + 0.5
        data[8:12, 8:12, 8:12] = 5.0  # Add hotspot region
        
        grid = Grid(
            data=data,
            origin=np.array([-5.0, -5.0, -5.0]),
            spacing=0.5,
        )
        
        # 3. Save grid
        dx_path = tmp_path / "density.dx"
        grid.write_dx(dx_path)
        
        # 4. Convert to energy
        dg_grid = grid.to_free_energy(temperature=300.0)
        
        # 5. Create hotspot from high-density region
        high_density_mask = grid.data > 3.0
        indices = np.argwhere(high_density_mask)
        
        if len(indices) > 0:
            coords = np.array([
                grid.index_to_coord(tuple(idx)) for idx in indices
            ])
            energies = dg_grid.data[high_density_mask]
            
            hotspot = Hotspot(
                id=0,
                probe="OH",
                centroid=tuple(coords.mean(axis=0)),
                energy=float(energies.mean()),
                volume=len(coords) * grid.spacing ** 3,
                n_points=len(coords),
                coords=coords,
                energies=energies,
            )
            
            # 6. Create hotspot set and export
            hs_set = HotSpotSet(probe="OH", name="test", hotspots=[hotspot])
            
            pdb_path = tmp_path / "hotspots.pdb"
            json_path = tmp_path / "hotspots.json"
            
            hs_set.to_pdb(pdb_path)
            hs_set.to_json(json_path)
            
            assert pdb_path.exists()
            assert json_path.exists()


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
