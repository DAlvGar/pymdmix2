"""Tests for HotSpotSet class."""

import numpy as np
import pytest

from pymdmix.analysis.hotspots import Hotspot, HotSpotSet


@pytest.fixture
def sample_hotspots():
    """Create sample hotspots for testing."""
    return [
        Hotspot(
            id=0,
            probe="OH",
            centroid=(0.0, 0.0, 0.0),
            energy=-2.0,
            volume=10.0,
            n_points=20,
            coords=np.random.randn(20, 3),
            energies=np.random.randn(20) - 2.0,
        ),
        Hotspot(
            id=1,
            probe="OH",
            centroid=(1.0, 1.0, 1.0),
            energy=-1.5,
            volume=8.0,
            n_points=15,
            coords=np.random.randn(15, 3) + 1,
            energies=np.random.randn(15) - 1.5,
        ),
        Hotspot(
            id=2,
            probe="OH",
            centroid=(10.0, 10.0, 10.0),
            energy=-0.5,
            volume=5.0,
            n_points=10,
            coords=np.random.randn(10, 3) + 10,
            energies=np.random.randn(10) - 0.5,
        ),
        Hotspot(
            id=3,
            probe="OH",
            centroid=(11.0, 10.0, 10.0),
            energy=-0.3,
            volume=3.0,
            n_points=5,
            coords=np.random.randn(5, 3) + np.array([11, 10, 10]),
            energies=np.random.randn(5) - 0.3,
        ),
    ]


class TestHotSpotSet:
    """Tests for HotSpotSet class."""

    def test_create_empty(self):
        """Test creating empty set."""
        hs_set = HotSpotSet(probe="OH", name="test")
        assert len(hs_set) == 0
        assert hs_set.n_hotspots == 0

    def test_create_with_hotspots(self, sample_hotspots):
        """Test creating with hotspots."""
        hs_set = HotSpotSet(probe="OH", name="test", hotspots=sample_hotspots)
        assert len(hs_set) == 4

    def test_add_hotspots_single(self, sample_hotspots):
        """Test adding single hotspot."""
        hs_set = HotSpotSet(probe="OH", name="test")
        hs_set.add_hotspots(sample_hotspots[0])
        assert len(hs_set) == 1

    def test_add_hotspots_list(self, sample_hotspots):
        """Test adding list of hotspots."""
        hs_set = HotSpotSet(probe="OH", name="test")
        hs_set.add_hotspots(sample_hotspots)
        assert len(hs_set) == 4

    def test_sorted_by_energy(self, sample_hotspots):
        """Test that hotspots are sorted by energy."""
        hs_set = HotSpotSet(probe="OH", name="test", hotspots=sample_hotspots)
        energies = [h.energy for h in hs_set]
        assert energies == sorted(energies)

    def test_get_by_id(self, sample_hotspots):
        """Test getting hotspots by ID."""
        hs_set = HotSpotSet(probe="OH", name="test", hotspots=sample_hotspots)
        result = hs_set.get_by_id(1)
        assert len(result) == 1
        assert result[0].id == 1

    def test_get_by_id_multiple(self, sample_hotspots):
        """Test getting multiple hotspots by ID."""
        hs_set = HotSpotSet(probe="OH", name="test", hotspots=sample_hotspots)
        result = hs_set.get_by_id([0, 2])
        assert len(result) == 2

    def test_centroids(self, sample_hotspots):
        """Test centroids property."""
        hs_set = HotSpotSet(probe="OH", name="test", hotspots=sample_hotspots)
        centroids = hs_set.centroids
        assert centroids.shape == (4, 3)

    def test_energies(self, sample_hotspots):
        """Test energies property."""
        hs_set = HotSpotSet(probe="OH", name="test", hotspots=sample_hotspots)
        energies = hs_set.energies
        assert len(energies) == 4

    def test_volumes(self, sample_hotspots):
        """Test volumes property."""
        hs_set = HotSpotSet(probe="OH", name="test", hotspots=sample_hotspots)
        volumes = hs_set.volumes
        assert len(volumes) == 4


class TestHotSpotSetFiltering:
    """Tests for HotSpotSet filtering methods."""

    def test_prune_by_energy(self, sample_hotspots):
        """Test pruning by energy."""
        hs_set = HotSpotSet(probe="OH", name="test", hotspots=sample_hotspots)
        pruned = hs_set.prune_by_energy(-1.0)
        assert len(pruned) == 2  # Only -2.0 and -1.5

    def test_prune_by_volume(self, sample_hotspots):
        """Test pruning by volume."""
        hs_set = HotSpotSet(probe="OH", name="test", hotspots=sample_hotspots)
        pruned = hs_set.prune_by_volume(6.0)
        assert len(pruned) == 2  # Only 10.0 and 8.0

    def test_prune_by_n_points(self, sample_hotspots):
        """Test pruning by number of points."""
        hs_set = HotSpotSet(probe="OH", name="test", hotspots=sample_hotspots)
        pruned = hs_set.prune_by_n_points(12)
        assert len(pruned) == 2  # Only 20 and 15


class TestHotSpotSetClustering:
    """Tests for HotSpotSet clustering."""

    def test_compute_distance_matrix(self, sample_hotspots):
        """Test distance matrix computation."""
        hs_set = HotSpotSet(probe="OH", name="test", hotspots=sample_hotspots)
        dm = hs_set.compute_distance_matrix()
        assert dm.shape == (4, 4)
        assert np.allclose(dm, dm.T)  # Symmetric
        assert np.allclose(np.diag(dm), 0)  # Zero diagonal

    def test_cluster(self, sample_hotspots):
        """Test clustering."""
        hs_set = HotSpotSet(probe="OH", name="test", hotspots=sample_hotspots)
        labels = hs_set.cluster(cutoff=5.0)
        assert len(labels) == 4
        # First two should cluster together, last two should cluster together
        # (distance ~1.7 vs ~1.0)

    def test_n_clusters(self, sample_hotspots):
        """Test cluster count."""
        hs_set = HotSpotSet(probe="OH", name="test", hotspots=sample_hotspots)
        hs_set.cluster(cutoff=3.0)
        # Two clusters: (0,1) near origin and (2,3) near (10,10,10)
        assert hs_set.n_clusters == 2

    def test_get_cluster_representatives(self, sample_hotspots):
        """Test getting cluster representatives."""
        hs_set = HotSpotSet(probe="OH", name="test", hotspots=sample_hotspots)
        reps = hs_set.get_cluster_representatives(cutoff=3.0)
        # Should have 2 representatives (one per cluster)
        assert len(reps) == 2


class TestHotSpotSetNearest:
    """Tests for finding nearest hotspot."""

    def test_find_nearest(self, sample_hotspots):
        """Test finding nearest hotspot."""
        hs_set = HotSpotSet(probe="OH", name="test", hotspots=sample_hotspots)
        nearest = hs_set.find_nearest((0.5, 0.5, 0.5))
        assert nearest is not None
        # Should be hotspot at (0,0,0) or (1,1,1)
        assert nearest.id in [0, 1]

    def test_find_nearest_beyond_max(self, sample_hotspots):
        """Test finding nearest with max distance."""
        hs_set = HotSpotSet(probe="OH", name="test", hotspots=sample_hotspots)
        nearest = hs_set.find_nearest((100.0, 100.0, 100.0), max_distance=5.0)
        assert nearest is None


class TestHotSpotSetIO:
    """Tests for HotSpotSet I/O."""

    def test_to_pdb(self, sample_hotspots, tmp_path):
        """Test PDB output."""
        hs_set = HotSpotSet(probe="OH", name="test", hotspots=sample_hotspots)
        pdb_path = tmp_path / "hotspots.pdb"
        hs_set.to_pdb(pdb_path)

        assert pdb_path.exists()
        content = pdb_path.read_text()
        assert "HETATM" in content
        assert "OH" in content

    def test_to_json(self, sample_hotspots, tmp_path):
        """Test JSON output."""
        hs_set = HotSpotSet(probe="OH", name="test", hotspots=sample_hotspots)
        json_path = tmp_path / "hotspots.json"
        hs_set.to_json(json_path)

        assert json_path.exists()
        import json

        with open(json_path) as f:
            data = json.load(f)
        assert data["n_hotspots"] == 4
        assert data["probe"] == "OH"

    def test_summary(self, sample_hotspots):
        """Test summary output."""
        hs_set = HotSpotSet(probe="OH", name="test", hotspots=sample_hotspots)
        summary = hs_set.summary()
        assert "HotSpotSet" in summary
        assert "OH" in summary
        assert "Energy range" in summary


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
