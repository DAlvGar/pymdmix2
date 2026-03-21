"""
Tests for pymdmix.core.trajectory module.
"""

import pytest
import numpy as np

from pymdmix.core.trajectory import (
    Frame,
    open_trajectory,
    get_available_backends,
    FrameSliceReader,
)


class TestFrame:
    """Test Frame dataclass."""
    
    def test_frame_creation(self, sample_coordinates):
        """Test creating a Frame."""
        frame = Frame(coordinates=sample_coordinates)
        
        assert frame.coordinates.shape == sample_coordinates.shape
        assert frame.time is None
        assert frame.box is None
    
    def test_frame_with_metadata(self, sample_coordinates):
        """Test Frame with time and box."""
        box = np.array([50.0, 50.0, 50.0, 90.0, 90.0, 90.0])
        frame = Frame(
            coordinates=sample_coordinates,
            time=100.0,
            box=box,
        )
        
        assert frame.time == 100.0
        assert np.array_equal(frame.box, box)
    
    def test_frame_n_atoms(self, sample_coordinates):
        """Test n_atoms property."""
        frame = Frame(coordinates=sample_coordinates)
        assert frame.n_atoms == len(sample_coordinates)


class TestMockTrajectory:
    """Test mock trajectory functionality."""
    
    def test_mock_trajectory_iteration(self, mock_trajectory):
        """Test iterating over mock trajectory."""
        frames = list(mock_trajectory)
        
        assert len(frames) == 10
        for frame in frames:
            assert hasattr(frame, 'coordinates')
            assert frame.coordinates.shape == (100, 3)
    
    def test_mock_trajectory_properties(self, mock_trajectory):
        """Test mock trajectory properties."""
        assert mock_trajectory.n_frames == 10
        assert mock_trajectory.n_atoms == 100
        assert len(mock_trajectory) == 10


class TestFrameSliceReader:
    """Test trajectory slicing."""
    
    def test_slice_start(self, mock_trajectory):
        """Test slicing from start."""
        sliced = FrameSliceReader(mock_trajectory, start=5)
        
        assert sliced.n_frames == 5
        frames = list(sliced)
        assert len(frames) == 5
    
    def test_slice_stop(self, mock_trajectory):
        """Test slicing with stop."""
        sliced = FrameSliceReader(mock_trajectory, stop=5)
        
        assert sliced.n_frames == 5
    
    def test_slice_step(self, mock_trajectory):
        """Test slicing with step."""
        sliced = FrameSliceReader(mock_trajectory, step=2)
        
        assert sliced.n_frames == 5
        frames = list(sliced)
        assert len(frames) == 5
    
    def test_slice_combined(self, mock_trajectory):
        """Test combined slicing."""
        sliced = FrameSliceReader(mock_trajectory, start=2, stop=8, step=2)
        
        assert sliced.n_frames == 3  # frames 2, 4, 6


class TestBackendAvailability:
    """Test backend availability detection."""
    
    def test_get_available_backends(self):
        """Test getting available backends."""
        backends = get_available_backends()
        
        assert isinstance(backends, list)
        # At minimum, should have amber if netCDF4 is installed
        # or mdanalysis if that's installed


class TestOpenTrajectory:
    """Test trajectory opening (requires real files or mocking)."""
    
    def test_open_nonexistent_file(self, tmp_output_dir):
        """Test opening non-existent file raises error."""
        # May raise FileNotFoundError or ImportError depending on backend availability
        with pytest.raises((FileNotFoundError, ImportError)):
            open_trajectory(
                tmp_output_dir / "nonexistent.prmtop",
                tmp_output_dir / "nonexistent.nc",
            )
    
    def test_unknown_backend(self, tmp_output_dir):
        """Test unknown backend raises error."""
        # Create dummy files
        (tmp_output_dir / "test.prmtop").write_text("")
        (tmp_output_dir / "test.nc").write_text("")
        
        with pytest.raises(ValueError, match="Unknown backend"):
            open_trajectory(
                tmp_output_dir / "test.prmtop",
                tmp_output_dir / "test.nc",
                backend="nonexistent",
            )
