"""Tests for pymdmix.io module."""
import struct
import tempfile
from pathlib import Path

import numpy as np
import pytest

from pymdmix.core.grid import Grid
from pymdmix.io.grids import (
    read_dx,
    write_dx,
    read_mrc,
    write_mrc,
    convert_grid,
    GridFormat,
    grid_info,
)


# =============================================================================
# Fixtures
# =============================================================================

@pytest.fixture
def sample_grid():
    """Create a sample Grid for testing."""
    # Create grid using from_bounds
    min_coord = np.array([0.0, 0.0, 0.0])
    max_coord = np.array([10.0, 12.0, 8.0])
    grid = Grid.from_bounds(min_coord, max_coord, spacing=1.0)
    # Add some random data
    grid.data = np.random.randn(*grid.shape).astype(np.float64)
    return grid


@pytest.fixture
def sample_dx_file(sample_grid, tmp_path):
    """Create a sample DX file."""
    dx_path = tmp_path / "test.dx"
    sample_grid.write_dx(dx_path)
    return dx_path


# =============================================================================
# GridFormat Tests
# =============================================================================

class TestGridFormat:
    def test_from_path_dx(self):
        assert GridFormat.from_path("file.dx") == GridFormat.DX
        assert GridFormat.from_path("file.opendx") == GridFormat.DX
    
    def test_from_path_mrc(self):
        assert GridFormat.from_path("file.mrc") == GridFormat.MRC
        assert GridFormat.from_path("file.map") == GridFormat.MRC
        assert GridFormat.from_path("file.ccp4") == GridFormat.MRC
    
    def test_from_path_unknown(self):
        with pytest.raises(ValueError, match="Unknown grid format"):
            GridFormat.from_path("file.xyz")


# =============================================================================
# DX I/O Tests
# =============================================================================

class TestDXIO:
    def test_read_write_roundtrip(self, sample_grid, tmp_path):
        """Test that write -> read preserves data."""
        dx_path = tmp_path / "test.dx"
        
        write_dx(sample_grid, dx_path)
        loaded = read_dx(dx_path)
        
        assert loaded.shape == sample_grid.shape
        np.testing.assert_allclose(loaded.origin, sample_grid.origin)
        assert loaded.spacing == sample_grid.spacing
        np.testing.assert_allclose(loaded.data, sample_grid.data, rtol=1e-5)
    
    def test_read_dx(self, sample_dx_file, sample_grid):
        """Test reading existing DX file."""
        loaded = read_dx(sample_dx_file)
        assert loaded.shape == sample_grid.shape


# =============================================================================
# MRC I/O Tests
# =============================================================================

class TestMRCIO:
    def test_write_mrc_creates_file(self, sample_grid, tmp_path):
        """Test that write_mrc creates a file."""
        mrc_path = tmp_path / "test.mrc"
        write_mrc(sample_grid, mrc_path)
        assert mrc_path.exists()
    
    def test_write_mrc_header_size(self, sample_grid, tmp_path):
        """MRC header should be 1024 bytes."""
        mrc_path = tmp_path / "test.mrc"
        write_mrc(sample_grid, mrc_path)
        
        with open(mrc_path, "rb") as f:
            header = f.read(1024)
        
        # Check header exists
        assert len(header) == 1024
        
        # Check dimensions
        nx, ny, nz = struct.unpack('<iii', header[0:12])
        assert nx == sample_grid.shape[0]
        assert ny == sample_grid.shape[1]
        assert nz == sample_grid.shape[2]
        
        # Check mode (should be 2 for float32)
        mode = struct.unpack('<i', header[12:16])[0]
        assert mode == 2
        
        # Check MAP marker
        assert header[104:108] == b'MAP '
    
    def test_read_write_roundtrip(self, sample_grid, tmp_path):
        """Test that write -> read preserves data."""
        mrc_path = tmp_path / "test.mrc"
        
        write_mrc(sample_grid, mrc_path)
        loaded = read_mrc(mrc_path)
        
        assert loaded.shape == sample_grid.shape
        # MRC uses float32, so precision is lower
        np.testing.assert_allclose(loaded.data, sample_grid.data, rtol=1e-4)
    
    def test_read_mrc_with_offset_origin(self, tmp_path):
        """Test MRC with non-zero origin."""
        min_coord = np.array([10.0, 20.0, 30.0])
        max_coord = np.array([15.0, 25.0, 35.0])
        grid = Grid.from_bounds(min_coord, max_coord, spacing=1.0)
        grid.data = np.ones(grid.shape)
        
        mrc_path = tmp_path / "offset.mrc"
        write_mrc(grid, mrc_path)
        loaded = read_mrc(mrc_path)
        
        # Origin should be approximately preserved
        np.testing.assert_allclose(loaded.origin, grid.origin, atol=0.1)


# =============================================================================
# Format Conversion Tests
# =============================================================================

class TestConversion:
    def test_dx_to_mrc(self, sample_grid, tmp_path):
        """Test converting DX to MRC."""
        dx_path = tmp_path / "input.dx"
        mrc_path = tmp_path / "output.mrc"
        
        sample_grid.write_dx(dx_path)
        convert_grid(dx_path, mrc_path)
        
        assert mrc_path.exists()
        loaded = read_mrc(mrc_path)
        assert loaded.shape == sample_grid.shape
    
    def test_mrc_to_dx(self, sample_grid, tmp_path):
        """Test converting MRC to DX."""
        mrc_path = tmp_path / "input.mrc"
        dx_path = tmp_path / "output.dx"
        
        write_mrc(sample_grid, mrc_path)
        convert_grid(mrc_path, dx_path)
        
        assert dx_path.exists()
        loaded = read_dx(dx_path)
        assert loaded.shape == sample_grid.shape
    
    def test_explicit_format(self, sample_grid, tmp_path):
        """Test specifying output format explicitly."""
        dx_path = tmp_path / "input.dx"
        out_path = tmp_path / "output.dat"  # Unusual extension
        
        sample_grid.write_dx(dx_path)
        convert_grid(dx_path, out_path, output_format=GridFormat.MRC)
        
        # Should be readable as MRC
        loaded = read_mrc(out_path)
        assert loaded.shape == sample_grid.shape


# =============================================================================
# Grid Info Tests
# =============================================================================

class TestGridInfo:
    def test_dx_info(self, sample_dx_file, sample_grid):
        """Test getting info from DX file."""
        info = grid_info(sample_dx_file)
        
        assert info["format"] == "DX"
        assert info["shape"] == sample_grid.shape
    
    def test_mrc_info(self, sample_grid, tmp_path):
        """Test getting info from MRC file."""
        mrc_path = tmp_path / "test.mrc"
        write_mrc(sample_grid, mrc_path)
        
        info = grid_info(mrc_path)
        
        assert info["format"] == "MRC"
        assert info["shape"] == sample_grid.shape
        assert "min" in info
        assert "max" in info
        assert "mean" in info


# =============================================================================
# Edge Cases
# =============================================================================

class TestEdgeCases:
    def test_small_grid(self, tmp_path):
        """Test with very small grid."""
        min_coord = np.array([0.0, 0.0, 0.0])
        max_coord = np.array([2.0, 2.0, 2.0])
        grid = Grid.from_bounds(min_coord, max_coord, spacing=1.0)
        grid.data = np.array([[[1, 2], [3, 4]], [[5, 6], [7, 8]]], dtype=np.float64)
        
        mrc_path = tmp_path / "small.mrc"
        write_mrc(grid, mrc_path)
        loaded = read_mrc(mrc_path)
        
        np.testing.assert_allclose(loaded.data, grid.data, rtol=1e-4)
    
    def test_anisotropic_spacing_via_bounds(self, tmp_path):
        """Test grid with effective anisotropic cell dimensions."""
        # Create rectangular grid (different dimensions, same spacing)
        min_coord = np.array([0.0, 0.0, 0.0])
        max_coord = np.array([10.0, 20.0, 30.0])
        grid = Grid.from_bounds(min_coord, max_coord, spacing=1.0)
        grid.data = np.random.randn(*grid.shape)
        
        mrc_path = tmp_path / "rect.mrc"
        write_mrc(grid, mrc_path)
        loaded = read_mrc(mrc_path)
        
        # Shape should be preserved
        assert loaded.shape == grid.shape
    
    def test_negative_values(self, tmp_path):
        """Test grid with negative values (energy grids)."""
        min_coord = np.array([0.0, 0.0, 0.0])
        max_coord = np.array([5.0, 5.0, 5.0])
        grid = Grid.from_bounds(min_coord, max_coord, spacing=1.0)
        grid.data = np.random.randn(*grid.shape) - 2.0  # Centered around -2
        
        mrc_path = tmp_path / "negative.mrc"
        write_mrc(grid, mrc_path)
        loaded = read_mrc(mrc_path)
        
        np.testing.assert_allclose(loaded.data, grid.data, rtol=1e-4)
