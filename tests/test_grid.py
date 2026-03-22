"""
Tests for pymdmix.core.grid module.
"""

import numpy as np
import pytest

from pymdmix.core.grid import Grid


class TestGridCreation:
    """Test Grid creation methods."""

    def test_create_from_bounds(self):
        """Test grid creation from bounding box."""
        grid = Grid.from_bounds(
            min_coord=(0, 0, 0),
            max_coord=(10, 10, 10),
            spacing=1.0,
        )

        assert grid.shape == (11, 11, 11)
        assert grid.origin == (0, 0, 0)
        assert grid.spacing == 1.0
        assert grid.data.dtype == np.float64
        assert np.all(grid.data == 0)

    def test_create_from_structure(self, sample_coordinates):
        """Test grid creation from structure coordinates."""
        grid = Grid.from_structure(
            sample_coordinates,
            spacing=0.5,
            padding=5.0,
        )

        # Grid should encompass all coordinates plus padding
        min_coord = sample_coordinates.min(axis=0) - 5.0
        max_coord = sample_coordinates.max(axis=0) + 5.0
        expected_shape = tuple(
            int(np.ceil((max_coord[i] - min_coord[i]) / 0.5)) + 1 for i in range(3)
        )

        assert grid.shape == expected_shape
        assert grid.spacing == 0.5

    def test_invalid_spacing(self):
        """Test that invalid spacing raises error."""
        with pytest.raises(ValueError, match="spacing must be positive"):
            Grid(
                data=np.zeros((10, 10, 10)),
                origin=(0, 0, 0),
                spacing=-1.0,
            )

    def test_invalid_dimensions(self):
        """Test that non-3D data raises error."""
        with pytest.raises(ValueError, match="must be 3D"):
            Grid(
                data=np.zeros((10, 10)),
                origin=(0, 0, 0),
                spacing=1.0,
            )


class TestGridOperations:
    """Test Grid operations."""

    def test_add_count_single(self):
        """Test adding single count."""
        grid = Grid.from_bounds((0, 0, 0), (10, 10, 10), spacing=1.0)

        # Add count at center
        result = grid.add_count(np.array([5.0, 5.0, 5.0]))

        assert result is True
        assert grid.data[5, 5, 5] == 1

    def test_add_count_outside(self):
        """Test adding count outside grid bounds."""
        grid = Grid.from_bounds((0, 0, 0), (10, 10, 10), spacing=1.0)

        # Add count outside bounds
        result = grid.add_count(np.array([100.0, 100.0, 100.0]))

        assert result is False
        assert np.sum(grid.data) == 0

    def test_add_counts_bulk(self, sample_coordinates):
        """Test bulk coordinate counting."""
        grid = Grid.from_structure(sample_coordinates, spacing=1.0, padding=5.0)

        # Add all coordinates
        n_added = grid.add_counts_bulk(sample_coordinates)

        assert n_added == len(sample_coordinates)
        assert np.sum(grid.data) == len(sample_coordinates)

    def test_coord_to_index(self):
        """Test coordinate to index conversion."""
        grid = Grid.from_bounds((0, 0, 0), (10, 10, 10), spacing=1.0)

        idx = grid.coord_to_index((5.0, 5.0, 5.0))
        assert idx == (5, 5, 5)

        idx = grid.coord_to_index((5.4, 5.4, 5.4))
        assert idx == (5, 5, 5)  # Should floor

    def test_index_to_coord(self):
        """Test index to coordinate conversion."""
        grid = Grid.from_bounds((0, 0, 0), (10, 10, 10), spacing=1.0)

        coord = grid.index_to_coord((5, 5, 5))
        assert coord == (5.0, 5.0, 5.0)

    def test_to_density(self):
        """Test conversion to density."""
        grid = Grid.from_bounds((0, 0, 0), (10, 10, 10), spacing=1.0)
        grid.data[5, 5, 5] = 100

        density = grid.to_density(n_frames=10)

        assert density.data[5, 5, 5] == 10.0
        assert density.origin == grid.origin
        assert density.spacing == grid.spacing

    def test_to_density_invalid_frames(self):
        """Test density conversion with invalid frame count."""
        grid = Grid.from_bounds((0, 0, 0), (10, 10, 10), spacing=1.0)

        with pytest.raises(ValueError, match="n_frames must be positive"):
            grid.to_density(n_frames=0)

    def test_to_free_energy(self):
        """Test conversion to free energy."""
        grid = Grid.from_bounds((0, 0, 0), (10, 10, 10), spacing=1.0)
        grid.data[:] = 1.0  # Uniform density
        grid.data[5, 5, 5] = 10.0  # Hotspot

        dg = grid.to_free_energy(temperature=300.0)

        # Hotspot should have negative (favorable) free energy relative to mean
        assert dg.data[5, 5, 5] < 0

    def test_grid_addition(self):
        """Test adding two grids."""
        grid1 = Grid.from_bounds((0, 0, 0), (10, 10, 10), spacing=1.0)
        grid2 = Grid.from_bounds((0, 0, 0), (10, 10, 10), spacing=1.0)

        grid1.data[5, 5, 5] = 5
        grid2.data[5, 5, 5] = 3

        result = grid1 + grid2

        assert result.data[5, 5, 5] == 8

    def test_grid_multiplication(self):
        """Test multiplying grid by scalar."""
        grid = Grid.from_bounds((0, 0, 0), (10, 10, 10), spacing=1.0)
        grid.data[5, 5, 5] = 5

        result = grid * 2.0

        assert result.data[5, 5, 5] == 10


class TestGridIO:
    """Test Grid file I/O."""

    def test_write_dx(self, tmp_output_dir):
        """Test writing DX file."""
        grid = Grid.from_bounds((0, 0, 0), (5, 5, 5), spacing=1.0)
        grid.data[2, 2, 2] = 100

        dx_path = tmp_output_dir / "test.dx"
        grid.write_dx(dx_path)

        assert dx_path.exists()
        content = dx_path.read_text()
        assert "gridpositions counts 6 6 6" in content
        assert "origin 0.000000 0.000000 0.000000" in content

    def test_read_dx(self, sample_dx_file):
        """Test reading DX file."""
        grid = Grid.read_dx(sample_dx_file)

        assert grid.shape == (3, 3, 3)
        assert grid.origin == (0.0, 0.0, 0.0)
        assert grid.spacing == 1.0
        assert grid.data[0, 0, 0] == 1.0

    def test_roundtrip_dx(self, tmp_output_dir, sample_coordinates):
        """Test write then read DX file."""
        # Create grid with data
        grid = Grid.from_structure(sample_coordinates, spacing=1.0, padding=2.0)
        grid.add_counts_bulk(sample_coordinates)

        # Write
        dx_path = tmp_output_dir / "roundtrip.dx"
        grid.write_dx(dx_path)

        # Read back
        loaded = Grid.read_dx(dx_path)

        assert loaded.shape == grid.shape
        # Compare origins with tolerance (DX format has limited precision)
        assert np.allclose(loaded.origin, grid.origin, rtol=1e-5)
        assert loaded.spacing == grid.spacing
        assert np.allclose(loaded.data, grid.data)


class TestGridProperties:
    """Test Grid property methods."""

    def test_dimensions(self):
        """Test grid physical dimensions."""
        grid = Grid.from_bounds((0, 0, 0), (10, 10, 10), spacing=0.5)

        dims = grid.dimensions
        # Shape should be (21, 21, 21), so dims = (10.5, 10.5, 10.5)
        assert all(d > 10 for d in dims)

    def test_center(self):
        """Test grid center."""
        grid = Grid.from_bounds((0, 0, 0), (10, 10, 10), spacing=1.0)

        center = grid.center
        assert all(4 < c < 6 for c in center)

    def test_repr(self):
        """Test grid string representation."""
        grid = Grid.from_bounds((0, 0, 0), (10, 10, 10), spacing=1.0)
        grid.data[5, 5, 5] = 100

        repr_str = repr(grid)
        assert "Grid" in repr_str
        assert "shape=" in repr_str
        assert "11, 11, 11" in repr_str
