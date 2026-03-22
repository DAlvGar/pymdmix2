"""Tests for energy conversion and Boltzmann averaging."""

import numpy as np
import pytest

from pymdmix.analysis.energy import (
    _boltzmann_average_energy,
    boltzmann_average,
    calculate_expected_density,
    density_to_free_energy,
    normalize_grid,
)
from pymdmix.core.grid import Grid


class TestDensityToFreeEnergy:
    """Tests for density to free energy conversion."""

    def test_basic_conversion(self):
        """Test basic free energy conversion."""
        grid = Grid(
            data=np.ones((5, 5, 5)) * 2.0,  # 2x reference density
            origin=(0.0, 0.0, 0.0),
            spacing=1.0,
        )

        fe_grid = density_to_free_energy(grid, temperature=300.0, reference_density=1.0)

        # Higher density = negative free energy
        assert fe_grid.data.mean() < 0

    def test_low_density_positive_energy(self):
        """Test that low density gives positive free energy."""
        grid = Grid(
            data=np.ones((5, 5, 5)) * 0.5,  # 0.5x reference density
            origin=(0.0, 0.0, 0.0),
            spacing=1.0,
        )

        fe_grid = density_to_free_energy(grid, temperature=300.0, reference_density=1.0)

        # Lower density = positive free energy
        assert fe_grid.data.mean() > 0

    def test_auto_reference_density(self):
        """Test automatic reference density calculation."""
        data = np.ones((5, 5, 5))
        data[2, 2, 2] = 5.0  # High density spot

        grid = Grid(data=data, origin=(0.0, 0.0, 0.0), spacing=1.0)

        # Should use mean as reference
        fe_grid = density_to_free_energy(grid, temperature=300.0)

        # High density point should be more negative
        assert fe_grid.data[2, 2, 2] < fe_grid.data[0, 0, 0]

    def test_temperature_dependence(self):
        """Test that higher temperature gives larger absolute energy values."""
        grid = Grid(
            data=np.ones((5, 5, 5)) * 2.0,
            origin=(0.0, 0.0, 0.0),
            spacing=1.0,
        )

        fe_300 = density_to_free_energy(grid, temperature=300.0, reference_density=1.0)
        fe_600 = density_to_free_energy(grid, temperature=600.0, reference_density=1.0)

        # ΔG = -RT ln(ρ/ρ₀), higher T = larger absolute ΔG
        assert abs(fe_600.data.mean()) > abs(fe_300.data.mean())


class TestBoltzmannAverage:
    """Tests for Boltzmann averaging."""

    def test_single_grid(self):
        """Test that single grid returns unchanged."""
        grid = Grid(
            data=np.random.rand(5, 5, 5),
            origin=(0.0, 0.0, 0.0),
            spacing=1.0,
        )

        result = boltzmann_average([grid])

        np.testing.assert_array_equal(result.data, grid.data)

    def test_equal_grids_average(self):
        """Test averaging identical grids."""
        data = np.random.rand(5, 5, 5)
        grids = [
            Grid(data=data.copy(), origin=(0.0, 0.0, 0.0), spacing=1.0),
            Grid(data=data.copy(), origin=(0.0, 0.0, 0.0), spacing=1.0),
        ]

        result = boltzmann_average(grids)

        np.testing.assert_array_almost_equal(result.data, data)

    def test_weighted_average_density(self):
        """Test weighted average of density grids."""
        grid1 = Grid(
            data=np.ones((5, 5, 5)) * 1.0,
            origin=(0.0, 0.0, 0.0),
            spacing=1.0,
        )
        grid2 = Grid(
            data=np.ones((5, 5, 5)) * 3.0,
            origin=(0.0, 0.0, 0.0),
            spacing=1.0,
        )

        result = boltzmann_average([grid1, grid2], weights=[0.75, 0.25])

        # 0.75 * 1 + 0.25 * 3 = 1.5
        np.testing.assert_array_almost_equal(result.data, np.ones((5, 5, 5)) * 1.5)

    def test_from_file_paths(self, tmp_path):
        """Test averaging from file paths."""
        grid1 = Grid(
            data=np.ones((5, 5, 5)) * 1.0,
            origin=(0.0, 0.0, 0.0),
            spacing=1.0,
        )
        grid2 = Grid(
            data=np.ones((5, 5, 5)) * 2.0,
            origin=(0.0, 0.0, 0.0),
            spacing=1.0,
        )

        path1 = tmp_path / "grid1.dx"
        path2 = tmp_path / "grid2.dx"
        grid1.write_dx(path1)
        grid2.write_dx(path2)

        result = boltzmann_average([path1, path2])

        # Equal weights: (1 + 2) / 2 = 1.5
        np.testing.assert_array_almost_equal(result.data, np.ones((5, 5, 5)) * 1.5)

    def test_incompatible_shapes_raises(self):
        """Test that incompatible grid shapes raise error."""
        grid1 = Grid(data=np.ones((5, 5, 5)), origin=(0.0, 0.0, 0.0), spacing=1.0)
        grid2 = Grid(data=np.ones((6, 6, 6)), origin=(0.0, 0.0, 0.0), spacing=1.0)

        with pytest.raises(ValueError, match="shape"):
            boltzmann_average([grid1, grid2])


class TestExpectedDensity:
    """Tests for expected density calculation."""

    def test_expected_density(self):
        """Test expected density calculation."""
        result = calculate_expected_density(
            n_probe_atoms=100,
            box_volume=1000.0,  # 10^3 Å³
            n_frames=10,
            voxel_volume=1.0,  # 1 Å³
        )

        # 100 atoms / 1000 Å³ = 0.1 per voxel per frame
        # 0.1 * 10 frames = 1.0 expected total
        assert result == pytest.approx(1.0)

    def test_normalize_grid(self):
        """Test grid normalization."""
        grid = Grid(
            data=np.ones((5, 5, 5)) * 100.0,
            origin=(0.0, 0.0, 0.0),
            spacing=1.0,
        )

        normalized = normalize_grid(grid, expected_density=50.0)

        # 100 / 50 = 2.0
        np.testing.assert_array_almost_equal(normalized.data, np.ones((5, 5, 5)) * 2.0)


class TestBoltzmannAverageEnergy:
    """Tests for Boltzmann energy averaging."""

    def test_energy_average_stability(self):
        """Test numerical stability of energy averaging."""
        # Create energy grids with large negative values
        data1 = np.ones((5, 5, 5)) * -5.0
        data2 = np.ones((5, 5, 5)) * -3.0

        grids = [
            Grid(data=data1, origin=(0.0, 0.0, 0.0), spacing=1.0),
            Grid(data=data2, origin=(0.0, 0.0, 0.0), spacing=1.0),
        ]

        result = _boltzmann_average_energy(grids, [0.5, 0.5], temperature=300.0)

        # Result should be between -5 and -3, weighted toward lower energy
        assert -5.0 < result.mean() < -3.0
        # Lower energy should dominate
        assert result.mean() < -4.0
