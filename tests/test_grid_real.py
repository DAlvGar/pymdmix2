"""
Tests for Grid I/O and analysis using pre-computed reference DX files.

These tests load the legacy ``ETA_CT.dx`` density grid computed from the
real test trajectory and verify that the Grid class, DX round-trips, free
energy conversion, and hotspot detection produce correct results.

Marked ``@pytest.mark.real_data`` — automatically skipped when Git LFS data
is not checked out.  Run locally after::

    git lfs pull
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from pymdmix.core.grid import Grid

# ---------------------------------------------------------------------------
# Expected characteristics of ETA_CT.dx
# ---------------------------------------------------------------------------

# Grid was computed from a 20-frame WAT-solvated peptide trajectory with
# ETA (ethanol) probe CT (carbon tail) atom.
_EXPECTED_SHAPE = (132, 132, 132)
_EXPECTED_ORIGIN = (-33.0, -33.0, -33.0)
_EXPECTED_SPACING = 0.5  # Å
_EXPECTED_MAX_COUNT = 853.0  # maximum raw count in the grid
_EXPECTED_MEAN_APPROX = 4.2  # approximate mean (±1.0 tolerance)


# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------


def _is_real_dx(path: Path) -> bool:
    """Return False if the file is an LFS pointer rather than a real DX file."""
    return path.exists() and path.stat().st_size > 1_000_000


# ===========================================================================
# 1.  Loading the pre-computed ETA_CT.dx
# ===========================================================================


@pytest.mark.real_data
class TestEtaCTDxLoad:
    """Load tests/data/grids/ETA_CT.dx and verify its properties."""

    def test_dx_file_loads_without_error(self, eta_ct_dx_path):
        """Grid.read_dx() parses ETA_CT.dx without raising exceptions."""
        if not _is_real_dx(eta_ct_dx_path):
            pytest.skip("ETA_CT.dx is an LFS pointer — run: git lfs pull")
        grid = Grid.read_dx(eta_ct_dx_path)
        assert isinstance(grid, Grid)

    def test_grid_shape_matches_expected(self, eta_ct_dx_path):
        """Grid shape is (132, 132, 132)."""
        if not _is_real_dx(eta_ct_dx_path):
            pytest.skip("ETA_CT.dx is an LFS pointer — run: git lfs pull")
        grid = Grid.read_dx(eta_ct_dx_path)
        assert grid.data.shape == _EXPECTED_SHAPE

    def test_grid_origin_matches_expected(self, eta_ct_dx_path):
        """Grid origin is at (-33, -33, -33) Å."""
        if not _is_real_dx(eta_ct_dx_path):
            pytest.skip("ETA_CT.dx is an LFS pointer — run: git lfs pull")
        grid = Grid.read_dx(eta_ct_dx_path)
        np.testing.assert_allclose(grid.origin, _EXPECTED_ORIGIN, atol=0.01)

    def test_grid_spacing_matches_expected(self, eta_ct_dx_path):
        """Grid spacing is 0.5 Å."""
        if not _is_real_dx(eta_ct_dx_path):
            pytest.skip("ETA_CT.dx is an LFS pointer — run: git lfs pull")
        grid = Grid.read_dx(eta_ct_dx_path)
        sp = grid.spacing if isinstance(grid.spacing, float) else grid.spacing[0]
        assert abs(sp - _EXPECTED_SPACING) < 0.001

    def test_grid_max_value(self, eta_ct_dx_path):
        """Maximum count value matches the reference."""
        if not _is_real_dx(eta_ct_dx_path):
            pytest.skip("ETA_CT.dx is an LFS pointer — run: git lfs pull")
        grid = Grid.read_dx(eta_ct_dx_path)
        assert abs(float(grid.data.max()) - _EXPECTED_MAX_COUNT) < 1.0

    def test_grid_mean_value(self, eta_ct_dx_path):
        """Mean count value is close to the expected reference value."""
        if not _is_real_dx(eta_ct_dx_path):
            pytest.skip("ETA_CT.dx is an LFS pointer — run: git lfs pull")
        grid = Grid.read_dx(eta_ct_dx_path)
        assert abs(float(grid.data.mean()) - _EXPECTED_MEAN_APPROX) < 1.0

    def test_grid_data_non_negative(self, eta_ct_dx_path):
        """All grid values are ≥ 0 (count data)."""
        if not _is_real_dx(eta_ct_dx_path):
            pytest.skip("ETA_CT.dx is an LFS pointer — run: git lfs pull")
        grid = Grid.read_dx(eta_ct_dx_path)
        assert float(grid.data.min()) >= 0.0

    def test_grid_has_nonzero_cells(self, eta_ct_dx_path):
        """The majority of grid cells are non-zero (dense solvent grid)."""
        if not _is_real_dx(eta_ct_dx_path):
            pytest.skip("ETA_CT.dx is an LFS pointer — run: git lfs pull")
        grid = Grid.read_dx(eta_ct_dx_path)
        assert np.count_nonzero(grid.data) > 1_000_000, (
            "Expected > 1M non-zero cells in ETA_CT.dx"
        )

    def test_grid_metadata_present(self, eta_ct_dx_path):
        """Grid metadata dict is populated after loading."""
        if not _is_real_dx(eta_ct_dx_path):
            pytest.skip("ETA_CT.dx is an LFS pointer — run: git lfs pull")
        grid = Grid.read_dx(eta_ct_dx_path)
        assert isinstance(grid.metadata, dict)


# ===========================================================================
# 2.  DX round-trip (write and re-read)
# ===========================================================================


@pytest.mark.real_data
class TestDXRoundTrip:
    """Write the loaded grid back to DX and verify data is preserved."""

    def test_write_and_reload_preserves_shape(self, eta_ct_dx_path, tmp_output_dir):
        """Re-loaded grid has the same shape as the original."""
        if not _is_real_dx(eta_ct_dx_path):
            pytest.skip("ETA_CT.dx is an LFS pointer — run: git lfs pull")
        grid = Grid.read_dx(eta_ct_dx_path)
        out = tmp_output_dir / "roundtrip.dx"
        grid.write_dx(out)
        reloaded = Grid.read_dx(out)
        assert reloaded.data.shape == grid.data.shape

    def test_write_and_reload_preserves_values(self, eta_ct_dx_path, tmp_output_dir):
        """Re-loaded grid values match the original within float32 tolerance."""
        if not _is_real_dx(eta_ct_dx_path):
            pytest.skip("ETA_CT.dx is an LFS pointer — run: git lfs pull")
        grid = Grid.read_dx(eta_ct_dx_path)
        out = tmp_output_dir / "roundtrip_vals.dx"
        grid.write_dx(out)
        reloaded = Grid.read_dx(out)
        np.testing.assert_allclose(reloaded.data, grid.data, rtol=1e-5)

    def test_write_and_reload_preserves_origin(self, eta_ct_dx_path, tmp_output_dir):
        """Re-loaded grid has the same origin."""
        if not _is_real_dx(eta_ct_dx_path):
            pytest.skip("ETA_CT.dx is an LFS pointer — run: git lfs pull")
        grid = Grid.read_dx(eta_ct_dx_path)
        out = tmp_output_dir / "roundtrip_origin.dx"
        grid.write_dx(out)
        reloaded = Grid.read_dx(out)
        np.testing.assert_allclose(reloaded.origin, grid.origin, atol=0.01)

    def test_write_and_reload_preserves_spacing(self, eta_ct_dx_path, tmp_output_dir):
        """Re-loaded grid has the same spacing."""
        if not _is_real_dx(eta_ct_dx_path):
            pytest.skip("ETA_CT.dx is an LFS pointer — run: git lfs pull")
        grid = Grid.read_dx(eta_ct_dx_path)
        out = tmp_output_dir / "roundtrip_spacing.dx"
        grid.write_dx(out)
        reloaded = Grid.read_dx(out)
        orig_sp = grid.spacing if isinstance(grid.spacing, float) else grid.spacing[0]
        new_sp = reloaded.spacing if isinstance(reloaded.spacing, float) else reloaded.spacing[0]
        assert abs(orig_sp - new_sp) < 0.001


# ===========================================================================
# 3.  Free energy conversion on the real grid
# ===========================================================================


@pytest.mark.real_data
class TestFreeEnergyFromRealGrid:
    """Convert the ETA_CT.dx density to free energy and verify properties."""

    def test_free_energy_conversion_shape(self, eta_ct_dx_path):
        """FE grid has the same shape as the density grid."""
        if not _is_real_dx(eta_ct_dx_path):
            pytest.skip("ETA_CT.dx is an LFS pointer — run: git lfs pull")
        grid = Grid.read_dx(eta_ct_dx_path)
        fe = grid.to_free_energy()
        assert fe.data.shape == grid.data.shape

    def test_free_energy_has_negative_minimum(self, eta_ct_dx_path):
        """The most favourable FE cell is negative (attractive region exists)."""
        if not _is_real_dx(eta_ct_dx_path):
            pytest.skip("ETA_CT.dx is an LFS pointer — run: git lfs pull")
        grid = Grid.read_dx(eta_ct_dx_path)
        fe = grid.to_free_energy()
        populated = grid.data > 0
        assert float(fe.data[populated].min()) < 0.0, (
            "Expected some negative free energy values in populated cells"
        )

    def test_free_energy_finite_in_populated_cells(self, eta_ct_dx_path):
        """FE values are finite wherever the density is non-zero."""
        if not _is_real_dx(eta_ct_dx_path):
            pytest.skip("ETA_CT.dx is an LFS pointer — run: git lfs pull")
        grid = Grid.read_dx(eta_ct_dx_path)
        fe = grid.to_free_energy()
        populated = grid.data > 0
        assert np.all(np.isfinite(fe.data[populated])), (
            "Non-finite FE values in populated density cells"
        )

    def test_free_energy_below_threshold_exists(self, eta_ct_dx_path):
        """Some cells have FE ≤ −1.5 kcal/mol after Boltzmann inversion of density."""
        if not _is_real_dx(eta_ct_dx_path):
            pytest.skip("ETA_CT.dx is an LFS pointer — run: git lfs pull")
        grid = Grid.read_dx(eta_ct_dx_path)
        fe = grid.to_free_energy()
        n_hotspot = int(np.sum(fe.data <= -1.5))
        assert n_hotspot > 0, "No cells below −1.5 kcal/mol in ETA_CT.dx"


# ===========================================================================
# 4.  Hotspot detection on the real grid
# ===========================================================================


@pytest.mark.real_data
class TestHotspotsFromRealGrid:
    """Run HotspotAction on the ETA_CT.dx grid and check output."""

    def test_hotspot_action_produces_output(self, eta_ct_dx_path, tmp_output_dir):
        """HotspotAction writes at least one output file for the real density grid."""
        if not _is_real_dx(eta_ct_dx_path):
            pytest.skip("ETA_CT.dx is an LFS pointer — run: git lfs pull")
        from pymdmix.analysis.hotspots import HotspotAction

        # HotspotAction.run() converts density → FE internally; pass raw counts.
        density_grid = Grid.read_dx(eta_ct_dx_path)

        action = HotspotAction()
        result = action.run(
            grids={"ETA_CT": density_grid},
            energy_cutoff=-2.0,  # tight cutoff avoids memory-intensive clustering
            output_dir=tmp_output_dir,
        )
        assert result.success
        assert len(result.output_files) >= 1

    def test_hotspot_metadata_contains_energy(self, eta_ct_dx_path, tmp_output_dir):
        """Hotspot metadata includes the energy_cutoff key."""
        if not _is_real_dx(eta_ct_dx_path):
            pytest.skip("ETA_CT.dx is an LFS pointer — run: git lfs pull")
        from pymdmix.analysis.hotspots import HotspotAction

        density_grid = Grid.read_dx(eta_ct_dx_path)

        action = HotspotAction()
        result = action.run(
            grids={"ETA_CT": density_grid},
            energy_cutoff=-2.0,
            output_dir=tmp_output_dir,
        )
        assert "energy_cutoff" in result.metadata
