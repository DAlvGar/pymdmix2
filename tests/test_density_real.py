"""
Tests for probe density calculation using the real Amber NetCDF test data.

These tests run :class:`~pymdmix.analysis.density.DensityAction` (pure-Python)
and, when cpptraj is available,
:class:`~pymdmix.analysis.density.CpptrajDensityAction` against the real
solvated-peptide trajectory.

Marked ``@pytest.mark.real_data`` — automatically skipped when Git LFS data
is not checked out.  Run locally after::

    git lfs pull
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from pymdmix.analysis.density import DensityAction
from pymdmix.core.grid import Grid
from pymdmix.core.structure import load_structure
from pymdmix.core.trajectory import open_trajectory

# ---------------------------------------------------------------------------
# Expected constants for the bundled system
# ---------------------------------------------------------------------------

_EXPECTED_FRAMES = 20
_EXPECTED_ATOMS = 7079
# WAT oxygen atom count derived at fixture setup
_WAT_RESNAME = "WAT"


# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------


def _is_real_trajectory(path: Path) -> bool:
    """Return False if the file is an LFS pointer rather than real binary data."""
    return path.exists() and path.stat().st_size > 100_000


def _get_wat_o_indices(topology_path: Path) -> np.ndarray:
    """Return atom indices for water oxygen atoms."""
    import MDAnalysis as mda

    u = mda.Universe(str(topology_path))
    return u.select_atoms(f"resname {_WAT_RESNAME} and name O").indices.astype(np.int64)


# ===========================================================================
# 1.  DensityAction (pure Python) with real trajectory
# ===========================================================================


@pytest.mark.real_data
class TestDensityActionRealTrajectory:
    """Run the pure-Python DensityAction on the real solvated-peptide trajectory."""

    def test_density_produces_dx_file(
        self, solvated_pep_pdb_path, amber_traj_nc_path, tmp_output_dir
    ):
        """DensityAction writes a .dx grid file to disk."""
        if not _is_real_trajectory(amber_traj_nc_path):
            pytest.skip("traj.nc is an LFS pointer — run: git lfs pull")
        traj = open_trajectory(solvated_pep_pdb_path, amber_traj_nc_path)
        ref = load_structure(solvated_pep_pdb_path)
        wat_o = _get_wat_o_indices(solvated_pep_pdb_path)

        action = DensityAction()
        result = action.run(
            trajectory=traj,
            reference=ref,
            probe_indices={"WAT_O": wat_o},
            spacing=1.0,
            output_dir=tmp_output_dir,
        )

        assert len(result.output_files) >= 1
        for f in result.output_files:
            assert Path(f).exists(), f"Expected output file {f} does not exist"

    def test_density_grid_shape_is_3d(
        self, solvated_pep_pdb_path, amber_traj_nc_path, tmp_output_dir
    ):
        """Density grid is a 3-D array."""
        if not _is_real_trajectory(amber_traj_nc_path):
            pytest.skip("traj.nc is an LFS pointer — run: git lfs pull")
        traj = open_trajectory(solvated_pep_pdb_path, amber_traj_nc_path)
        ref = load_structure(solvated_pep_pdb_path)
        wat_o = _get_wat_o_indices(solvated_pep_pdb_path)

        action = DensityAction()
        result = action.run(
            trajectory=traj,
            reference=ref,
            probe_indices={"WAT_O": wat_o},
            spacing=1.0,
            output_dir=tmp_output_dir,
        )
        grid = Grid.read_dx(result.output_files[0])
        assert grid.data.ndim == 3

    def test_density_values_non_negative(
        self, solvated_pep_pdb_path, amber_traj_nc_path, tmp_output_dir
    ):
        """All density values are non-negative (raw counts or normalised)."""
        if not _is_real_trajectory(amber_traj_nc_path):
            pytest.skip("traj.nc is an LFS pointer — run: git lfs pull")
        traj = open_trajectory(solvated_pep_pdb_path, amber_traj_nc_path)
        ref = load_structure(solvated_pep_pdb_path)
        wat_o = _get_wat_o_indices(solvated_pep_pdb_path)

        action = DensityAction()
        result = action.run(
            trajectory=traj,
            reference=ref,
            probe_indices={"WAT_O": wat_o},
            spacing=1.0,
            output_dir=tmp_output_dir,
        )
        grid = Grid.read_dx(result.output_files[0])
        assert float(grid.data.min()) >= 0.0, "Density values must be non-negative"

    def test_density_has_nonzero_values(
        self, solvated_pep_pdb_path, amber_traj_nc_path, tmp_output_dir
    ):
        """Density grid contains some non-zero cells (water was actually counted)."""
        if not _is_real_trajectory(amber_traj_nc_path):
            pytest.skip("traj.nc is an LFS pointer — run: git lfs pull")
        traj = open_trajectory(solvated_pep_pdb_path, amber_traj_nc_path)
        ref = load_structure(solvated_pep_pdb_path)
        wat_o = _get_wat_o_indices(solvated_pep_pdb_path)

        action = DensityAction()
        result = action.run(
            trajectory=traj,
            reference=ref,
            probe_indices={"WAT_O": wat_o},
            spacing=1.0,
            output_dir=tmp_output_dir,
        )
        grid = Grid.read_dx(result.output_files[0])
        assert np.count_nonzero(grid.data) > 0, "Density grid has no non-zero cells"

    def test_density_grid_spacing_matches_request(
        self, solvated_pep_pdb_path, amber_traj_nc_path, tmp_output_dir
    ):
        """Grid spacing in the output DX file matches the requested spacing."""
        if not _is_real_trajectory(amber_traj_nc_path):
            pytest.skip("traj.nc is an LFS pointer — run: git lfs pull")
        traj = open_trajectory(solvated_pep_pdb_path, amber_traj_nc_path)
        ref = load_structure(solvated_pep_pdb_path)
        wat_o = _get_wat_o_indices(solvated_pep_pdb_path)

        requested_spacing = 1.0
        action = DensityAction()
        result = action.run(
            trajectory=traj,
            reference=ref,
            probe_indices={"WAT_O": wat_o},
            spacing=requested_spacing,
            output_dir=tmp_output_dir,
        )
        grid = Grid.read_dx(result.output_files[0])
        sp = grid.spacing if isinstance(grid.spacing, float) else grid.spacing[0]
        assert abs(sp - requested_spacing) < 0.01

    def test_density_to_free_energy_has_finite_values(
        self, solvated_pep_pdb_path, amber_traj_nc_path, tmp_output_dir
    ):
        """Free energy conversion of real density grid produces finite values in populated cells."""
        if not _is_real_trajectory(amber_traj_nc_path):
            pytest.skip("traj.nc is an LFS pointer — run: git lfs pull")
        traj = open_trajectory(solvated_pep_pdb_path, amber_traj_nc_path)
        ref = load_structure(solvated_pep_pdb_path)
        wat_o = _get_wat_o_indices(solvated_pep_pdb_path)

        action = DensityAction()
        result = action.run(
            trajectory=traj,
            reference=ref,
            probe_indices={"WAT_O": wat_o},
            spacing=1.0,
            output_dir=tmp_output_dir,
        )
        density_grid = Grid.read_dx(result.output_files[0])
        fe_grid = density_grid.to_free_energy()

        # Populated cells must have finite FE
        populated = density_grid.data > 0
        assert np.all(np.isfinite(fe_grid.data[populated])), (
            "Non-finite FE values in populated density cells"
        )

    def test_density_result_metadata(
        self, solvated_pep_pdb_path, amber_traj_nc_path, tmp_output_dir
    ):
        """ActionResult metadata includes n_frames equal to the trajectory length."""
        if not _is_real_trajectory(amber_traj_nc_path):
            pytest.skip("traj.nc is an LFS pointer — run: git lfs pull")
        traj = open_trajectory(solvated_pep_pdb_path, amber_traj_nc_path)
        ref = load_structure(solvated_pep_pdb_path)
        wat_o = _get_wat_o_indices(solvated_pep_pdb_path)

        action = DensityAction()
        result = action.run(
            trajectory=traj,
            reference=ref,
            probe_indices={"WAT_O": wat_o},
            spacing=1.0,
            output_dir=tmp_output_dir,
        )
        assert result.success
        assert "n_frames" in result.metadata
        assert result.metadata["n_frames"] == _EXPECTED_FRAMES


# ===========================================================================
# 2.  MDAnalysis selection-based probes
# ===========================================================================


@pytest.mark.real_data
class TestDensityWithSelectionProbes:
    """Use probe_selections (MDAnalysis syntax) instead of pre-computed indices."""

    def test_density_with_selection_string(
        self, solvated_pep_pdb_path, amber_traj_nc_path, tmp_output_dir
    ):
        """probe_selections dictionary produces the same probe as probe_indices."""
        if not _is_real_trajectory(amber_traj_nc_path):
            pytest.skip("traj.nc is an LFS pointer — run: git lfs pull")
        traj = open_trajectory(solvated_pep_pdb_path, amber_traj_nc_path, backend="mdanalysis")
        ref = load_structure(solvated_pep_pdb_path)

        action = DensityAction()
        result = action.run(
            trajectory=traj,
            reference=ref,
            probe_selections={"WAT_O": "resname WAT and name O"},
            spacing=1.0,
            output_dir=tmp_output_dir,
        )
        assert len(result.output_files) >= 1
        grid = Grid.read_dx(result.output_files[0])
        assert np.count_nonzero(grid.data) > 0


# ===========================================================================
# 3.  CpptrajDensityAction with real trajectory (requires cpptraj)
# ===========================================================================


@pytest.mark.real_data
@pytest.mark.cpptraj
class TestCpptrajDensityRealTrajectory:
    """Run CpptrajDensityAction on the real solvated-peptide trajectory.

    The real trajectory ``tests/data/amber/traj.nc`` contains 7079 atoms
    (solvated peptide + 2319 WAT molecules).  cpptraj requires a topology
    that matches this atom count, so we use ``pep_WAT_WAT_1.pdb`` as the
    topology file (cpptraj accepts PDB as a topology source).

    Grid: 50 × 50 × 50 cells at 1 Å spacing centred near the peptide
    (peptide centre-of-geometry is close to the origin in the reference frame).
    """

    # Amber mask for WAT oxygen atoms
    _WAT_MASK = ":WAT@O"
    # Grid parameters matching the solvated box (~47 Å side)
    _GRID_DIM = (50, 50, 50)
    _GRID_ORIGIN = (-25.0, -25.0, -25.0)
    _GRID_SPACING = 1.0

    def test_cpptraj_density_produces_dx(
        self, solvated_pep_pdb_path, amber_traj_nc_path, tmp_output_dir
    ):
        """CpptrajDensityAction writes a DX grid file when cpptraj is present."""
        if not _is_real_trajectory(amber_traj_nc_path):
            pytest.skip("traj.nc is an LFS pointer — run: git lfs pull")
        from pymdmix.analysis.density import CpptrajDensityAction

        action = CpptrajDensityAction()
        result = action.run(
            topology=solvated_pep_pdb_path,
            trajectory_pattern=[str(amber_traj_nc_path)],
            probe_masks={"WAT_O": self._WAT_MASK},
            grid_dimensions=self._GRID_DIM,
            grid_origin=self._GRID_ORIGIN,
            grid_spacing=self._GRID_SPACING,
            output_dir=tmp_output_dir,
            output_prefix="real_",
        )
        assert result.success, f"CpptrajDensityAction failed: {result.error}"
        assert len(result.output_files) >= 1
        for f in result.output_files:
            assert Path(f).exists(), f"Expected output file does not exist: {f}"

    def test_cpptraj_density_grid_shape(
        self, solvated_pep_pdb_path, amber_traj_nc_path, tmp_output_dir
    ):
        """Grid shape matches the requested dimensions."""
        if not _is_real_trajectory(amber_traj_nc_path):
            pytest.skip("traj.nc is an LFS pointer — run: git lfs pull")
        from pymdmix.analysis.density import CpptrajDensityAction

        action = CpptrajDensityAction()
        result = action.run(
            topology=solvated_pep_pdb_path,
            trajectory_pattern=[str(amber_traj_nc_path)],
            probe_masks={"WAT_O": self._WAT_MASK},
            grid_dimensions=self._GRID_DIM,
            grid_origin=self._GRID_ORIGIN,
            grid_spacing=self._GRID_SPACING,
            output_dir=tmp_output_dir,
            output_prefix="shape_",
        )
        assert result.success, f"CpptrajDensityAction failed: {result.error}"
        grid = Grid.read_dx(result.output_files[0])
        assert grid.data.shape == self._GRID_DIM

    def test_cpptraj_density_grid_non_negative(
        self, solvated_pep_pdb_path, amber_traj_nc_path, tmp_output_dir
    ):
        """CpptrajDensityAction grid values are non-negative (raw counts)."""
        if not _is_real_trajectory(amber_traj_nc_path):
            pytest.skip("traj.nc is an LFS pointer — run: git lfs pull")
        from pymdmix.analysis.density import CpptrajDensityAction

        action = CpptrajDensityAction()
        result = action.run(
            topology=solvated_pep_pdb_path,
            trajectory_pattern=[str(amber_traj_nc_path)],
            probe_masks={"WAT_O": self._WAT_MASK},
            grid_dimensions=self._GRID_DIM,
            grid_origin=self._GRID_ORIGIN,
            grid_spacing=self._GRID_SPACING,
            output_dir=tmp_output_dir,
            output_prefix="nonneg_",
        )
        assert result.success, f"CpptrajDensityAction failed: {result.error}"
        grid = Grid.read_dx(result.output_files[0])
        assert float(grid.data.min()) >= 0.0

    def test_cpptraj_density_grid_has_nonzero_counts(
        self, solvated_pep_pdb_path, amber_traj_nc_path, tmp_output_dir
    ):
        """WAT oxygen atoms are counted in the grid (non-trivial density)."""
        if not _is_real_trajectory(amber_traj_nc_path):
            pytest.skip("traj.nc is an LFS pointer — run: git lfs pull")
        from pymdmix.analysis.density import CpptrajDensityAction

        action = CpptrajDensityAction()
        result = action.run(
            topology=solvated_pep_pdb_path,
            trajectory_pattern=[str(amber_traj_nc_path)],
            probe_masks={"WAT_O": self._WAT_MASK},
            grid_dimensions=self._GRID_DIM,
            grid_origin=self._GRID_ORIGIN,
            grid_spacing=self._GRID_SPACING,
            output_dir=tmp_output_dir,
            output_prefix="counts_",
        )
        assert result.success, f"CpptrajDensityAction failed: {result.error}"
        grid = Grid.read_dx(result.output_files[0])
        assert np.count_nonzero(grid.data) > 0, (
            "Expected non-zero WAT counts in the density grid"
        )
