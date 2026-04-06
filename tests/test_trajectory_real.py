"""
Tests for real trajectory reading using the Amber NetCDF test data.

These tests exercise both the MDAnalysis and native Amber NetCDF backends
against the real solvated-peptide trajectory included in ``tests/data/amber/``.

Marked ``@pytest.mark.real_data`` — automatically skipped when Git LFS data
is not checked out (``traj.nc`` absent).  Run locally after::

    git lfs pull
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from pymdmix.core.trajectory import (
    AmberNetCDFReader,
    Frame,
    FrameSliceReader,
    MDAnalysisReader,
    TrajectoryReader,
    open_trajectory,
)

# ---------------------------------------------------------------------------
# Expected constants for the bundled trajectory (pep_WAT_WAT_1 + traj.nc)
# ---------------------------------------------------------------------------

_EXPECTED_FRAMES = 20
_EXPECTED_ATOMS = 7079  # 8-residue peptide + 2319 WAT + H atoms
_EXPECTED_PROTEIN_ATOMS = 122


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _is_real_trajectory(path: Path) -> bool:
    """Return False if LFS pointer (text) instead of real binary."""
    return path.exists() and path.stat().st_size > 100_000


# ===========================================================================
# 1.  MDAnalysis backend
# ===========================================================================


@pytest.mark.real_data
class TestMDAnalysisBackend:
    """Open traj.nc via MDAnalysis and verify frame / atom counts."""

    def test_open_trajectory_mdanalysis(self, solvated_pep_pdb_path, amber_traj_nc_path):
        """Trajectory opens without error using the MDAnalysis backend."""
        if not _is_real_trajectory(amber_traj_nc_path):
            pytest.skip("traj.nc is an LFS pointer — run: git lfs pull")
        traj = open_trajectory(solvated_pep_pdb_path, amber_traj_nc_path, backend="mdanalysis")
        assert isinstance(traj, MDAnalysisReader)

    def test_frame_count(self, solvated_pep_pdb_path, amber_traj_nc_path):
        """Trajectory has exactly 20 frames."""
        if not _is_real_trajectory(amber_traj_nc_path):
            pytest.skip("traj.nc is an LFS pointer — run: git lfs pull")
        traj = open_trajectory(solvated_pep_pdb_path, amber_traj_nc_path, backend="mdanalysis")
        assert traj.n_frames == _EXPECTED_FRAMES
        assert len(traj) == _EXPECTED_FRAMES

    def test_atom_count(self, solvated_pep_pdb_path, amber_traj_nc_path):
        """Trajectory has the expected number of atoms (peptide + WAT)."""
        if not _is_real_trajectory(amber_traj_nc_path):
            pytest.skip("traj.nc is an LFS pointer — run: git lfs pull")
        traj = open_trajectory(solvated_pep_pdb_path, amber_traj_nc_path, backend="mdanalysis")
        assert traj.n_atoms == _EXPECTED_ATOMS

    def test_iterate_all_frames(self, solvated_pep_pdb_path, amber_traj_nc_path):
        """Iterating over all frames yields Frame objects with correct shapes."""
        if not _is_real_trajectory(amber_traj_nc_path):
            pytest.skip("traj.nc is an LFS pointer — run: git lfs pull")
        traj = open_trajectory(solvated_pep_pdb_path, amber_traj_nc_path, backend="mdanalysis")
        frames = list(traj)
        assert len(frames) == _EXPECTED_FRAMES
        for frame in frames:
            assert isinstance(frame, Frame)
            assert frame.coordinates.shape == (_EXPECTED_ATOMS, 3)
            assert frame.coordinates.dtype == np.float64

    def test_frame_coordinates_are_finite(self, solvated_pep_pdb_path, amber_traj_nc_path):
        """No NaN or Inf in any frame's coordinates."""
        if not _is_real_trajectory(amber_traj_nc_path):
            pytest.skip("traj.nc is an LFS pointer — run: git lfs pull")
        traj = open_trajectory(solvated_pep_pdb_path, amber_traj_nc_path, backend="mdanalysis")
        for frame in traj:
            assert np.all(np.isfinite(frame.coordinates)), "Non-finite coordinates in frame"

    def test_frame_has_box_dimensions(self, solvated_pep_pdb_path, amber_traj_nc_path):
        """Each frame carries periodic-box information."""
        if not _is_real_trajectory(amber_traj_nc_path):
            pytest.skip("traj.nc is an LFS pointer — run: git lfs pull")
        traj = open_trajectory(solvated_pep_pdb_path, amber_traj_nc_path, backend="mdanalysis")
        frame = next(iter(traj))
        assert frame.box is not None
        assert frame.box.shape == (6,)
        # Box lengths must be positive
        assert np.all(frame.box[:3] > 0)

    def test_coordinates_vary_across_frames(self, solvated_pep_pdb_path, amber_traj_nc_path):
        """Frames are not identical (trajectory is not a static snapshot)."""
        if not _is_real_trajectory(amber_traj_nc_path):
            pytest.skip("traj.nc is an LFS pointer — run: git lfs pull")
        traj = open_trajectory(solvated_pep_pdb_path, amber_traj_nc_path, backend="mdanalysis")
        frames = list(traj)
        assert not np.allclose(frames[0].coordinates, frames[-1].coordinates), (
            "First and last frames are identical — trajectory may not have loaded correctly"
        )

    def test_select_atoms_protein(self, solvated_pep_pdb_path, amber_traj_nc_path):
        """MDAnalysisReader.select_atoms returns correct protein indices."""
        if not _is_real_trajectory(amber_traj_nc_path):
            pytest.skip("traj.nc is an LFS pointer — run: git lfs pull")
        traj = open_trajectory(solvated_pep_pdb_path, amber_traj_nc_path, backend="mdanalysis")
        prot_idx = traj.select_atoms("protein")
        assert len(prot_idx) == _EXPECTED_PROTEIN_ATOMS

    def test_repr(self, solvated_pep_pdb_path, amber_traj_nc_path):
        """repr() returns a human-readable string."""
        if not _is_real_trajectory(amber_traj_nc_path):
            pytest.skip("traj.nc is an LFS pointer — run: git lfs pull")
        traj = open_trajectory(solvated_pep_pdb_path, amber_traj_nc_path, backend="mdanalysis")
        s = repr(traj)
        assert "traj.nc" in s


# ===========================================================================
# 2.  Native Amber NetCDF backend
# ===========================================================================


@pytest.mark.real_data
class TestAmberNetCDFBackend:
    """Open traj.nc with the minimal netCDF4-based reader."""

    def test_open_trajectory_amber(self, solvated_pep_pdb_path, amber_traj_nc_path):
        """Trajectory opens without error using the Amber backend."""
        if not _is_real_trajectory(amber_traj_nc_path):
            pytest.skip("traj.nc is an LFS pointer — run: git lfs pull")
        traj = open_trajectory(solvated_pep_pdb_path, amber_traj_nc_path, backend="amber")
        assert isinstance(traj, AmberNetCDFReader)

    def test_frame_count_amber(self, solvated_pep_pdb_path, amber_traj_nc_path):
        """Amber reader reports the same frame count as MDAnalysis."""
        if not _is_real_trajectory(amber_traj_nc_path):
            pytest.skip("traj.nc is an LFS pointer — run: git lfs pull")
        traj = open_trajectory(solvated_pep_pdb_path, amber_traj_nc_path, backend="amber")
        assert traj.n_frames == _EXPECTED_FRAMES

    def test_atom_count_amber(self, solvated_pep_pdb_path, amber_traj_nc_path):
        """Amber reader reports the correct atom count."""
        if not _is_real_trajectory(amber_traj_nc_path):
            pytest.skip("traj.nc is an LFS pointer — run: git lfs pull")
        traj = open_trajectory(solvated_pep_pdb_path, amber_traj_nc_path, backend="amber")
        assert traj.n_atoms == _EXPECTED_ATOMS

    def test_iterate_frames_amber(self, solvated_pep_pdb_path, amber_traj_nc_path):
        """Amber reader yields frames with correct shape."""
        if not _is_real_trajectory(amber_traj_nc_path):
            pytest.skip("traj.nc is an LFS pointer — run: git lfs pull")
        traj = open_trajectory(solvated_pep_pdb_path, amber_traj_nc_path, backend="amber")
        for frame in traj:
            assert frame.coordinates.shape == (_EXPECTED_ATOMS, 3)

    def test_amber_and_mdanalysis_agree_on_coordinates(
        self, solvated_pep_pdb_path, amber_traj_nc_path
    ):
        """Both backends return identical first-frame coordinates."""
        if not _is_real_trajectory(amber_traj_nc_path):
            pytest.skip("traj.nc is an LFS pointer — run: git lfs pull")
        traj_mda = open_trajectory(solvated_pep_pdb_path, amber_traj_nc_path, backend="mdanalysis")
        traj_nc = open_trajectory(solvated_pep_pdb_path, amber_traj_nc_path, backend="amber")
        frame_mda = next(iter(traj_mda))
        frame_nc = next(iter(traj_nc))
        np.testing.assert_allclose(frame_mda.coordinates, frame_nc.coordinates, atol=1e-3)

    def test_frame_box_amber(self, solvated_pep_pdb_path, amber_traj_nc_path):
        """Amber reader provides periodic box from cell_lengths/cell_angles."""
        if not _is_real_trajectory(amber_traj_nc_path):
            pytest.skip("traj.nc is an LFS pointer — run: git lfs pull")
        traj = open_trajectory(solvated_pep_pdb_path, amber_traj_nc_path, backend="amber")
        frame = next(iter(traj))
        assert frame.box is not None
        assert frame.box.shape == (6,)
        assert np.all(frame.box[:3] > 0)


# ===========================================================================
# 3.  Auto-detect backend
# ===========================================================================


@pytest.mark.real_data
class TestAutoBackend:
    """open_trajectory() auto-selects a backend based on available libraries."""

    def test_auto_opens_nc(self, solvated_pep_pdb_path, amber_traj_nc_path):
        """Auto mode opens a .nc file without an explicit backend specification."""
        if not _is_real_trajectory(amber_traj_nc_path):
            pytest.skip("traj.nc is an LFS pointer — run: git lfs pull")
        traj = open_trajectory(solvated_pep_pdb_path, amber_traj_nc_path)
        assert isinstance(traj, TrajectoryReader)
        assert traj.n_frames == _EXPECTED_FRAMES

    def test_auto_satisfies_protocol(self, solvated_pep_pdb_path, amber_traj_nc_path):
        """Auto-detected reader satisfies the TrajectoryReader protocol."""
        if not _is_real_trajectory(amber_traj_nc_path):
            pytest.skip("traj.nc is an LFS pointer — run: git lfs pull")
        traj = open_trajectory(solvated_pep_pdb_path, amber_traj_nc_path)
        assert isinstance(traj, TrajectoryReader)


# ===========================================================================
# 4.  FrameSliceReader
# ===========================================================================


@pytest.mark.real_data
class TestFrameSliceReader:
    """FrameSliceReader correctly subsets a real trajectory."""

    def test_slice_first_half(self, solvated_pep_pdb_path, amber_traj_nc_path):
        """Slicing first 10 frames yields exactly 10 frames."""
        if not _is_real_trajectory(amber_traj_nc_path):
            pytest.skip("traj.nc is an LFS pointer — run: git lfs pull")
        base = open_trajectory(solvated_pep_pdb_path, amber_traj_nc_path)
        sliced = FrameSliceReader(base, start=0, stop=10)
        assert sliced.n_frames == 10
        assert len(list(sliced)) == 10

    def test_slice_with_step(self, solvated_pep_pdb_path, amber_traj_nc_path):
        """Step=2 on 20 frames gives 10 frames."""
        if not _is_real_trajectory(amber_traj_nc_path):
            pytest.skip("traj.nc is an LFS pointer — run: git lfs pull")
        base = open_trajectory(solvated_pep_pdb_path, amber_traj_nc_path)
        sliced = FrameSliceReader(base, step=2)
        assert sliced.n_frames == 10
        assert len(list(sliced)) == 10

    def test_slice_inherits_n_atoms(self, solvated_pep_pdb_path, amber_traj_nc_path):
        """Sliced reader inherits atom count from underlying reader."""
        if not _is_real_trajectory(amber_traj_nc_path):
            pytest.skip("traj.nc is an LFS pointer — run: git lfs pull")
        base = open_trajectory(solvated_pep_pdb_path, amber_traj_nc_path)
        sliced = FrameSliceReader(base, stop=5)
        assert sliced.n_atoms == _EXPECTED_ATOMS
