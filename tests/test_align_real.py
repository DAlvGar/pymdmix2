"""
Tests for trajectory alignment using the real Amber NetCDF test data.

These tests exercise ``align_trajectory`` (MDAnalysis backend) against the
real solvated-peptide trajectory in ``tests/data/amber/``.

Marked ``@pytest.mark.real_data`` — automatically skipped when Git LFS data
is not checked out.  Run locally after::

    git lfs pull
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from pymdmix.analysis.align import AlignmentResult, _convert_mask_to_mda, align_trajectory

# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------


def _is_real_trajectory(path: Path) -> bool:
    """Return False if the file is an LFS pointer rather than real binary data."""
    return path.exists() and path.stat().st_size > 100_000


# ===========================================================================
# 1.  Alignment with real trajectory — MDAnalysis backend
# ===========================================================================


@pytest.mark.real_data
class TestAlignmentMDAnalysis:
    """Align traj.nc to pep_WAT_WAT_1.pdb and verify statistics."""

    def test_align_returns_result(
        self, solvated_pep_pdb_path, amber_traj_nc_path, tmp_output_dir
    ):
        """align_trajectory returns an AlignmentResult instance."""
        if not _is_real_trajectory(amber_traj_nc_path):
            pytest.skip("traj.nc is an LFS pointer — run: git lfs pull")
        out = tmp_output_dir / "aligned.nc"
        result = align_trajectory(
            topology=solvated_pep_pdb_path,
            trajectory=amber_traj_nc_path,
            output=out,
            reference=solvated_pep_pdb_path,
            mask="protein and name CA",
            method="mdanalysis",
        )
        assert isinstance(result, AlignmentResult)

    def test_aligned_output_file_created(
        self, solvated_pep_pdb_path, amber_traj_nc_path, tmp_output_dir
    ):
        """Output trajectory file is written to disk."""
        if not _is_real_trajectory(amber_traj_nc_path):
            pytest.skip("traj.nc is an LFS pointer — run: git lfs pull")
        out = tmp_output_dir / "aligned.nc"
        align_trajectory(
            topology=solvated_pep_pdb_path,
            trajectory=amber_traj_nc_path,
            output=out,
            reference=solvated_pep_pdb_path,
            mask="protein and name CA",
            method="mdanalysis",
        )
        assert out.exists()
        assert out.stat().st_size > 0

    def test_frame_count_preserved(
        self, solvated_pep_pdb_path, amber_traj_nc_path, tmp_output_dir
    ):
        """Aligned trajectory has the same number of frames as the input."""
        if not _is_real_trajectory(amber_traj_nc_path):
            pytest.skip("traj.nc is an LFS pointer — run: git lfs pull")
        out = tmp_output_dir / "aligned.nc"
        result = align_trajectory(
            topology=solvated_pep_pdb_path,
            trajectory=amber_traj_nc_path,
            output=out,
            reference=solvated_pep_pdb_path,
            mask="protein and name CA",
            method="mdanalysis",
        )
        assert result.n_frames == 20

    def test_rmsd_mean_is_positive_finite(
        self, solvated_pep_pdb_path, amber_traj_nc_path, tmp_output_dir
    ):
        """Mean backbone RMSD is a positive finite number."""
        if not _is_real_trajectory(amber_traj_nc_path):
            pytest.skip("traj.nc is an LFS pointer — run: git lfs pull")
        out = tmp_output_dir / "aligned.nc"
        result = align_trajectory(
            topology=solvated_pep_pdb_path,
            trajectory=amber_traj_nc_path,
            output=out,
            reference=solvated_pep_pdb_path,
            mask="protein and name CA",
            method="mdanalysis",
        )
        assert np.isfinite(result.rmsd_mean)
        assert result.rmsd_mean > 0

    def test_rmsd_within_sensible_range(
        self, solvated_pep_pdb_path, amber_traj_nc_path, tmp_output_dir
    ):
        """
        Mean Cα RMSD < 5 Å for a short peptide around its reference.

        This acts as a regression guard — unexpectedly large RMSD would
        indicate a bug in the alignment code or a wrong reference.
        """
        if not _is_real_trajectory(amber_traj_nc_path):
            pytest.skip("traj.nc is an LFS pointer — run: git lfs pull")
        out = tmp_output_dir / "aligned.nc"
        result = align_trajectory(
            topology=solvated_pep_pdb_path,
            trajectory=amber_traj_nc_path,
            output=out,
            reference=solvated_pep_pdb_path,
            mask="protein and name CA",
            method="mdanalysis",
        )
        assert result.rmsd_mean < 5.0, (
            f"RMSD mean {result.rmsd_mean:.3f} Å is unexpectedly large"
        )

    def test_rmsd_std_is_non_negative(
        self, solvated_pep_pdb_path, amber_traj_nc_path, tmp_output_dir
    ):
        """RMSD standard deviation is non-negative."""
        if not _is_real_trajectory(amber_traj_nc_path):
            pytest.skip("traj.nc is an LFS pointer — run: git lfs pull")
        out = tmp_output_dir / "aligned.nc"
        result = align_trajectory(
            topology=solvated_pep_pdb_path,
            trajectory=amber_traj_nc_path,
            output=out,
            reference=solvated_pep_pdb_path,
            mask="protein and name CA",
            method="mdanalysis",
        )
        assert result.rmsd_std >= 0

    def test_method_label_is_mdanalysis(
        self, solvated_pep_pdb_path, amber_traj_nc_path, tmp_output_dir
    ):
        """Result records the 'mdanalysis' method label."""
        if not _is_real_trajectory(amber_traj_nc_path):
            pytest.skip("traj.nc is an LFS pointer — run: git lfs pull")
        out = tmp_output_dir / "aligned.nc"
        result = align_trajectory(
            topology=solvated_pep_pdb_path,
            trajectory=amber_traj_nc_path,
            output=out,
            reference=solvated_pep_pdb_path,
            mask="protein and name CA",
            method="mdanalysis",
        )
        assert result.method == "mdanalysis"

    def test_mean_rmsd_alias_equals_rmsd_mean(
        self, solvated_pep_pdb_path, amber_traj_nc_path, tmp_output_dir
    ):
        """AlignmentResult.mean_rmsd is an alias for rmsd_mean."""
        if not _is_real_trajectory(amber_traj_nc_path):
            pytest.skip("traj.nc is an LFS pointer — run: git lfs pull")
        out = tmp_output_dir / "aligned.nc"
        result = align_trajectory(
            topology=solvated_pep_pdb_path,
            trajectory=amber_traj_nc_path,
            output=out,
            reference=solvated_pep_pdb_path,
            mask="protein and name CA",
            method="mdanalysis",
        )
        assert result.mean_rmsd == result.rmsd_mean

    def test_align_without_explicit_reference(
        self, solvated_pep_pdb_path, amber_traj_nc_path, tmp_output_dir
    ):
        """Alignment uses first frame as reference when reference=None."""
        if not _is_real_trajectory(amber_traj_nc_path):
            pytest.skip("traj.nc is an LFS pointer — run: git lfs pull")
        out = tmp_output_dir / "aligned_noref.nc"
        result = align_trajectory(
            topology=solvated_pep_pdb_path,
            trajectory=amber_traj_nc_path,
            output=out,
            reference=None,
            mask="protein and name CA",
            method="mdanalysis",
        )
        # First frame RMSD to itself is 0
        assert result.n_frames == 20
        assert result.rmsd_mean >= 0

    def test_align_subset_frames(
        self, solvated_pep_pdb_path, amber_traj_nc_path, tmp_output_dir
    ):
        """start/stop parameters restrict the frames that are aligned."""
        if not _is_real_trajectory(amber_traj_nc_path):
            pytest.skip("traj.nc is an LFS pointer — run: git lfs pull")
        out = tmp_output_dir / "aligned_sub.nc"
        result = align_trajectory(
            topology=solvated_pep_pdb_path,
            trajectory=amber_traj_nc_path,
            output=out,
            reference=solvated_pep_pdb_path,
            mask="protein and name CA",
            method="mdanalysis",
            start=0,
            stop=10,
        )
        assert result.n_frames == 10

    def test_aligned_trajectory_is_readable(
        self, solvated_pep_pdb_path, amber_traj_nc_path, tmp_output_dir
    ):
        """The output aligned trajectory can be re-opened with open_trajectory."""
        if not _is_real_trajectory(amber_traj_nc_path):
            pytest.skip("traj.nc is an LFS pointer — run: git lfs pull")
        from pymdmix.core.trajectory import open_trajectory

        out = tmp_output_dir / "aligned_reread.nc"
        align_trajectory(
            topology=solvated_pep_pdb_path,
            trajectory=amber_traj_nc_path,
            output=out,
            reference=solvated_pep_pdb_path,
            mask="protein and name CA",
            method="mdanalysis",
        )
        traj = open_trajectory(solvated_pep_pdb_path, out)
        assert traj.n_frames == 20


# ===========================================================================
# 2.  Mask conversion (unit tests — no real data required)
# ===========================================================================


class TestMaskConversionExtended:
    """Additional mask conversion tests covering edge cases."""

    def test_multi_atom_amber_mask(self):
        assert _convert_mask_to_mda("@CA,C,N,O") == "name CA C N O"

    def test_colon_single_residue(self):
        assert _convert_mask_to_mda(":WAT") == "resname WAT"

    def test_colon_multiple_residues(self):
        result = _convert_mask_to_mda(":WAT,ETA,MAM")
        assert result == "resname WAT ETA MAM"

    def test_backbone_keyword_contains_ca(self):
        result = _convert_mask_to_mda("backbone")
        assert "CA" in result

    def test_mda_passthrough_with_and(self):
        sel = "protein and name CA"
        assert _convert_mask_to_mda(sel) == sel

    def test_mda_passthrough_with_or(self):
        sel = "resname WAT or resname ETA"
        assert _convert_mask_to_mda(sel) == sel

    def test_name_passthrough(self):
        sel = "name O1 O2"
        assert _convert_mask_to_mda(sel) == sel
