"""Tests for trajectory alignment."""
import pytest

from pymdmix.analysis.align import (
    _convert_mask_to_mda,
    AlignmentResult,
)


class TestMaskConversion:
    """Tests for Amber mask to MDAnalysis conversion."""
    
    def test_atom_mask(self):
        """Test converting atom mask."""
        result = _convert_mask_to_mda("@CA,C,N")
        assert result == "name CA C N"
    
    def test_single_atom(self):
        """Test single atom mask."""
        result = _convert_mask_to_mda("@CA")
        assert result == "name CA"
    
    def test_residue_mask(self):
        """Test residue mask."""
        result = _convert_mask_to_mda(":WAT")
        assert result == "resname WAT"
    
    def test_multiple_residues(self):
        """Test multiple residue mask."""
        result = _convert_mask_to_mda(":WAT,ETA")
        assert result == "resname WAT ETA"
    
    def test_backbone_keyword(self):
        """Test backbone keyword."""
        result = _convert_mask_to_mda("backbone")
        assert "CA" in result
        assert "C" in result
        assert "N" in result
    
    def test_mda_syntax_passthrough(self):
        """Test that MDAnalysis syntax passes through unchanged."""
        result = _convert_mask_to_mda("protein and name CA")
        assert result == "protein and name CA"
        
        result = _convert_mask_to_mda("name CA C N")
        assert result == "name CA C N"


class TestAlignmentResult:
    """Tests for AlignmentResult dataclass."""
    
    def test_repr(self):
        """Test string representation."""
        from pathlib import Path
        
        result = AlignmentResult(
            output_trajectory=Path("aligned.nc"),
            n_frames=100,
            rmsd_mean=1.234,
            rmsd_std=0.567,
            method="mdanalysis",
        )
        
        s = repr(result)
        assert "100" in s
        assert "1.234" in s
        assert "0.567" in s


# Note: Full alignment tests would require trajectory files
# These are covered in integration tests with mock data
class TestAlignmentFunction:
    """Tests for alignment functions that don't need files."""
    
    def test_alignment_result_fields(self):
        """Test AlignmentResult has all required fields."""
        from pathlib import Path
        
        result = AlignmentResult(
            output_trajectory=Path("out.nc"),
            n_frames=50,
            rmsd_mean=2.0,
            rmsd_std=0.5,
            method="cpptraj",
        )
        
        assert result.output_trajectory == Path("out.nc")
        assert result.n_frames == 50
        assert result.rmsd_mean == 2.0
        assert result.rmsd_std == 0.5
        assert result.method == "cpptraj"
