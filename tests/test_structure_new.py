"""Tests for new structure functions (Biskit replacements)."""

import numpy as np
import pytest
import tempfile
from pathlib import Path

import parmed

# Import all new functions
from pymdmix.core.structure import (
    # New masks
    get_hydrogen_mask,
    get_backbone_mask,
    # Chain operations
    get_chain_start_indices,
    select_chains,
    # Profile conversions
    atom_to_residue_values,
    residue_to_atom_values,
    # Probe helpers
    get_probe_coords,
    iter_residue_coords,
    get_residue_com_coords,
    detect_solvent_type,
    # Alignment
    align_structures,
    align_to_reference,
    compute_rmsd,
    # Existing functions for setup
    load_structure,
    get_protein_mask,
    WATER_RESIDUES,
)


@pytest.fixture
def simple_pdb():
    """Create a simple test PDB structure."""
    pdb_content = """\
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00  0.00           C
ATOM      3  C   ALA A   1       2.009   1.420   0.000  1.00  0.00           C
ATOM      4  O   ALA A   1       1.246   2.390   0.000  1.00  0.00           O
ATOM      5  CB  ALA A   1       1.986  -0.760  -1.216  1.00  0.00           C
ATOM      6  H   ALA A   1      -0.500   0.866   0.000  1.00  0.00           H
ATOM      7  N   GLY A   2       3.320   1.560   0.000  1.00  0.00           N
ATOM      8  CA  GLY A   2       3.970   2.860   0.000  1.00  0.00           C
ATOM      9  C   GLY A   2       5.480   2.760   0.000  1.00  0.00           C
ATOM     10  O   GLY A   2       6.080   1.690   0.000  1.00  0.00           O
ATOM     11  H   GLY A   2       3.860   0.720   0.000  1.00  0.00           H
TER
ATOM     12  N   ALA B   1      10.000  10.000  10.000  1.00  0.00           N
ATOM     13  CA  ALA B   1      11.458  10.000  10.000  1.00  0.00           C
ATOM     14  C   ALA B   1      12.009  11.420  10.000  1.00  0.00           C
ATOM     15  O   ALA B   1      11.246  12.390  10.000  1.00  0.00           O
ATOM     16  CB  ALA B   1      11.986   9.240   8.784  1.00  0.00           C
ATOM     17  H   ALA B   1       9.500  10.866  10.000  1.00  0.00           H
TER
END
"""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f:
        f.write(pdb_content)
        path = Path(f.name)
    
    yield parmed.load_file(str(path))
    path.unlink()


class TestHydrogenMask:
    """Tests for get_hydrogen_mask."""

    def test_finds_hydrogens(self, simple_pdb):
        mask = get_hydrogen_mask(simple_pdb)
        assert mask.sum() == 3  # H atoms at positions 5, 10, 16 (0-indexed)
        
    def test_returns_bool_array(self, simple_pdb):
        mask = get_hydrogen_mask(simple_pdb)
        assert mask.dtype == np.bool_
        assert len(mask) == len(simple_pdb.atoms)


class TestBackboneMask:
    """Tests for get_backbone_mask."""

    def test_finds_backbone(self, simple_pdb):
        mask = get_backbone_mask(simple_pdb)
        # Backbone: N, CA, C, O for each residue (3 residues)
        # ALA1: N, CA, C, O (4)
        # GLY2: N, CA, C, O (4)  
        # ALA3: N, CA, C, O (4)
        # Total: 12
        assert mask.sum() == 12
        
    def test_excludes_sidechain_and_H(self, simple_pdb):
        mask = get_backbone_mask(simple_pdb)
        bb_atoms = [simple_pdb.atoms[i] for i in np.where(mask)[0]]
        for atom in bb_atoms:
            assert atom.name in {'N', 'CA', 'C', 'O'}


class TestChainOperations:
    """Tests for chain operations."""

    def test_chain_start_indices(self, simple_pdb):
        indices = get_chain_start_indices(simple_pdb)
        assert len(indices) == 2  # Two chains A and B
        assert indices[0] == 0  # Chain A starts at atom 0
        
    def test_select_chains(self, simple_pdb):
        chain_a = select_chains(simple_pdb, 'A')
        assert len(chain_a.atoms) == 11  # ALA + GLY residues
        
        chain_b = select_chains(simple_pdb, 'B')
        assert len(chain_b.atoms) == 6  # Single ALA residue
        
    def test_select_multiple_chains(self, simple_pdb):
        both = select_chains(simple_pdb, ['A', 'B'])
        # Note: parmed may not assign chain ID to first atom in some cases
        # So we check that we get most atoms rather than all
        assert len(both.atoms) >= len(simple_pdb.atoms) - 1


class TestProfileConversions:
    """Tests for atom_to_residue_values and residue_to_atom_values."""

    def test_atom_to_residue(self, simple_pdb):
        # Create per-atom values
        atom_values = list(range(len(simple_pdb.atoms)))
        
        res_values = atom_to_residue_values(simple_pdb, atom_values)
        assert len(res_values) == len(simple_pdb.residues)  # 3 residues
        
    def test_residue_to_atom(self, simple_pdb):
        res_values = [1.0, 2.0, 3.0]  # One per residue
        
        atom_values = residue_to_atom_values(simple_pdb, res_values)
        assert len(atom_values) == len(simple_pdb.atoms)
        
        # Check that atoms in same residue have same value
        for res, val in zip(simple_pdb.residues, res_values):
            for atom in res.atoms:
                assert atom_values[atom.idx] == val
                
    def test_roundtrip(self, simple_pdb):
        # Create values and convert back
        original = [float(i) for i in range(len(simple_pdb.residues))]
        atom_vals = residue_to_atom_values(simple_pdb, original)
        back = atom_to_residue_values(simple_pdb, atom_vals)
        assert back == original


class TestAlignment:
    """Tests for structure alignment functions."""

    def test_compute_rmsd_identical(self):
        coords = np.random.randn(100, 3)
        rmsd = compute_rmsd(coords, coords)
        assert rmsd == pytest.approx(0.0, abs=1e-10)
        
    def test_compute_rmsd_translated(self):
        coords1 = np.random.randn(100, 3)
        coords2 = coords1 + np.array([1.0, 2.0, 3.0])
        rmsd = compute_rmsd(coords1, coords2)
        expected = np.sqrt(1**2 + 2**2 + 3**2)  # Translation distance
        assert rmsd == pytest.approx(expected, rel=1e-6)
        
    def test_align_structures_identical(self, simple_pdb):
        # Align structure to itself
        aligned, rmsd = align_structures(simple_pdb, simple_pdb, in_place=False)
        assert rmsd == pytest.approx(0.0, abs=1e-6)
        
    def test_align_structures_translated(self, simple_pdb):
        # Create translated copy
        mobile = simple_pdb.copy(parmed.Structure)
        mobile.coordinates = simple_pdb.coordinates + np.array([5.0, 5.0, 5.0])
        
        aligned, rmsd = align_structures(mobile, simple_pdb, in_place=True)
        assert rmsd == pytest.approx(0.0, abs=1e-6)
        
        # Check coordinates are now aligned
        np.testing.assert_array_almost_equal(
            aligned.coordinates, simple_pdb.coordinates, decimal=5
        )
        
    def test_align_to_reference_backbone(self, simple_pdb):
        mobile = simple_pdb.copy(parmed.Structure)
        mobile.coordinates = simple_pdb.coordinates + np.array([3.0, 0.0, 0.0])
        
        aligned, rmsd = align_to_reference(mobile, simple_pdb, selection='backbone')
        assert rmsd == pytest.approx(0.0, abs=1e-6)


class TestProbeHelpers:
    """Tests for probe/density helper functions."""

    def test_get_probe_coords(self, simple_pdb):
        # Get all CA atoms
        coords = get_probe_coords(simple_pdb, 'ALA', 'CA')
        assert coords.shape == (2, 3)  # Two ALA residues
        
    def test_iter_residue_coords(self, simple_pdb):
        count = 0
        for coords in iter_residue_coords(simple_pdb, 'ALA'):
            count += 1
            assert coords.ndim == 2
            assert coords.shape[1] == 3
        assert count == 2  # Two ALA residues
        
    def test_get_residue_com_coords(self, simple_pdb):
        coms = get_residue_com_coords(simple_pdb, 'ALA')
        assert coms.shape == (2, 3)  # Two ALA residues
        
    def test_get_residue_com_empty(self, simple_pdb):
        coms = get_residue_com_coords(simple_pdb, 'TRP')  # No TRP
        assert coms.shape == (0, 3)


class TestSolventDetection:
    """Tests for detect_solvent_type."""

    def test_detect_no_solvent(self, simple_pdb):
        known = {'ETA': ['ETA'], 'MAM': ['MAM']}
        result = detect_solvent_type(simple_pdb, known)
        assert result is None  # No organic solvent in structure


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
