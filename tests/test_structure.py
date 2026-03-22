"""
Tests for pymdmix.core.structure module.
"""

import numpy as np
import pytest

from pymdmix.core.structure import (
    PROTEIN_RESIDUES,
    WATER_RESIDUES,
    count_residues,
    find_disulfides,
    get_atom_mask,
    get_heavy_atom_mask,
    get_protein_mask,
    get_residue_mask,
    get_water_mask,
    load_structure,
)


class TestResidueSets:
    """Test residue set definitions."""

    def test_protein_residues_complete(self):
        """Test protein residues include standard amino acids."""
        standard_aa = {
            "ALA",
            "ARG",
            "ASN",
            "ASP",
            "CYS",
            "GLN",
            "GLU",
            "GLY",
            "HIS",
            "ILE",
            "LEU",
            "LYS",
            "MET",
            "PHE",
            "PRO",
            "SER",
            "THR",
            "TRP",
            "TYR",
            "VAL",
        }

        assert standard_aa.issubset(PROTEIN_RESIDUES)

    def test_protein_residues_include_variants(self):
        """Test protein residues include Amber variants."""
        amber_variants = {"CYX", "HID", "HIE", "HIP", "ACE", "NME"}

        assert amber_variants.issubset(PROTEIN_RESIDUES)

    def test_water_residues(self):
        """Test water residue names."""
        common_water = {"WAT", "HOH", "TIP3"}

        assert common_water.issubset(WATER_RESIDUES)


class TestLoadStructure:
    """Test structure loading."""

    def test_load_pdb(self, sample_pdb_file):
        """Test loading PDB file."""
        struct = load_structure(sample_pdb_file)

        assert struct is not None
        assert len(struct.atoms) == 13

    def test_load_nonexistent(self, tmp_output_dir):
        """Test loading non-existent file raises error."""
        with pytest.raises(FileNotFoundError):
            load_structure(tmp_output_dir / "nonexistent.pdb")


class TestMasks:
    """Test atom mask functions."""

    def test_protein_mask(self, sample_pdb_file):
        """Test protein mask generation."""
        struct = load_structure(sample_pdb_file)
        mask = get_protein_mask(struct)

        assert mask.dtype == np.bool_
        assert len(mask) == len(struct.atoms)
        # Should match ALA and GLY atoms (10 atoms)
        assert mask.sum() == 10

    def test_water_mask(self, sample_pdb_file):
        """Test water mask generation."""
        struct = load_structure(sample_pdb_file)
        mask = get_water_mask(struct)

        # Should match WAT atoms (3 atoms)
        assert mask.sum() == 3

    def test_residue_mask_single(self, sample_pdb_file):
        """Test residue mask for single residue."""
        struct = load_structure(sample_pdb_file)
        mask = get_residue_mask(struct, "ALA")

        # ALA has 5 atoms
        assert mask.sum() == 5

    def test_residue_mask_multiple(self, sample_pdb_file):
        """Test residue mask for multiple residues."""
        struct = load_structure(sample_pdb_file)
        mask = get_residue_mask(struct, ["ALA", "GLY"])

        # ALA (5) + GLY (5) = 10
        assert mask.sum() == 10

    def test_atom_mask(self, sample_pdb_file):
        """Test atom name mask."""
        struct = load_structure(sample_pdb_file)
        mask = get_atom_mask(struct, "CA")

        # 2 CA atoms (ALA + GLY)
        assert mask.sum() == 2

    def test_heavy_atom_mask(self, sample_pdb_file):
        """Test heavy atom mask."""
        struct = load_structure(sample_pdb_file)
        mask = get_heavy_atom_mask(struct)

        # Should exclude H atoms
        h_count = sum(1 for a in struct.atoms if a.name.startswith("H"))
        assert mask.sum() == len(struct.atoms) - h_count


class TestDisulfides:
    """Test disulfide bond detection."""

    def test_find_disulfides_none(self, sample_pdb_file):
        """Test finding disulfides when none present."""
        struct = load_structure(sample_pdb_file)
        disulfides = find_disulfides(struct)

        # Sample PDB has no CYS
        assert disulfides == []

    # Note: Would need a PDB with CYS to test actual disulfide detection


class TestCountResidues:
    """Test residue counting."""

    def test_count_residues(self, sample_pdb_file):
        """Test counting residues."""
        struct = load_structure(sample_pdb_file)
        counts = count_residues(struct)

        assert counts["ALA"] == 1
        assert counts["GLY"] == 1
        assert counts["WAT"] == 1
