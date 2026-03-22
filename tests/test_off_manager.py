"""
Tests for pymdmix.io.off_manager module.
"""

import numpy as np
import pytest

from pymdmix.io.off_manager import (
    Atom,
    OFFManager,
    OFFManagerError,
    OFFSectionError,
    Residue,
)

# Sample OFF file content for testing
SAMPLE_OFF_CONTENT = """!!index array str
 "ETA"
 "WAT"
!entry.ETA.unit.atoms table  str name  str type  int typex  int reession  int residuePdbSequenceNumber  int flags  int sequence  int elmnt  dbl charge
 "C1" "CT" 0 1 1 131072 1 6 -0.1824
 "C2" "CT" 0 1 1 131072 2 6 0.0177
 "O1" "OH" 0 1 1 131072 3 8 -0.6600
 "H1" "H1" 0 1 1 131072 4 1 0.0600
 "H2" "H1" 0 1 1 131072 5 1 0.0600
 "H3" "HC" 0 1 1 131072 6 1 0.0600
 "H4" "HC" 0 1 1 131072 7 1 0.0600
 "H5" "HC" 0 1 1 131072 8 1 0.0600
 "HO" "HO" 0 1 1 131072 9 1 0.4047
!entry.ETA.unit.boundbox array dbl
 -1.000000
 0.0
 0.0
 0.0
 0.0
!entry.ETA.unit.connectivity table  int atom1x  int atom2x  int flags
 1 2 1
 1 6 1
 1 7 1
 1 8 1
 2 3 1
 2 4 1
 2 5 1
 3 9 1
!entry.ETA.unit.hierarchy table  str abession  int residuence  str atomna
 "U" 0 "UNK"
 "R" 1 "ETA"
 "A" 1 "C1"
 "A" 2 "C2"
 "A" 3 "O1"
 "A" 4 "H1"
 "A" 5 "H2"
 "A" 6 "H3"
 "A" 7 "H4"
 "A" 8 "H5"
 "A" 9 "HO"
!entry.ETA.unit.name single str
 "ETA"
!entry.ETA.unit.positions table  dbl x  dbl y  dbl z
 1.000000 0.000000 0.000000
 2.500000 0.000000 0.000000
 3.000000 1.200000 0.000000
 0.500000 1.000000 0.000000
 0.500000 -0.500000 0.866000
 0.500000 -0.500000 -0.866000
 2.800000 -0.500000 0.866000
 2.800000 -0.500000 -0.866000
 3.800000 1.300000 0.500000
!entry.ETA.unit.residueconnect table  int c1x  int c2x  int c3x  int c4x  int c5x  int c6x
 0 0 0 0 0 0
!entry.ETA.unit.residues table  str name  int seq  int childseq  int startato  str restype  int imageid
 "ETA" 1 10 1 "?" 0
!entry.ETA.unit.residuesPdbSequenceNumber array int
 0
!entry.ETA.unit.solventcap array dbl
 -1.000000
 0.0
 0.0
 0.0
 0.0
!entry.ETA.unit.velocities table  dbl x  dbl y  dbl z
 0.0 0.0 0.0
 0.0 0.0 0.0
 0.0 0.0 0.0
 0.0 0.0 0.0
 0.0 0.0 0.0
 0.0 0.0 0.0
 0.0 0.0 0.0
 0.0 0.0 0.0
 0.0 0.0 0.0
!entry.WAT.unit.atoms table  str name  str type  int typex  int reession  int residuePdbSequenceNumber  int flags  int sequence  int elmnt  dbl charge
 "O" "OW" 0 1 1 131072 1 8 -0.8340
 "H1" "HW" 0 1 1 131072 2 1 0.4170
 "H2" "HW" 0 1 1 131072 3 1 0.4170
!entry.WAT.unit.boundbox array dbl
 -1.000000
 0.0
 0.0
 0.0
 0.0
!entry.WAT.unit.connectivity table  int atom1x  int atom2x  int flags
 1 2 1
 1 3 1
!entry.WAT.unit.hierarchy table  str abession  int residuence  str atomna
 "U" 0 "UNK"
 "R" 1 "WAT"
 "A" 1 "O"
 "A" 2 "H1"
 "A" 3 "H2"
!entry.WAT.unit.name single str
 "WAT"
!entry.WAT.unit.positions table  dbl x  dbl y  dbl z
 0.000000 0.000000 0.117000
 0.756000 0.000000 -0.469000
 -0.756000 0.000000 -0.469000
!entry.WAT.unit.residueconnect table  int c1x  int c2x  int c3x  int c4x  int c5x  int c6x
 0 0 0 0 0 0
!entry.WAT.unit.residues table  str name  int seq  int childseq  int startato  str restype  int imageid
 "WAT" 1 4 1 "?" 0
!entry.WAT.unit.residuesPdbSequenceNumber array int
 0
!entry.WAT.unit.solventcap array dbl
 -1.000000
 0.0
 0.0
 0.0
 0.0
!entry.WAT.unit.velocities table  dbl x  dbl y  dbl z
 0.0 0.0 0.0
 0.0 0.0 0.0
 0.0 0.0 0.0
"""

# OFF with box dimensions for volume testing
OFF_WITH_BOX = """!!index array str
 "SOLV"
!entry.SOLV.unit.atoms table  str name  str type  int typex  int reession  int residuePdbSequenceNumber  int flags  int sequence  int elmnt  dbl charge
 "O" "OW" 0 1 1 131072 1 8 -0.8340
!entry.SOLV.unit.boundbox array dbl
 1.000000
 90.0
 20.0
 20.0
 20.0
!entry.SOLV.unit.connectivity table  int atom1x  int atom2x  int flags
!entry.SOLV.unit.hierarchy table  str abession  int residuence  str atomna
 "U" 0 "UNK"
!entry.SOLV.unit.name single str
 "SOLV"
!entry.SOLV.unit.positions table  dbl x  dbl y  dbl z
 0.000000 0.000000 0.000000
!entry.SOLV.unit.residueconnect table  int c1x  int c2x  int c3x  int c4x  int c5x  int c6x
 0 0 0 0 0 0
!entry.SOLV.unit.residues table  str name  int seq  int childseq  int startato  str restype  int imageid
 "WAT" 1 2 1 "?" 0
!entry.SOLV.unit.residuesPdbSequenceNumber array int
 0
!entry.SOLV.unit.solventcap array dbl
 -1.000000
 0.0
 0.0
 0.0
 0.0
!entry.SOLV.unit.velocities table  dbl x  dbl y  dbl z
 0.0 0.0 0.0
"""


class TestAtom:
    """Test Atom dataclass."""

    def test_atom_creation(self):
        """Test creating an atom."""
        atom = Atom(id=1, name="C1", atom_type="CT", element=6, charge=-0.1824)

        assert atom.id == 1
        assert atom.name == "C1"
        assert atom.atom_type == "CT"
        assert atom.element == 6
        assert atom.charge == -0.1824

    def test_atom_repr(self):
        """Test atom string representation."""
        atom = Atom(id=1, name="C1", atom_type="CT", element=6, charge=-0.1824)
        repr_str = repr(atom)

        assert "Atom" in repr_str
        assert "C1" in repr_str
        assert "CT" in repr_str

    def test_atom_equality_by_name(self):
        """Test atom comparison by name."""
        atom1 = Atom(id=1, name="C1", atom_type="CT", element=6, charge=-0.1824)
        atom2 = Atom(id=2, name="C1", atom_type="CT", element=6, charge=0.0)

        assert atom1 == atom2
        assert atom1 == "C1"
        assert not (atom1 == "C2")


class TestResidue:
    """Test Residue dataclass."""

    def test_residue_creation(self):
        """Test creating a residue."""
        atoms = [
            Atom(id=1, name="C1", atom_type="CT", element=6, charge=-0.5),
            Atom(id=2, name="C2", atom_type="CT", element=6, charge=0.5),
        ]
        xyz = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]], dtype=np.float32)

        residue = Residue(
            name="TST",
            atoms=atoms,
            connectivity=[(1, 2)],
            xyz=xyz,
        )

        assert residue.name == "TST"
        assert len(residue.atoms) == 2
        assert residue.num_atoms == 2

    def test_residue_charge(self):
        """Test total charge calculation."""
        atoms = [
            Atom(id=1, name="C1", atom_type="CT", element=6, charge=-0.5),
            Atom(id=2, name="C2", atom_type="CT", element=6, charge=0.5),
        ]
        residue = Residue(name="TST", atoms=atoms)

        assert residue.charge == 0.0

    def test_residue_charge_nonzero(self):
        """Test charge calculation with non-zero total."""
        atoms = [
            Atom(id=1, name="N", atom_type="N3", element=7, charge=0.5),
            Atom(id=2, name="H1", atom_type="H", element=1, charge=0.25),
            Atom(id=3, name="H2", atom_type="H", element=1, charge=0.25),
        ]
        residue = Residue(name="NH3", atoms=atoms)

        assert residue.charge == 1.0

    def test_residue_center(self):
        """Test center of mass calculation."""
        xyz = np.array(
            [
                [0.0, 0.0, 0.0],
                [2.0, 0.0, 0.0],
                [1.0, 2.0, 0.0],
            ],
            dtype=np.float32,
        )
        residue = Residue(name="TST", xyz=xyz)

        center = residue.center()
        np.testing.assert_array_almost_equal(center, [1.0, 2 / 3, 0.0])

    def test_residue_atom_lookup(self):
        """Test atom lookup by id and name."""
        atoms = [
            Atom(id=1, name="C1", atom_type="CT", element=6, charge=0.0),
            Atom(id=2, name="O1", atom_type="OH", element=8, charge=0.0),
        ]
        residue = Residue(name="TST", atoms=atoms)

        assert residue.get_atom_by_id(1).name == "C1"
        assert residue.get_atom_by_name("O1").id == 2
        assert residue.get_atom_by_id(99) is None

    def test_residue_equality(self):
        """Test residue comparison by name."""
        res1 = Residue(name="ETA")
        res2 = Residue(name="ETA")
        res3 = Residue(name="WAT")

        assert res1 == res2
        assert res1 == "ETA"
        assert not (res1 == res3)


class TestOFFManager:
    """Test OFFManager class."""

    def test_from_string(self):
        """Test creating from string content."""
        off = OFFManager.from_string(SAMPLE_OFF_CONTENT)

        assert off.content == SAMPLE_OFF_CONTENT

    def test_from_file(self, tmp_path):
        """Test creating from file."""
        off_file = tmp_path / "test.off"
        off_file.write_text(SAMPLE_OFF_CONTENT)

        off = OFFManager.from_file(off_file)

        assert off.content == SAMPLE_OFF_CONTENT

    def test_from_file_not_found(self):
        """Test error when file doesn't exist."""
        with pytest.raises(OFFManagerError, match="not found"):
            OFFManager.from_file("/nonexistent/path.off")

    def test_get_units(self):
        """Test getting unit names."""
        off = OFFManager.from_string(SAMPLE_OFF_CONTENT)

        units = off.get_units()

        assert "ETA" in units
        assert "WAT" in units
        assert len(units) == 2

    def test_has_unit(self):
        """Test checking for unit existence."""
        off = OFFManager.from_string(SAMPLE_OFF_CONTENT)

        assert off.has_unit("ETA")
        assert off.has_unit("WAT")
        assert not off.has_unit("UNKNOWN")

    def test_get_atoms(self):
        """Test getting atom information."""
        off = OFFManager.from_string(SAMPLE_OFF_CONTENT)

        atoms = off.get_atoms("ETA")

        assert len(atoms) == 9
        assert any(a.name == "C1" for a in atoms)
        assert any(a.name == "O1" for a in atoms)

        # Check specific atom properties
        c1 = next(a for a in atoms if a.name == "C1")
        assert c1.atom_type == "CT"
        assert c1.element == 6
        assert c1.charge == pytest.approx(-0.1824)

    def test_get_atoms_skip_hydrogen(self):
        """Test skipping hydrogens."""
        off = OFFManager.from_string(SAMPLE_OFF_CONTENT)

        atoms = off.get_atoms("ETA", skip_h=True)

        # Only C1, C2, O1 (no H atoms)
        assert len(atoms) == 3
        assert all(a.element != 1 for a in atoms)

    def test_get_coords(self):
        """Test getting coordinates."""
        off = OFFManager.from_string(SAMPLE_OFF_CONTENT)

        coords = off.get_coords("ETA")

        assert coords.shape == (9, 3)
        assert coords.dtype == np.float32
        # First atom at (1, 0, 0)
        np.testing.assert_array_almost_equal(coords[0], [1.0, 0.0, 0.0])

    def test_get_connectivity(self):
        """Test getting bond connectivity."""
        off = OFFManager.from_string(SAMPLE_OFF_CONTENT)

        connectivity = off.get_connectivity("ETA")

        # Each bond listed in both directions
        assert (1, 2) in connectivity
        assert (2, 1) in connectivity
        assert (2, 3) in connectivity

    def test_get_residue(self):
        """Test getting a full residue."""
        off = OFFManager.from_string(SAMPLE_OFF_CONTENT)

        residue = off.get_residue("ETA")

        assert residue.name == "ETA"
        assert len(residue.atoms) == 9
        assert residue.xyz.shape == (9, 3)
        assert len(residue.connectivity) > 0

    def test_get_residue_list(self):
        """Test getting residue names in a unit."""
        off = OFFManager.from_string(SAMPLE_OFF_CONTENT)

        # Unique list
        res_unique = off.get_residue_list("ETA", unique=True)
        assert "ETA" in res_unique

        # Full list
        res_all = off.get_residue_list("ETA", unique=False)
        assert "ETA" in res_all

    def test_read_off_section(self):
        """Test reading a specific section."""
        off = OFFManager.from_string(SAMPLE_OFF_CONTENT)

        atoms_section = off.read_off_section("ETA", "atoms")

        assert len(atoms_section) > 0
        assert any("C1" in line for line in atoms_section)

    def test_read_off_section_not_found(self):
        """Test error when section doesn't exist."""
        off = OFFManager.from_string(SAMPLE_OFF_CONTENT)

        with pytest.raises(OFFSectionError):
            off.read_off_section("ETA", "nonexistent")

    def test_get_num_res(self):
        """Test counting residues."""
        off = OFFManager.from_string(SAMPLE_OFF_CONTENT)

        count = off.get_num_res("ETA", "ETA")
        assert count == 1

    def test_get_box_dimensions(self):
        """Test getting box dimensions."""
        off = OFFManager.from_string(OFF_WITH_BOX)

        dims = off.get_box_dimensions("SOLV")

        # Should be [20, 20, 20] (skipping angle and type)
        assert len(dims) == 3
        assert dims == [20.0, 20.0, 20.0]

    def test_get_volume(self):
        """Test volume calculation."""
        off = OFFManager.from_string(OFF_WITH_BOX)

        volume = off.get_volume("SOLV")

        assert volume == pytest.approx(8000.0)  # 20 * 20 * 20

    def test_get_volume_no_box(self):
        """Test volume when no box defined."""
        off = OFFManager.from_string(SAMPLE_OFF_CONTENT)

        # ETA has no meaningful box dimensions
        volume = off.get_volume("ETA")
        # Should return a volume (product of whatever is there) or 0
        assert volume is not None

    def test_write_and_read_back(self, tmp_path):
        """Test writing OFF content to file."""
        off = OFFManager.from_string(SAMPLE_OFF_CONTENT)
        outfile = tmp_path / "output.off"

        written_path = off.write(outfile)

        assert written_path.exists()
        assert written_path.read_text() == SAMPLE_OFF_CONTENT

    def test_write_tmp(self):
        """Test writing to temporary file."""
        off = OFFManager.from_string(SAMPLE_OFF_CONTENT)

        tmp_path = off.write_tmp()

        assert tmp_path.exists()
        assert tmp_path.read_text() == SAMPLE_OFF_CONTENT

        # Clean up
        assert off.clean_tmp()
        assert not tmp_path.exists()

    def test_clean_tmp_no_file(self):
        """Test clean_tmp when no file created."""
        off = OFFManager.from_string(SAMPLE_OFF_CONTENT)

        # Should not raise, just return False
        assert not off.clean_tmp()

    def test_repr(self):
        """Test string representation."""
        off = OFFManager.from_string(SAMPLE_OFF_CONTENT)

        repr_str = repr(off)

        assert "OFFManager" in repr_str
        assert "ETA" in repr_str
        assert "WAT" in repr_str


class TestOFFManagerWater:
    """Test OFFManager with water unit."""

    def test_water_atoms(self):
        """Test getting water atoms."""
        off = OFFManager.from_string(SAMPLE_OFF_CONTENT)

        atoms = off.get_atoms("WAT")

        assert len(atoms) == 3
        oxygen = next(a for a in atoms if a.name == "O")
        assert oxygen.element == 8
        assert oxygen.charge == pytest.approx(-0.8340)

    def test_water_connectivity(self):
        """Test water connectivity."""
        off = OFFManager.from_string(SAMPLE_OFF_CONTENT)

        connectivity = off.get_connectivity("WAT")

        assert (1, 2) in connectivity
        assert (1, 3) in connectivity

    def test_water_charge(self):
        """Test that water is neutral."""
        off = OFFManager.from_string(SAMPLE_OFF_CONTENT)

        residue = off.get_residue("WAT")

        assert residue.charge == pytest.approx(0.0)
