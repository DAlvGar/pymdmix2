"""
Tests for pymdmix.core.system module.
"""

import pytest
import pickle
import tempfile
import os
from pathlib import Path

from pymdmix.core.system import (
    System,
    SolvatedSystem,
    SystemError,
    BadFileError,
    FileLock,
    load_system,
)


# Sample OFF content for testing
SAMPLE_OFF = '''!!index array str
 "TST"
!entry.TST.unit.atoms table  str name  str type  int typex  int reession  int residuePdbSequenceNumber  int flags  int sequence  int elmnt  dbl charge
 "C1" "CT" 0 1 1 131072 1 6 -0.1824
 "C2" "CT" 0 1 1 131072 2 6 0.0177
 "O1" "OH" 0 1 1 131072 3 8 -0.6600
!entry.TST.unit.boundbox array dbl
 -1.000000
 0.0
 0.0
 0.0
 0.0
!entry.TST.unit.connectivity table  int atom1x  int atom2x  int flags
 1 2 1
 2 3 1
!entry.TST.unit.hierarchy table  str abession  int residuence  str atomna
 "U" 0 "UNK"
!entry.TST.unit.name single str
 "TST"
!entry.TST.unit.positions table  dbl x  dbl y  dbl z
 0.000000 0.000000 0.000000
 1.500000 0.000000 0.000000
 2.500000 1.000000 0.000000
!entry.TST.unit.residueconnect table  int c1x  int c2x  int c3x  int c4x  int c5x  int c6x
 0 0 0 0 0 0
!entry.TST.unit.residues table  str name  int seq  int childseq  int startato  str restype  int imageid
 "TST" 1 4 1 "?" 0
!entry.TST.unit.residuesPdbSequenceNumber array int
 0
!entry.TST.unit.solventcap array dbl
 -1.000000
 0.0
 0.0
 0.0
 0.0
!entry.TST.unit.velocities table  dbl x  dbl y  dbl z
 0.0 0.0 0.0
 0.0 0.0 0.0
 0.0 0.0 0.0
'''

# Sample topology content (minimal)
SAMPLE_PRMTOP = '''%VERSION  VERSION_STAMP = V0001.000  DATE = 01/01/24  00:00:00
%FLAG TITLE
%FORMAT(20a4)
TEST TOPOLOGY
%FLAG POINTERS
%FORMAT(10I8)
       3       2       0       0       0       0       0       0       0       0
'''

# Sample coordinate content (minimal)
SAMPLE_PRMCRD = '''TEST
    3
  0.0000000   0.0000000   0.0000000   1.5000000   0.0000000   0.0000000
  2.5000000   1.0000000   0.0000000
'''


class TestFileLock:
    """Test FileLock class."""

    def test_lock_acquire_release(self, tmp_path):
        """Test acquiring and releasing a lock."""
        test_file = tmp_path / "test.txt"
        test_file.write_text("test")

        lock = FileLock(str(test_file))

        assert not lock.is_locked
        lock.acquire()
        assert lock.is_locked

        # Lock file should exist
        assert Path(lock.lockfile).exists()

        lock.release()
        assert not lock.is_locked
        assert not Path(lock.lockfile).exists()

    def test_lock_context_manager(self, tmp_path):
        """Test lock as context manager."""
        test_file = tmp_path / "test.txt"
        test_file.write_text("test")

        with FileLock(str(test_file)) as lock:
            assert lock.is_locked

        assert not lock.is_locked

    def test_lock_automatic_cleanup(self, tmp_path):
        """Test that lock is cleaned up on deletion."""
        test_file = tmp_path / "test.txt"
        test_file.write_text("test")

        lock = FileLock(str(test_file))
        lock.acquire()
        lockfile = lock.lockfile

        assert Path(lockfile).exists()

        del lock
        assert not Path(lockfile).exists()


class TestSystem:
    """Test System class."""

    def test_system_default_name(self):
        """Test that system gets a random name by default."""
        system = System()

        assert system.name.startswith("mdmix_system_")
        assert len(system.name) == len("mdmix_system_") + 3

    def test_system_custom_name(self):
        """Test creating system with custom name."""
        system = System(name="myprotein")

        assert system.name == "myprotein"
        assert system.sys_file_path == "myprotein.msys"

    def test_system_with_ff(self):
        """Test system with force fields."""
        system = System(name="test", ff=["ff14SB", "gaff2"])

        assert "ff14SB" in system.ff
        assert "gaff2" in system.ff

    def test_system_with_extra_residues(self):
        """Test system with extra residue list."""
        system = System(name="test", extra_res_list=["LIG", "ION"])

        assert "LIG" in system.extra_res_list
        assert "ION" in system.extra_res_list

    def test_system_set_off(self, tmp_path):
        """Test setting OFF file."""
        off_file = tmp_path / "test.off"
        off_file.write_text(SAMPLE_OFF)

        system = System(name="test")
        system.set_off(off_file)

        assert system.off_manager is not None
        assert system.unit_name == "TST"

    def test_system_set_off_with_unit_name(self, tmp_path):
        """Test setting OFF file with explicit unit name."""
        off_file = tmp_path / "test.off"
        off_file.write_text(SAMPLE_OFF)

        system = System(name="test")
        system.set_off(off_file, unit_name="TST")

        assert system.unit_name == "TST"

    def test_system_from_off(self, tmp_path):
        """Test creating system from OFF file."""
        off_file = tmp_path / "test.off"
        off_file.write_text(SAMPLE_OFF)

        system = System.from_off(off_file, name="mytest")

        assert system.name == "mytest"
        assert system.off_manager is not None
        assert system.unit_name == "TST"

    def test_system_write_off(self, tmp_path):
        """Test writing OFF content to file."""
        off_file = tmp_path / "test.off"
        off_file.write_text(SAMPLE_OFF)

        system = System.from_off(off_file, name="test")
        output = tmp_path / "output.off"
        result = system.write_off(output)

        assert result == output.absolute()
        assert output.read_text() == SAMPLE_OFF

    def test_system_write_off_no_off(self):
        """Test write_off returns None when no OFF loaded."""
        system = System(name="test")

        result = system.write_off("output.off")

        assert result is None

    def test_system_copy(self, tmp_path):
        """Test copying a system."""
        off_file = tmp_path / "test.off"
        off_file.write_text(SAMPLE_OFF)

        original = System.from_off(off_file, name="original", ff=["ff14SB"])
        copy = original.copy()

        assert copy.name == original.name
        assert copy.ff == original.ff
        assert copy is not original

        # Modify copy, original should be unchanged
        copy.name = "modified"
        assert original.name == "original"

    def test_system_write_and_load(self, tmp_path):
        """Test saving and loading a system."""
        off_file = tmp_path / "test.off"
        off_file.write_text(SAMPLE_OFF)

        original = System.from_off(off_file, name="savetest", ff=["ff14SB"])
        save_file = tmp_path / "test.msys"
        original.write(save_file)

        assert save_file.exists()

        loaded = System.from_file(save_file)

        assert loaded.name == "savetest"
        assert loaded.ff == ["ff14SB"]
        assert loaded.unit_name == "TST"

    def test_system_load_not_found(self):
        """Test loading nonexistent file."""
        system = System(name="test")

        with pytest.raises(BadFileError, match="not found"):
            system.load("/nonexistent/path.msys")

    def test_system_solvate_no_off(self):
        """Test solvation fails without OFF file."""
        system = System(name="test")

        with pytest.raises(SystemError, match="Cannot solvate"):
            system.solvate("WAT")

    def test_system_repr(self):
        """Test string representation."""
        system = System(name="myprotein")

        assert "myprotein" in repr(system)
        assert "System" in repr(system)

    def test_system_str(self):
        """Test string output."""
        system = System(name="myprotein", ff=["ff14SB"], extra_res_list=["LIG"])

        s = str(system)

        assert "SYSTEM: myprotein" in s
        assert "ff14SB" in s
        assert "LIG" in s


class TestSolvatedSystem:
    """Test SolvatedSystem class."""

    def test_solvated_system_creation(self):
        """Test creating a solvated system."""
        system = SolvatedSystem(name="solvated", solvent="ETA")

        assert system.name == "solvated"
        assert system.solvent == "ETA"

    def test_solvated_system_from_top_crd(self, tmp_path):
        """Test creating from topology and coordinates."""
        top_file = tmp_path / "test.prmtop"
        crd_file = tmp_path / "test.prmcrd"
        top_file.write_text(SAMPLE_PRMTOP)
        crd_file.write_text(SAMPLE_PRMCRD)

        system = SolvatedSystem.from_top_crd(
            name="test",
            top=top_file,
            crd=crd_file,
            solvent="WAT",
        )

        assert system.name == "test"
        assert system.top == SAMPLE_PRMTOP
        assert system.crd == SAMPLE_PRMCRD
        assert system.solvent == "WAT"

    def test_set_top_crd(self, tmp_path):
        """Test setting topology and coordinate files."""
        top_file = tmp_path / "test.prmtop"
        crd_file = tmp_path / "test.prmcrd"
        top_file.write_text(SAMPLE_PRMTOP)
        crd_file.write_text(SAMPLE_PRMCRD)

        system = SolvatedSystem(name="test")
        system.set_top_crd(top_file, crd_file)

        assert system.top == SAMPLE_PRMTOP
        assert system.crd == SAMPLE_PRMCRD

    def test_set_top_crd_file_not_found(self, tmp_path):
        """Test error when files not found."""
        system = SolvatedSystem(name="test")

        with pytest.raises(BadFileError):
            system.set_top_crd(
                tmp_path / "nonexistent.prmtop",
                tmp_path / "test.prmcrd",
            )

    def test_save_top_crd(self, tmp_path):
        """Test saving topology and coordinates."""
        top_file = tmp_path / "test.prmtop"
        crd_file = tmp_path / "test.prmcrd"
        top_file.write_text(SAMPLE_PRMTOP)
        crd_file.write_text(SAMPLE_PRMCRD)

        system = SolvatedSystem.from_top_crd(
            name="test",
            top=top_file,
            crd=crd_file,
        )

        output_prefix = tmp_path / "output"
        result = system.save_top_crd(output_prefix)

        assert result
        assert (tmp_path / "output.prmtop").exists()
        assert (tmp_path / "output.prmcrd").exists()

    def test_save_top_crd_empty(self):
        """Test save fails when no topology loaded."""
        system = SolvatedSystem(name="test")

        result = system.save_top_crd("output")

        assert not result

    def test_save_pdb(self, tmp_path):
        """Test saving PDB content."""
        system = SolvatedSystem(name="test")
        system.pdb = "ATOM      1  C   TST     1       0.000   0.000   0.000"

        output = tmp_path / "output.pdb"
        result = system.save_pdb(output)

        assert result
        assert output.exists()
        assert "ATOM" in output.read_text()

    def test_get_tmp_top_crd_files(self, tmp_path):
        """Test getting temporary topology/coordinate files."""
        top_file = tmp_path / "test.prmtop"
        crd_file = tmp_path / "test.prmcrd"
        top_file.write_text(SAMPLE_PRMTOP)
        crd_file.write_text(SAMPLE_PRMCRD)

        system = SolvatedSystem.from_top_crd(
            name="test",
            top=top_file,
            crd=crd_file,
        )

        tmp_top, tmp_crd = system.get_tmp_top_crd_files()

        assert tmp_top.exists()
        assert tmp_crd.exists()
        assert tmp_top.read_text() == SAMPLE_PRMTOP

        system.clean_tmp()

        assert not tmp_top.exists()
        assert not tmp_crd.exists()

    def test_get_tmp_pdb_file(self, tmp_path):
        """Test getting temporary PDB file."""
        system = SolvatedSystem(name="test")
        system.pdb = "ATOM      1  C   TST     1       0.000   0.000   0.000"

        tmp_pdb = system.get_tmp_pdb_file()

        assert tmp_pdb.exists()
        assert "ATOM" in tmp_pdb.read_text()

        system.clean_tmp()

        assert not tmp_pdb.exists()

    def test_solvate_raises_error(self):
        """Test that solvate raises error on already solvated system."""
        system = SolvatedSystem(name="test", solvent="WAT")

        with pytest.raises(SystemError, match="Cannot solvate"):
            system.solvate(solvent="ETA")

    def test_repr(self):
        """Test string representation."""
        system = SolvatedSystem(name="test", solvent="ETA")

        repr_str = repr(system)

        assert "SolvatedSystem" in repr_str
        assert "test" in repr_str
        assert "ETA" in repr_str


class TestLoadSystem:
    """Test load_system function."""

    def test_load_system_from_path(self, tmp_path):
        """Test loading system from specific path."""
        # Create a saved system
        off_file = tmp_path / "test.off"
        off_file.write_text(SAMPLE_OFF)

        original = System.from_off(off_file, name="loadtest")
        save_file = tmp_path / "test.msys"
        original.write(save_file)

        loaded = load_system(save_file)

        assert loaded.name == "loadtest"

    def test_load_system_auto_find(self, tmp_path, monkeypatch):
        """Test auto-finding .msys file."""
        monkeypatch.chdir(tmp_path)

        off_file = tmp_path / "test.off"
        off_file.write_text(SAMPLE_OFF)

        original = System.from_off(off_file, name="autoload")
        save_file = tmp_path / "autoload.msys"
        original.write(save_file)

        loaded = load_system()

        assert loaded.name == "autoload"

    def test_load_system_no_file(self, tmp_path, monkeypatch):
        """Test error when no .msys file found."""
        monkeypatch.chdir(tmp_path)

        with pytest.raises(SystemError, match="No .msys file found"):
            load_system()

    def test_load_system_multiple_files(self, tmp_path, monkeypatch):
        """Test error when multiple .msys files found."""
        monkeypatch.chdir(tmp_path)

        (tmp_path / "file1.msys").write_text("")
        (tmp_path / "file2.msys").write_text("")

        with pytest.raises(SystemError, match="Multiple system files"):
            load_system()
