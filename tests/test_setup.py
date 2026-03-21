"""
Tests for pymdmix.setup module.
"""

import pytest
import numpy as np
from pathlib import Path

from pymdmix.setup.prepare import (
    prepare_structure,
    add_caps,
    find_and_fix_disulfides,
    renumber_residues,
    PrepareResult,
)
from pymdmix.setup.solvate import (
    BoxConfig,
    IonConfig,
    generate_leap_script,
    SolvateResult,
)
from pymdmix.core.structure import load_structure
from pymdmix.core.solvent import Solvent


class TestPrepareResult:
    """Test PrepareResult dataclass."""
    
    def test_default_result(self):
        """Test default prepare result."""
        import parmed
        struct = parmed.Structure()
        result = PrepareResult(structure=struct)
        
        assert result.n_caps_added == 0
        assert result.n_disulfides == 0
        assert result.warnings == []


class TestPrepareStructure:
    """Test structure preparation."""
    
    def test_prepare_simple_structure(self, sample_pdb_file):
        """Test preparing a simple structure."""
        struct = load_structure(sample_pdb_file)
        result = prepare_structure(struct)
        
        assert isinstance(result, PrepareResult)
        assert result.structure is not None
    
    def test_prepare_remove_waters(self, sample_pdb_file):
        """Test removing waters during preparation."""
        struct = load_structure(sample_pdb_file)
        
        # Sample has 3 water atoms
        result = prepare_structure(struct, remove_waters=True)
        
        assert result.n_waters_removed == 3
    
    def test_prepare_with_options(self, sample_pdb_file):
        """Test preparation with various options."""
        struct = load_structure(sample_pdb_file)
        
        result = prepare_structure(
            struct,
            cap_termini=True,
            fix_disulfides=True,
            remove_waters=False,
        )
        
        assert isinstance(result, PrepareResult)


class TestDisulfides:
    """Test disulfide handling."""
    
    def test_find_no_disulfides(self, sample_pdb_file):
        """Test finding disulfides when none present."""
        struct = load_structure(sample_pdb_file)
        disulfides = find_and_fix_disulfides(struct)
        
        # Sample PDB has no CYS
        assert disulfides == []


class TestRenumberResidues:
    """Test residue renumbering."""
    
    def test_renumber_from_one(self, sample_pdb_file):
        """Test renumbering residues from 1."""
        struct = load_structure(sample_pdb_file)
        renumber_residues(struct, start=1)
        
        assert struct.residues[0].number == 1
        assert struct.residues[1].number == 2
    
    def test_renumber_from_custom(self, sample_pdb_file):
        """Test renumbering from custom start."""
        struct = load_structure(sample_pdb_file)
        renumber_residues(struct, start=100)
        
        assert struct.residues[0].number == 100


class TestBoxConfig:
    """Test BoxConfig dataclass."""
    
    def test_default_config(self):
        """Test default box configuration."""
        config = BoxConfig()
        
        assert config.shape == "truncated_octahedron"
        assert config.buffer == 12.0
    
    def test_leap_command_oct(self):
        """Test LEaP command for truncated octahedron."""
        config = BoxConfig(shape="truncated_octahedron", buffer=10.0)
        cmd = config.to_leap_command("sys", "ETABOX")
        
        assert "solvateoct" in cmd
        assert "sys" in cmd
        assert "ETABOX" in cmd
        assert "10.0" in cmd
    
    def test_leap_command_box(self):
        """Test LEaP command for rectangular box."""
        config = BoxConfig(shape="box", buffer=15.0)
        cmd = config.to_leap_command("sys", "WATBOX")
        
        assert "solvatebox" in cmd


class TestIonConfig:
    """Test IonConfig dataclass."""
    
    def test_default_config(self):
        """Test default ion configuration."""
        config = IonConfig()
        
        assert config.neutralize is True
        assert config.cation == "Na+"
        assert config.anion == "Cl-"
    
    def test_leap_commands(self):
        """Test LEaP command generation."""
        config = IonConfig(neutralize=True)
        commands = config.to_leap_commands("sys")
        
        assert len(commands) >= 2
        assert any("Na+" in cmd for cmd in commands)


class TestLeapScript:
    """Test LEaP script generation."""
    
    def test_generate_basic_script(self, tmp_output_dir):
        """Test generating a basic LEaP script."""
        pdb_path = tmp_output_dir / "test.pdb"
        pdb_path.write_text("ATOM      1  CA  ALA A   1       0.0   0.0   0.0  1.0  0.0\nEND\n")
        
        solvent = Solvent(
            name="ETA",
            description="Ethanol",
        )
        
        script = generate_leap_script(
            pdb_path=pdb_path,
            solvent=solvent,
            output_prefix="system",
        )
        
        assert "source leaprc" in script
        assert "loadpdb" in script
        assert "saveamberparm" in script
        assert "quit" in script
    
    def test_script_with_custom_ff(self, tmp_output_dir):
        """Test script with custom force fields."""
        pdb_path = tmp_output_dir / "test.pdb"
        pdb_path.write_text("END\n")
        
        solvent = Solvent(name="WAT")
        
        script = generate_leap_script(
            pdb_path=pdb_path,
            solvent=solvent,
            output_prefix="system",
            force_fields=["leaprc.protein.ff14SB", "leaprc.water.tip3p"],
        )
        
        assert "ff14SB" in script
        assert "tip3p" in script


class TestSolvateResult:
    """Test SolvateResult dataclass."""
    
    def test_default_result(self):
        """Test default solvate result."""
        result = SolvateResult()
        
        assert result.success is False
        assert result.topology is None
        assert result.error is None
    
    def test_successful_result(self, tmp_output_dir):
        """Test successful solvation result."""
        result = SolvateResult(
            success=True,
            topology=tmp_output_dir / "system.prmtop",
            coordinates=tmp_output_dir / "system.inpcrd",
            n_solvent_residues=1000,
        )
        
        assert result.success is True
        assert result.n_solvent_residues == 1000
