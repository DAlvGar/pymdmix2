"""Tests for GROMACS engine interface."""

import pytest

from pymdmix.engines.gromacs import (
    GromacsCheck,
    GromacsCheckResult,
    GromacsConfig,
    GromacsEngine,
)


class TestGromacsConfig:
    """Tests for GROMACS configuration."""

    def test_default_config(self):
        """Test default configuration values."""
        config = GromacsConfig()

        assert config.timestep == 2.0
        assert config.temperature == 300.0
        assert config.pressure == 1.0
        assert config.cutoff == 0.9

    def test_custom_config(self):
        """Test custom configuration."""
        config = GromacsConfig(
            timestep=4.0,
            temperature=310.0,
            num_threads=8,
        )

        assert config.timestep == 4.0
        assert config.temperature == 310.0
        assert config.num_threads == 8


class TestGromacsEngine:
    """Tests for GROMACS engine."""

    @pytest.fixture
    def engine(self):
        """Create engine instance."""
        return GromacsEngine()

    @pytest.fixture
    def mock_gro(self, tmp_path):
        """Create a mock GRO file."""
        gro = tmp_path / "system.gro"
        gro.write_text("""Test system
    6
    1ALA      N    1   0.000   0.000   0.000
    1ALA     CA    2   0.100   0.000   0.000
    1ALA      C    3   0.200   0.000   0.000
    1ALA      O    4   0.300   0.000   0.000
    1WAT     OW    5   0.500   0.500   0.500
    1WAT    HW1    6   0.550   0.500   0.500
   5.00000   5.00000   5.00000
""")
        return gro

    def test_write_minimization(self, engine, tmp_path):
        """Test writing minimization MDP."""
        output = engine.write_minimization(tmp_path / "min.mdp", steps=1000)

        assert output.exists()
        content = output.read_text()

        assert "integrator" in content
        assert "steep" in content
        assert "1000" in content
        assert "emtol" in content

    def test_write_minimization_no_restraints(self, engine, tmp_path):
        """Test minimization without restraints."""
        output = engine.write_minimization(
            tmp_path / "min.mdp",
            restraints=False,
        )

        content = output.read_text()
        assert "POSRES" not in content

    def test_write_nvt_equilibration(self, engine, tmp_path):
        """Test writing NVT equilibration MDP."""
        output = engine.write_nvt_equilibration(
            tmp_path / "nvt.mdp",
            steps=50000,
        )

        assert output.exists()
        content = output.read_text()

        assert "V-rescale" in content
        assert "gen_vel" in content
        assert "pcoupl" in content
        assert "no" in content  # No pressure coupling

    def test_write_npt_equilibration(self, engine, tmp_path):
        """Test writing NPT equilibration MDP."""
        output = engine.write_npt_equilibration(
            tmp_path / "npt.mdp",
            steps=100000,
        )

        assert output.exists()
        content = output.read_text()

        assert "Berendsen" in content
        assert "ref_p" in content

    def test_write_production_nvt(self, engine, tmp_path):
        """Test writing NVT production MDP."""
        output = engine.write_production(
            tmp_path / "prod.mdp",
            steps=500000,
            ensemble="NVT",
        )

        assert output.exists()
        content = output.read_text()

        assert "500000" in content
        assert "pcoupl                  = no" in content

    def test_write_production_npt(self, engine, tmp_path):
        """Test writing NPT production MDP."""
        output = engine.write_production(
            tmp_path / "prod.mdp",
            steps=500000,
            ensemble="NPT",
        )

        assert output.exists()
        content = output.read_text()

        assert "Parrinello-Rahman" in content

    def test_write_restraints_itp(self, engine, mock_gro, tmp_path):
        """Test writing position restraints ITP."""
        output = engine.write_restraints_itp(
            mock_gro,
            tmp_path / "posre.itp",
            selection="protein",
        )

        assert output.exists()
        content = output.read_text()

        assert "position_restraints" in content
        # Should have protein atoms (not water)
        assert "1000" in content  # Force constant

    def test_write_restraints_backbone_only(self, engine, mock_gro, tmp_path):
        """Test restraints for backbone atoms only."""
        output = engine.write_restraints_itp(
            mock_gro,
            tmp_path / "posre.itp",
            selection="backbone",
        )

        content = output.read_text()
        lines = [line for line in content.split("\n") if line.strip() and not line.startswith(";")]

        # Should have fewer atoms than full protein
        assert len(lines) < 6

    def test_create_index_groups(self, engine, mock_gro, tmp_path):
        """Test creating index file."""
        output = engine.create_index_groups(
            mock_gro,
            tmp_path / "index.ndx",
        )

        assert output.exists()
        content = output.read_text()

        assert "System" in content
        assert "Protein" in content
        assert "Non-Protein" in content

    def test_config_in_mdp(self, tmp_path):
        """Test that config values appear in MDP files."""
        config = GromacsConfig(
            timestep=4.0,
            temperature=310.0,
            cutoff=1.2,
        )
        engine = GromacsEngine(config)

        output = engine.write_production(tmp_path / "prod.mdp")
        content = output.read_text()

        assert "0.0040" in content  # 4.0 fs = 0.004 ps
        assert "310.0" in content
        assert "1.2" in content


class TestGromacsCheck:
    """Tests for GROMACS simulation checking."""

    def test_check_result_complete(self):
        """Test check result completeness."""
        result = GromacsCheckResult(
            minimization=True,
            equilibration=True,
            production=True,
        )
        assert result.complete

    def test_check_result_incomplete(self):
        """Test incomplete check result."""
        result = GromacsCheckResult(
            minimization=True,
            equilibration=True,
            production=False,
        )
        assert not result.complete

    def test_check_minimization(self, tmp_path):
        """Test minimization check."""
        checker = GromacsCheck()

        # No log file
        assert not checker.check_minimization(tmp_path)

        # Create log file with success message
        log = tmp_path / "min2.log"
        log.write_text("Some output\nFinished mdrun\nMore output")

        assert checker.check_minimization(tmp_path)

    def test_check_equilibration(self, tmp_path):
        """Test equilibration check."""
        checker = GromacsCheck()

        # No log files
        assert not checker.check_equilibration(tmp_path)

        # Create both nvt and npt logs
        (tmp_path / "nvt.log").write_text("Finished mdrun")
        assert not checker.check_equilibration(tmp_path)  # Still missing npt

        (tmp_path / "npt.log").write_text("Finished mdrun")
        assert checker.check_equilibration(tmp_path)

    def test_check_production(self, tmp_path):
        """Test production check."""
        checker = GromacsCheck()

        assert not checker.check_production(tmp_path)

        (tmp_path / "prod.log").write_text("Finished mdrun")
        assert checker.check_production(tmp_path)

    def test_check_all(self, tmp_path):
        """Test checking all stages."""
        checker = GromacsCheck()

        result = checker.check_all(tmp_path)
        assert not result.complete

        # Create all success logs
        (tmp_path / "min2.log").write_text("Finished mdrun")
        (tmp_path / "nvt.log").write_text("Finished mdrun")
        (tmp_path / "npt.log").write_text("Finished mdrun")
        (tmp_path / "prod.log").write_text("Finished mdrun")

        result = checker.check_all(tmp_path)
        assert result.complete

    def test_get_box_volume(self, tmp_path):
        """Test extracting box volume from log."""
        checker = GromacsCheck()

        log = tmp_path / "prod.log"
        log.write_text("""
Some output
           Box-X          Box-Y          Box-Z
     8.00000e+00    8.00000e+00    8.00000e+00
More output
""")

        volume = checker.get_box_volume(log)

        # 8nm * 8nm * 8nm = 512 nm³ = 512000 Å³ * 0.77 ≈ 394240
        assert volume is not None
        assert 300000 < volume < 500000


class TestGromacsIntegration:
    """Integration tests for GROMACS workflow."""

    def test_full_mdp_set(self, tmp_path):
        """Test generating complete set of MDP files."""
        engine = GromacsEngine()

        files = [
            engine.write_minimization(tmp_path / "min1.mdp", algorithm="steep"),
            engine.write_minimization(tmp_path / "min2.mdp", algorithm="cg"),
            engine.write_nvt_equilibration(tmp_path / "nvt.mdp"),
            engine.write_npt_equilibration(tmp_path / "npt.mdp"),
            engine.write_production(tmp_path / "prod.mdp"),
        ]

        assert all(f.exists() for f in files)

        # Check each file has required parameters
        for f in files:
            content = f.read_text()
            assert "cutoff-scheme" in content
            assert "PME" in content
