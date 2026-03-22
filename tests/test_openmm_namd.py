"""Tests for OpenMM and NAMD engine interfaces."""

import pytest

from pymdmix.engines.namd import NAMDConfig, NAMDEngine
from pymdmix.engines.openmm import OpenMMConfig, OpenMMEngine


class TestOpenMMConfig:
    """Tests for OpenMM configuration."""

    def test_default_config(self):
        """Test default configuration values."""
        config = OpenMMConfig()

        assert config.platform == "CUDA"
        assert config.precision == "mixed"
        assert config.temperature == 300.0
        assert config.timestep == 2.0

    def test_custom_config(self):
        """Test custom configuration."""
        config = OpenMMConfig(
            platform="CPU",
            temperature=310.0,
            pressure=2.0,
        )

        assert config.platform == "CPU"
        assert config.temperature == 310.0
        assert config.pressure == 2.0


class TestOpenMMEngine:
    """Tests for OpenMM engine."""

    def test_write_minimization(self, tmp_path):
        """Test writing minimization script."""
        engine = OpenMMEngine()

        output = engine.write_minimization(
            topology="system.prmtop",
            coordinates="system.inpcrd",
            output=tmp_path / "min.py",
            steps=1000,
        )

        assert output.exists()
        content = output.read_text()

        assert "openmm" in content.lower() or "simtk" in content.lower()
        assert "minimizeEnergy" in content
        assert "1000" in content  # steps
        assert "system.prmtop" in content

    def test_write_equilibration(self, tmp_path):
        """Test writing equilibration script."""
        engine = OpenMMEngine()

        output = engine.write_equilibration(
            topology="system.prmtop",
            coordinates="min.rst7",
            output=tmp_path / "eq.py",
            steps=50000,
            ensemble="NPT",
        )

        assert output.exists()
        content = output.read_text()

        assert "LangevinMiddleIntegrator" in content
        assert "MonteCarloBarostat" in content  # NPT
        assert "50000" in content

    def test_write_production(self, tmp_path):
        """Test writing production script."""
        engine = OpenMMEngine()

        output = engine.write_production(
            topology="system.prmtop",
            coordinates="eq.rst7",
            output=tmp_path / "prod.py",
            steps=100000,
            ensemble="NVT",
        )

        assert output.exists()
        content = output.read_text()

        assert "DCDReporter" in content
        assert "100000" in content
        assert "MonteCarloBarostat" not in content  # NVT

    def test_restraint_code(self, tmp_path):
        """Test restraint code generation."""
        engine = OpenMMEngine()

        output = engine.write_minimization(
            topology="system.prmtop",
            coordinates="system.inpcrd",
            output=tmp_path / "min.py",
            restraint_mask="backbone",
        )

        content = output.read_text()
        assert "CustomExternalForce" in content
        assert "restraint" in content.lower()

    def test_script_is_executable(self, tmp_path):
        """Test that generated scripts are executable."""
        engine = OpenMMEngine()

        output = engine.write_minimization(
            topology="system.prmtop",
            coordinates="system.inpcrd",
            output=tmp_path / "min.py",
        )

        # Check executable permission
        import stat

        mode = output.stat().st_mode
        assert mode & stat.S_IXUSR  # User executable


class TestNAMDConfig:
    """Tests for NAMD configuration."""

    def test_default_config(self):
        """Test default configuration values."""
        config = NAMDConfig()

        assert config.cutoff == 9.0
        assert config.temperature == 300.0
        assert config.timestep == 2.0

    def test_custom_config(self):
        """Test custom configuration."""
        config = NAMDConfig(
            cutoff=12.0,
            temperature=310.0,
            langevin_damping=5.0,
        )

        assert config.cutoff == 12.0
        assert config.temperature == 310.0
        assert config.langevin_damping == 5.0


class TestNAMDEngine:
    """Tests for NAMD engine."""

    @pytest.fixture
    def mock_inpcrd(self, tmp_path):
        """Create a mock inpcrd file with box dimensions."""
        inpcrd = tmp_path / "system.inpcrd"
        inpcrd.write_text("""title
12345
   1.0   2.0   3.0
  50.000  50.000  50.000  90.000  90.000  90.000
""")
        return inpcrd

    def test_write_minimization(self, tmp_path, mock_inpcrd):
        """Test writing minimization config."""
        engine = NAMDEngine()

        output = engine.write_minimization(
            topology=tmp_path / "system.prmtop",
            coordinates=mock_inpcrd,
            output=tmp_path / "min.namd",
            steps=5000,
        )

        assert output.exists()
        content = output.read_text()

        assert "amber" in content
        assert "minimization" in content.lower()
        assert "5000" in content
        assert "PME" in content

    def test_write_heating(self, tmp_path, mock_inpcrd):
        """Test writing heating config."""
        engine = NAMDEngine()

        output = engine.write_heating(
            topology=tmp_path / "system.prmtop",
            coordinates=mock_inpcrd,
            output=tmp_path / "heat.namd",
            start_temp=0.0,
            target_temp=300.0,
        )

        assert output.exists()
        content = output.read_text()

        assert "langevin" in content.lower()
        assert "reassignFreq" in content
        assert "reassignIncr" in content

    def test_write_equilibration_npt(self, tmp_path, mock_inpcrd):
        """Test writing NPT equilibration config."""
        engine = NAMDEngine()

        output = engine.write_equilibration(
            topology=tmp_path / "system.prmtop",
            coordinates=mock_inpcrd,
            output=tmp_path / "eq.namd",
            ensemble="NPT",
        )

        assert output.exists()
        content = output.read_text()

        assert "langevinPiston" in content  # NPT barostat
        assert "langevinPistonTarget" in content

    def test_write_equilibration_nvt(self, tmp_path, mock_inpcrd):
        """Test writing NVT equilibration config."""
        engine = NAMDEngine()

        output = engine.write_equilibration(
            topology=tmp_path / "system.prmtop",
            coordinates=mock_inpcrd,
            output=tmp_path / "eq.namd",
            ensemble="NVT",
        )

        assert output.exists()
        content = output.read_text()

        assert "langevin" in content.lower()
        assert "langevinPiston" not in content  # No barostat for NVT

    def test_write_production(self, tmp_path, mock_inpcrd):
        """Test writing production config."""
        engine = NAMDEngine()

        output = engine.write_production(
            topology=tmp_path / "system.prmtop",
            coordinates=mock_inpcrd,
            output=tmp_path / "prod.namd",
            steps=500000,
            dcd_freq=500,
        )

        assert output.exists()
        content = output.read_text()

        assert "500000" in content
        assert "DCDfreq" in content
        assert "500" in content

    def test_box_dimensions_parsed(self, tmp_path, mock_inpcrd):
        """Test that box dimensions are correctly parsed."""
        engine = NAMDEngine()

        output = engine.write_minimization(
            topology=tmp_path / "system.prmtop",
            coordinates=mock_inpcrd,
            output=tmp_path / "min.namd",
        )

        content = output.read_text()

        # Should have cell basis vectors from the 50x50x50 box
        assert "cellBasisVector1" in content
        assert "50.000" in content

    def test_restraint_section(self, tmp_path, mock_inpcrd):
        """Test restraint configuration."""
        engine = NAMDEngine()

        output = engine.write_minimization(
            topology=tmp_path / "system.prmtop",
            coordinates=mock_inpcrd,
            output=tmp_path / "min.namd",
            restraint_pdb="restraints.pdb",
        )

        content = output.read_text()

        assert "constraints" in content.lower()
        assert "consref" in content
        assert "restraints.pdb" in content


class TestEngineIntegration:
    """Integration tests for engines."""

    def test_openmm_config_in_script(self, tmp_path):
        """Test that config values appear in generated script."""
        config = OpenMMConfig(
            platform="OpenCL",
            precision="double",
            temperature=310.0,
            nonbonded_cutoff=12.0,
        )
        engine = OpenMMEngine(config)

        output = engine.write_production(
            topology="test.prmtop",
            coordinates="test.rst7",
            output=tmp_path / "prod.py",
        )

        content = output.read_text()

        assert "OpenCL" in content
        assert "double" in content
        assert "310.0" in content
        assert "12.0" in content

    def test_namd_config_in_script(self, tmp_path):
        """Test that config values appear in generated config."""
        # Create mock inpcrd
        inpcrd = tmp_path / "test.inpcrd"
        inpcrd.write_text("title\n0\n  50.0  50.0  50.0\n")

        config = NAMDConfig(
            cutoff=12.0,
            temperature=310.0,
            timestep=1.0,
        )
        engine = NAMDEngine(config)

        output = engine.write_production(
            topology=tmp_path / "test.prmtop",
            coordinates=inpcrd,
            output=tmp_path / "prod.namd",
        )

        content = output.read_text()

        assert "cutoff                  12.0" in content
        assert "langevinTemp            310.0" in content
        assert "timestep                1.0" in content
