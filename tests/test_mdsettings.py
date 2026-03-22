"""Tests for MDSettings class."""

import pytest
import tempfile
from pathlib import Path

from pymdmix.project.settings import MDSettings, load_settings


class TestMDSettings:
    """Tests for MDSettings dataclass."""

    def test_defaults(self):
        """Test default values."""
        settings = MDSettings(solvent="WAT")
        assert settings.solvent == "WAT"
        assert settings.nanos == 20
        assert settings.temperature == 300.0
        assert settings.timestep == 2.0
        assert settings.restraint_mode == "FREE"
        assert settings.md_program == "AMBER"

    def test_auto_name(self):
        """Test automatic name generation."""
        settings = MDSettings(solvent="ETA")
        assert settings.name == "Settings_ETA_FREE"

        settings = MDSettings(solvent="MAM", restraint_mode="HA")
        assert settings.name == "Settings_MAM_HA"

    def test_custom_name(self):
        """Test custom name."""
        settings = MDSettings(solvent="WAT", name="my_custom_settings")
        assert settings.name == "my_custom_settings"

    def test_has_restraints(self):
        """Test has_restraints property."""
        free = MDSettings(solvent="WAT", restraint_mode="FREE")
        assert not free.has_restraints

        restrained = MDSettings(solvent="WAT", restraint_mode="HA")
        assert restrained.has_restraints

    def test_n_trajectory_files(self):
        """Test trajectory file count calculation."""
        # 50ns, 2fs timestep, 500k steps per file
        # total steps = 50 * 1e6 / 2 = 25M steps
        # files = 25M / 500k = 50
        settings = MDSettings(
            solvent="WAT",
            nanos=50,
            timestep=2.0,
            production_steps=500000,
        )
        assert settings.n_trajectory_files == 50

    def test_n_snapshots(self):
        """Test snapshot count calculation."""
        # 20ns, 500k steps/file, 500 freq = 1000 snaps/ns * 20 = 20000
        settings = MDSettings(
            solvent="WAT",
            nanos=20,
            production_steps=500000,
            trajectory_frequency=500,
        )
        assert settings.n_snapshots == 20000

    def test_to_dict(self):
        """Test dictionary conversion."""
        settings = MDSettings(solvent="ETA", nanos=40)
        d = settings.to_dict()
        assert d["solvent"] == "ETA"
        assert d["nanos"] == 40
        assert isinstance(d, dict)

    def test_from_dict(self):
        """Test creation from dictionary."""
        data = {
            "solvent": "MAM",
            "nanos": 30,
            "temperature": 310.0,
            "restraint_mode": "BB",
        }
        settings = MDSettings.from_dict(data)
        assert settings.solvent == "MAM"
        assert settings.nanos == 30
        assert settings.temperature == 310.0
        assert settings.restraint_mode == "BB"

    def test_from_dict_ignores_unknown(self):
        """Test that unknown keys are ignored."""
        data = {
            "solvent": "WAT",
            "unknown_key": "should_be_ignored",
        }
        settings = MDSettings.from_dict(data)
        assert settings.solvent == "WAT"
        assert not hasattr(settings, "unknown_key")

    def test_invalid_restraint_mode(self):
        """Test validation of restraint mode."""
        with pytest.raises(ValueError):
            MDSettings(solvent="WAT", restraint_mode="INVALID")

    def test_summary(self):
        """Test summary output."""
        settings = MDSettings(solvent="ETA", nanos=40, restraint_mode="HA")
        summary = settings.summary()
        assert "ETA" in summary
        assert "40 ns" in summary
        assert "HA" in summary

    def test_repr(self):
        """Test string representation."""
        settings = MDSettings(solvent="WAT")
        assert "WAT" in repr(settings)


class TestLoadSettings:
    """Tests for load_settings function."""

    def test_load_with_solvent(self, tmp_path):
        """Test loading with solvent override."""
        toml_content = """
[mdsettings]
nanos = 30
temperature = 310.0
"""
        config = tmp_path / "settings.toml"
        config.write_text(toml_content)

        settings = load_settings(config, solvent="ETA")
        assert settings.solvent == "ETA"
        assert settings.nanos == 30
        assert settings.temperature == 310.0

    def test_load_with_overrides(self, tmp_path):
        """Test loading with additional overrides."""
        toml_content = """
[mdsettings]
solvent = "WAT"
nanos = 20
"""
        config = tmp_path / "settings.toml"
        config.write_text(toml_content)

        settings = load_settings(config, nanos=50)
        assert settings.solvent == "WAT"
        assert settings.nanos == 50  # Override applied

    def test_load_missing_solvent(self, tmp_path):
        """Test error when solvent is missing."""
        toml_content = """
[mdsettings]
nanos = 20
"""
        config = tmp_path / "settings.toml"
        config.write_text(toml_content)

        with pytest.raises(ValueError, match="Solvent name is required"):
            load_settings(config)


class TestMDSettingsEquality:
    """Tests for settings equality."""

    def test_equal_settings(self):
        """Test equal settings."""
        s1 = MDSettings(solvent="WAT", name="test")
        s2 = MDSettings(solvent="WAT", name="test")
        assert s1 == s2

    def test_unequal_settings(self):
        """Test unequal settings."""
        s1 = MDSettings(solvent="WAT")
        s2 = MDSettings(solvent="ETA")
        assert s1 != s2


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
