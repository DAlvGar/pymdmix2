"""
Tests for pymdmix.io.parsers and pymdmix.utils.settings_parser modules.
"""

import os
import tempfile

import pytest

from pymdmix.io.parsers import (
    BadFile,
    BadSolvent,
    MDSettingsConfigFileParser,
    MDSettingsParserError,
    SystemConfigFileParser,
    SystemParserError,
)
from pymdmix.utils.settings_parser import (
    CaseSensitiveConfigParser,
    InvalidValue,
    Setting,
    SettingsManager,
    SettingsParser,
)

# =============================================================================
# SettingsParser Tests
# =============================================================================


class TestSetting:
    """Tests for Setting class."""

    def test_setting_creation(self):
        """Test basic Setting creation."""
        s = Setting(name="test", value="42", vtype=str)
        assert s.name == "test"
        assert s.value == "42"
        assert s.vtype == str

    def test_setting_type_cast_int(self):
        """Test type casting to int."""
        s = Setting(name="count", value="42")
        s.type_cast(int)
        assert s.value == 42
        assert s.vtype == int

    def test_setting_type_cast_float(self):
        """Test type casting to float."""
        s = Setting(name="rate", value="3.14")
        s.type_cast(float)
        assert s.value == pytest.approx(3.14)
        assert s.vtype == float

    def test_setting_type_cast_bool_true(self):
        """Test type casting to bool (true)."""
        for val in ["true", "True", "yes", "1", "on"]:
            s = Setting(name="flag", value=val)
            s.type_cast(bool)
            assert s.value is True

    def test_setting_type_cast_bool_false(self):
        """Test type casting to bool (false)."""
        for val in ["false", "False", "no", "0", "off"]:
            s = Setting(name="flag", value=val)
            s.type_cast(bool)
            assert s.value is False

    def test_setting_type_cast_list(self):
        """Test type casting to list."""
        s = Setting(name="items", value="one, two, three")
        s.type_cast(list)
        assert s.value == ["one", "two", "three"]
        assert s.vtype == list

    def test_setting_type_cast_invalid(self):
        """Test type cast with invalid value."""
        s = Setting(name="num", value="not_a_number")
        with pytest.raises(InvalidValue):
            s.type_cast(int)

    def test_setting_formatted(self):
        """Test formatted output."""
        s = Setting(name="count", value=42, vtype=int)
        result = s.formatted()
        assert "int-count" in result
        assert "42" in result

    def test_setting_formatted_with_comment(self):
        """Test formatted output with comment."""
        s = Setting(name="test", value="val", comment="A comment")
        result = s.formatted()
        assert "## A comment" in result

    def test_setting_comparison(self):
        """Test Setting comparison (sorting)."""
        s1 = Setting(name="alpha")
        s2 = Setting(name="beta")
        assert s1 < s2


class TestCaseSensitiveConfigParser:
    """Tests for CaseSensitiveConfigParser."""

    def test_preserves_case(self):
        """Test that option names preserve case."""
        cfg_content = """
[TEST]
CamelCase = value
lowercase = value
UPPERCASE = value
"""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".cfg", delete=False) as f:
            f.write(cfg_content)
            path = f.name

        try:
            parser = CaseSensitiveConfigParser()
            parser.read(path)
            options = dict(parser.items("TEST"))
            assert "CamelCase" in options
            assert "lowercase" in options
            assert "UPPERCASE" in options
        finally:
            os.unlink(path)


class TestSettingsParser:
    """Tests for SettingsParser."""

    def test_parse_simple(self):
        """Test parsing simple settings."""
        cfg_content = """
[GENERAL]
name = test_name
"""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".cfg", delete=False) as f:
            f.write(cfg_content)
            path = f.name

        try:
            parser = SettingsParser(path)
            result = parser.parse()
            assert "name" in result
            assert result["name"].value == "test_name"
        finally:
            os.unlink(path)

    def test_parse_with_types(self):
        """Test parsing with type prefixes."""
        cfg_content = """
[GENERAL]
int-count = 42
float-rate = 3.14
bool-enabled = true
list-items = one, two, three
"""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".cfg", delete=False) as f:
            f.write(cfg_content)
            path = f.name

        try:
            parser = SettingsParser(path)
            result = parser.parse()
            assert result["count"].value == 42
            assert result["rate"].value == pytest.approx(3.14)
            assert result["enabled"].value is True
            assert result["items"].value == ["one", "two", "three"]
        finally:
            os.unlink(path)

    def test_parse_with_comments(self):
        """Test parsing with comments."""
        cfg_content = """
[GENERAL]
name = test # This is a comment
"""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".cfg", delete=False) as f:
            f.write(cfg_content)
            path = f.name

        try:
            parser = SettingsParser(path)
            result = parser.parse()
            assert result["name"].value == "test"
            assert "comment" in result["name"].comment.lower()
        finally:
            os.unlink(path)

    def test_parse_keep_sections(self):
        """Test parsing with sections preserved."""
        cfg_content = """
[GENERAL]
name = test

[OTHER]
value = 42
"""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".cfg", delete=False) as f:
            f.write(cfg_content)
            path = f.name

        try:
            parser = SettingsParser(path)
            result = parser.parse(keep_sections=True)
            assert "GENERAL" in result
            assert "OTHER" in result
            assert "name" in result["GENERAL"]
            assert "value" in result["OTHER"]
        finally:
            os.unlink(path)

    def test_parse_file_not_found(self):
        """Test parsing non-existent file."""
        parser = SettingsParser("/nonexistent/file.cfg")
        with pytest.raises(IOError):
            parser.parse()


class TestSettingsManager:
    """Tests for SettingsManager."""

    def test_settings_to_dict(self):
        """Test settings_to_dict conversion."""
        cfg_content = """
[GENERAL]
int-count = 42
name = test
"""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".cfg", delete=False) as f:
            f.write(cfg_content)
            default_path = f.name

        with tempfile.NamedTemporaryFile(mode="w", suffix=".cfg", delete=False) as f:
            f.write("")  # Empty user config
            user_path = f.name

        try:
            manager = SettingsManager(default_path, user_path, verbose=False)
            manager.collect_settings()
            d = manager.settings_to_dict()
            assert d["count"] == 42
            assert d["name"] == "test"
        finally:
            os.unlink(default_path)
            os.unlink(user_path)

    def test_user_overrides_default(self):
        """Test that user settings override defaults."""
        default_cfg = """
[GENERAL]
int-count = 10
"""
        user_cfg = """
[GENERAL]
int-count = 99
"""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".cfg", delete=False) as f:
            f.write(default_cfg)
            default_path = f.name

        with tempfile.NamedTemporaryFile(mode="w", suffix=".cfg", delete=False) as f:
            f.write(user_cfg)
            user_path = f.name

        try:
            manager = SettingsManager(default_path, user_path, verbose=False)
            manager.collect_settings()
            d = manager.settings_to_dict()
            assert d["count"] == 99
        finally:
            os.unlink(default_path)
            os.unlink(user_path)


# =============================================================================
# SystemConfigFileParser Tests
# =============================================================================


class TestSystemConfigFileParser:
    """Tests for SystemConfigFileParser."""

    def test_parse_basic(self):
        """Test parsing basic system config."""
        cfg_content = """
[SYSTEM]
name = my_system
restrmask = :1-100@CA
alignmask = auto
"""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".cfg", delete=False) as f:
            f.write(cfg_content)
            path = f.name

        try:
            parser = SystemConfigFileParser()
            params = parser.parse(path)
            assert params["name"] == "my_system"
            assert params["restrain_mask"] == ":1-100@CA"
            assert params["align_mask"] == "auto"
        finally:
            os.unlink(path)

    def test_parse_with_extra_residues(self):
        """Test parsing with extra residues."""
        cfg_content = """
[SYSTEM]
name = system
extrares = LIG, ION, ACE
"""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".cfg", delete=False) as f:
            f.write(cfg_content)
            path = f.name

        try:
            parser = SystemConfigFileParser()
            params = parser.parse(path)
            assert params["extra_residues"] == ["LIG", "ION", "ACE"]
        finally:
            os.unlink(path)

    def test_parse_missing_system_section(self):
        """Test parsing without SYSTEM section."""
        cfg_content = """
[OTHER]
name = test
"""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".cfg", delete=False) as f:
            f.write(cfg_content)
            path = f.name

        try:
            parser = SystemConfigFileParser()
            with pytest.raises(SystemParserError):
                parser.parse(path)
        finally:
            os.unlink(path)

    def test_parse_file_not_found(self):
        """Test parsing non-existent file."""
        parser = SystemConfigFileParser()
        with pytest.raises(BadFile):
            parser.parse("/nonexistent/file.cfg")


# =============================================================================
# MDSettingsConfigFileParser Tests
# =============================================================================


class TestMDSettingsConfigFileParser:
    """Tests for MDSettingsConfigFileParser."""

    def test_parse_basic(self):
        """Test parsing basic MD settings."""
        cfg_content = """
[MDSETTINGS]
solvents = ETA, WAT
nrepl = 2
nanos = 10
restr = HA
temp = 300.0
"""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".cfg", delete=False) as f:
            f.write(cfg_content)
            path = f.name

        try:
            parser = MDSettingsConfigFileParser()
            settings = parser.parse(path)
            # 2 solvents x 2 replicas = 4 settings
            assert len(settings) == 4
            assert all(s["nanos"] == 10 for s in settings)
            assert all(s["restr_mode"] == "HA" for s in settings)
            assert all(s["temp"] == 300.0 for s in settings)
        finally:
            os.unlink(path)

    def test_parse_per_solvent_nrepl(self):
        """Test parsing with per-solvent replica counts."""
        cfg_content = """
[MDSETTINGS]
solvents = ETA, WAT
nrepl = 3, WAT:1
"""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".cfg", delete=False) as f:
            f.write(cfg_content)
            path = f.name

        try:
            parser = MDSettingsConfigFileParser()
            settings = parser.parse(path)
            # ETA: 3 replicas, WAT: 1 replica = 4 total
            assert len(settings) == 4
            eta_count = sum(1 for s in settings if s["solvent"] == "ETA")
            wat_count = sum(1 for s in settings if s["solvent"] == "WAT")
            assert eta_count == 3
            assert wat_count == 1
        finally:
            os.unlink(path)

    def test_parse_per_replica_settings(self):
        """Test parsing with per-replica override."""
        cfg_content = """
[MDSETTINGS]
solvents = ETA
nrepl = 2
restr = HA, ETA/1/FREE
"""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".cfg", delete=False) as f:
            f.write(cfg_content)
            path = f.name

        try:
            parser = MDSettingsConfigFileParser()
            settings = parser.parse(path)
            assert len(settings) == 2
            # First replica should be FREE
            restr_modes = [s["restr_mode"] for s in settings]
            assert "FREE" in restr_modes
            assert "HA" in restr_modes
        finally:
            os.unlink(path)

    def test_parse_no_solvent(self):
        """Test parse_no_solvent method."""
        cfg_content = """
[MDSETTINGS]
nanos = 20
restr = BB
temp = 310.0
force = 5.0
"""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".cfg", delete=False) as f:
            f.write(cfg_content)
            path = f.name

        try:
            parser = MDSettingsConfigFileParser()
            settings = parser.parse_no_solvent(path)
            assert settings["nanos"] == 20
            assert settings["restr_mode"] == "BB"
            assert settings["temp"] == 310.0
            assert settings["restr_force"] == 5.0
        finally:
            os.unlink(path)

    def test_parse_missing_mdsettings_section(self):
        """Test parsing without MDSETTINGS section."""
        cfg_content = """
[OTHER]
name = test
"""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".cfg", delete=False) as f:
            f.write(cfg_content)
            path = f.name

        try:
            parser = MDSettingsConfigFileParser()
            with pytest.raises(MDSettingsParserError):
                parser.parse(path)
        finally:
            os.unlink(path)

    def test_parse_missing_solvents(self):
        """Test parsing without solvents."""
        cfg_content = """
[MDSETTINGS]
nrepl = 2
"""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".cfg", delete=False) as f:
            f.write(cfg_content)
            path = f.name

        try:
            parser = MDSettingsConfigFileParser()
            with pytest.raises(MDSettingsParserError):
                parser.parse(path)
        finally:
            os.unlink(path)

    def test_parse_invalid_solvent(self):
        """Test parsing with invalid solvent."""
        cfg_content = """
[MDSETTINGS]
solvents = INVALID
nrepl = 1
"""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".cfg", delete=False) as f:
            f.write(cfg_content)
            path = f.name

        try:
            # With validation enabled
            parser = MDSettingsConfigFileParser(available_solvents=["ETA", "WAT"])
            with pytest.raises(BadSolvent):
                parser.parse(path)
        finally:
            os.unlink(path)

    def test_parse_multiple_sections(self):
        """Test parsing multiple MDSETTINGS sections."""
        cfg_content = """
[MDSETTINGS1]
solvents = ETA
nrepl = 1
nanos = 10

[MDSETTINGS2]
solvents = WAT
nrepl = 2
nanos = 20
"""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".cfg", delete=False) as f:
            f.write(cfg_content)
            path = f.name

        try:
            parser = MDSettingsConfigFileParser()
            settings = parser.parse(path)
            # ETA: 1 + WAT: 2 = 3 total
            assert len(settings) == 3
        finally:
            os.unlink(path)
