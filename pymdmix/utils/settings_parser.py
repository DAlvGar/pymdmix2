"""
Settings Parser
===============

Read and parse INI-style settings files with type conversion.

Used for project configuration, replica configuration, and general settings.

Migrated from the original pyMDMix/SettingsParser.py.
"""

from __future__ import annotations

import configparser
import logging
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from pymdmix.utils.tools import absfile, backup, valid_binary, valid_path

log = logging.getLogger(__name__)


# =============================================================================
# Exceptions
# =============================================================================


class SettingsError(Exception):
    """Base exception for settings parsing."""

    pass


class InvalidType(SettingsError):
    """Raised when a type conversion fails."""

    pass


class InvalidValue(SettingsError):
    """Raised when a value is invalid."""

    pass


class InvalidFile(SettingsError):
    """Raised when a settings file is invalid."""

    pass


class SettingsWarning(SettingsError):
    """Non-fatal settings issue."""

    pass


class InvalidPath(SettingsWarning):
    """Path validation warning."""

    pass


class InvalidBinary(SettingsWarning):
    """Binary validation warning."""

    pass


class WriteCfgError(SettingsError):
    """Error writing configuration file."""

    pass


# =============================================================================
# Case-Sensitive ConfigParser
# =============================================================================


class CaseSensitiveConfigParser(configparser.ConfigParser):
    """ConfigParser that preserves option name case."""

    def optionxform(self, optionstr: str) -> str:
        return optionstr


# =============================================================================
# Setting Container
# =============================================================================


@dataclass
class Setting:
    """
    Container for a single configuration parameter.

    Attributes
    ----------
    name : str
        Parameter name.
    value : Any
        Parameter value.
    vtype : type
        Value type (str, int, float, bool, list).
    comment : str
        Optional comment.
    error : str | None
        Error message if value is invalid.
    section : str
        Section this parameter belongs to.
    """

    # Section constants
    GENERAL: str = field(default="GENERAL", repr=False, init=False)
    BIN: str = field(default="BINARIES", repr=False, init=False)
    PATH: str = field(default="PATHS", repr=False, init=False)

    name: str = ""
    value: Any = None
    vtype: type = str
    comment: str = ""
    error: str | None = None
    section: str = "GENERAL"

    def type_cast(self, vtype: type) -> None:
        """
        Recast value to a new type.

        Parameters
        ----------
        vtype : type
            New type (str, int, float, bool, list).

        Raises
        ------
        InvalidValue
            If conversion fails.
        """
        try:
            if self.value is not None:
                if vtype is list and isinstance(self.value, str):
                    self.value = [s.strip() for s in self.value.split(",")]
                elif vtype is bool and isinstance(self.value, str):
                    # Handle bool specially
                    self.value = self.value.lower() in ("true", "yes", "1", "on")
                else:
                    self.value = vtype(self.value)
            self.vtype = vtype
        except ValueError as e:
            raise InvalidValue(
                f"{self.name}: cannot convert '{self.value}' to {vtype.__name__}"
            ) from e

    def formatted(self) -> str:
        """
        Format parameter for settings file.

        Returns
        -------
        str
            Formatted parameter line.
        """
        comment = ""
        error = ""
        name = self.name
        value = self.value if self.value is not None else ""

        # Prefix name with type if not str
        if self.vtype is not str:
            name = f"{self.vtype.__name__}-{name}"
            if self.vtype is list and isinstance(value, list):
                value = ",".join(str(v) for v in value)

        if self.comment:
            comment = f"\t## {self.comment}"

        if self.error:
            error = f"\t#! {self.error} !"

        return f"{name} = {value}{comment}{error}"

    def __repr__(self) -> str:
        error = "(!)" if self.error else ""
        return f"{error}{self.name} = {self.vtype.__name__}({self.value})"

    def __lt__(self, other: Setting) -> bool:
        return self.name < other.name if isinstance(other, Setting) else NotImplemented


# =============================================================================
# Settings Parser
# =============================================================================


class SettingsParser:
    """
    Parse INI-style settings files with type conversion.

    Performs:
    - Read ini-style settings file
    - Type-cast options (e.g., 'int-some_name' -> int)
    - Validate [PATHS] entries exist
    - Validate [BINARIES] entries are executable

    Parameters
    ----------
    ini_file : str | Path
        Path to configuration file.

    Examples
    --------
    >>> parser = SettingsParser('config.cfg')
    >>> settings = parser.parse()
    >>> settings['my_param'].value
    42
    """

    def __init__(self, ini_file: str | Path):
        self.log = logging.getLogger("SettingsParser")
        self.f_ini = absfile(ini_file)
        self.result: dict[str, Setting | dict[str, Setting]] = {}

    def _extract_type(self, option: str, default: type = str) -> tuple[type, str]:
        """
        Extract type prefix from option name.

        Parameters
        ----------
        option : str
            Option name (e.g., 'int-count', 'my_param').
        default : type
            Default type if no prefix.

        Returns
        -------
        tuple[type, str]
            (type, cleaned option name)

        Raises
        ------
        TypeError
            If type prefix is invalid.
        """
        t = default
        o = option

        if "-" in option:
            parts = option.split("-", 1)
            type_str = parts[0]
            o = parts[1]

            # Map string to type
            type_map = {
                "str": str,
                "int": int,
                "float": float,
                "bool": bool,
                "list": list,
            }

            if type_str in type_map:
                t = type_map[type_str]
            else:
                # Not a type prefix, restore original
                o = option

        return t, o

    def _process_option(self, option: str, value: str, section: str = "GENERAL") -> Setting:
        """
        Process a single option.

        Parameters
        ----------
        option : str
            Option name.
        value : str
            Option value string.
        section : str
            Section name.

        Returns
        -------
        Setting
            Processed setting.
        """
        s = Setting(section=section)

        try:
            # Split off comments
            parts = value.split("#", 1)
            s.value = parts[0].strip() or None

            if len(parts) > 1:
                s.comment = parts[1].strip()

            # Extract type and clean name
            vtype, s.name = self._extract_type(option)
            s.type_cast(vtype)

            # Validate paths
            if section == "PATHS" and s.value:
                try:
                    s.value = valid_path(s.value)
                except Exception as e:
                    s.error = str(e)

            # Validate binaries
            if section == "BINARIES" and s.value:
                try:
                    s.value = valid_binary(s.value)
                except Exception as e:
                    s.error = str(e)

        except SettingsWarning as e:
            s.error = str(e)

        return s

    def _process_section(
        self, items: list[tuple[str, str]], section: str = "GENERAL", verbose: bool = False
    ) -> dict[str, Setting]:
        """
        Process all items in a section.

        Parameters
        ----------
        items : list[tuple[str, str]]
            ConfigParser items.
        section : str
            Section name.
        verbose : bool
            Log warnings.

        Returns
        -------
        dict[str, Setting]
            Dictionary of settings.
        """
        result = {}

        for name, value in items:
            s = self._process_option(name, value, section)
            result[s.name] = s

            if verbose and s.error:
                self.log.warning(s.error)

        return result

    def parse(self, keep_sections: bool = False) -> dict[str, Any]:
        """
        Parse the configuration file.

        Parameters
        ----------
        keep_sections : bool
            If True, preserve section structure in result.

        Returns
        -------
        dict
            Dictionary of settings. If keep_sections=True, nested by section.

        Raises
        ------
        IOError
            If file not found.
        InvalidFile
            If parsing fails.
        """
        try:
            c = CaseSensitiveConfigParser()

            if not self.f_ini.exists():
                raise OSError(f"Settings file not found: {self.f_ini}")

            c.read(self.f_ini)

            for section in c.sections():
                section_result = self._process_section(c.items(section), section)

                if keep_sections:
                    self.result[section] = section_result
                else:
                    self.result.update(section_result)

        except configparser.Error as e:
            raise InvalidFile(f"Error parsing settings file {self.f_ini}: {e}") from e

        return self.result

    def __repr__(self) -> str:
        errors = sum(1 for s in self.result.values() if isinstance(s, Setting) and s.error)
        return f"SettingsParser({self.f_ini}, entries={len(self.result)}, errors={errors})"


# =============================================================================
# Settings Manager
# =============================================================================


class SettingsManager:
    """
    Merge default and user configuration files.

    Combines:
    1. Default configuration (package defaults)
    2. User configuration (~/.mdmix/)

    User settings override defaults where valid.

    Parameters
    ----------
    f_default : str | Path
        Default configuration file.
    f_user : str | Path
        User configuration file.
    create_missing : bool
        Create user config if missing.
    verbose : bool
        Print warnings.

    Examples
    --------
    >>> manager = SettingsManager('defaults.cfg', '~/.mdmix/settings.cfg')
    >>> settings = manager.settings_to_dict()
    """

    USER_HEADER = """
##     pyMDMix User Configuration File
##
##     Parameters here override defaults in %(f_default)s.
##
##     Type prefixes: int-, float-, bool-, list-
##     Example: int-nreplicas = 3
##              list-solvents = ETA, WAT, MAM

"""

    def __init__(
        self,
        f_default: str | Path,
        f_user: str | Path,
        create_missing: bool = False,
        verbose: bool = True,
    ):
        self.log = logging.getLogger("SettingsManager")
        self.verbose = verbose
        self.f_default = Path(f_default)
        self.f_user = absfile(f_user)
        self.create_missing = create_missing
        self.user_missing = not self.f_user.exists()

        self.settings: dict[str, Setting] = {}

    def _update(
        self, cfg_default: dict[str, Setting], cfg_user: dict[str, Setting]
    ) -> dict[str, Setting]:
        """
        Override defaults with valid user settings.

        Parameters
        ----------
        cfg_default : dict
            Default settings.
        cfg_user : dict
            User settings.

        Returns
        -------
        dict
            Merged settings.
        """
        result = {}

        for name, default in cfg_default.items():
            user_setting = cfg_user.get(name, default)

            # Only use user setting if not more broken than default
            if user_setting.error and not default.error:
                if self.verbose:
                    self.log.warning(
                        f"User setting '{name}' reset to default ({default.value!r}), "
                        f"reason: {user_setting.error}"
                    )
                user_setting = default

            result[name] = user_setting

        return result

    def collect_settings(self) -> None:
        """Parse and combine default and user config files."""
        try:
            # Parse defaults
            p_default = SettingsParser(self.f_default)
            c_default = p_default.parse()

            # Parse user config (if exists)
            c_user: dict[str, Setting] = {}
            if self.f_user.exists():
                try:
                    p_user = SettingsParser(self.f_user)
                    c_user = p_user.parse()
                except OSError:
                    if self.verbose:
                        self.log.warning(f"User settings not found: {self.f_user}")

            # Merge
            self.settings = self._update(c_default, c_user)

        except SettingsError as e:
            self.log.error(str(e))
            raise

    def write_user_settings(self, errors_only: bool = False) -> None:
        """
        Write user configuration file.

        Parameters
        ----------
        errors_only : bool
            Only write settings with errors.
        """
        try:
            backup(self.f_user)

            # Create parent directory
            self.f_user.parent.mkdir(parents=True, exist_ok=True)

            # Group by section
            sections = ["GENERAL", "PATHS", "BINARIES"]
            by_section: dict[str, list[Setting]] = {s: [] for s in sections}

            for s in self.settings.values():
                section = s.section if s.section in sections else "GENERAL"
                by_section[section].append(s)

            for section in by_section:
                by_section[section].sort()

            # Write file
            with open(self.f_user, "w") as f:
                f.write(self.USER_HEADER % {"f_default": self.f_default})

                for section in sections:
                    settings = by_section[section]
                    if not settings:
                        continue

                    f.write(f"[{section}]\n\n")

                    for param in settings:
                        if errors_only and not param.error:
                            f.write(f"## {param.formatted()}\n")
                        else:
                            f.write(f"{param.formatted()}\n")

                    f.write("\n")

        except OSError as e:
            raise WriteCfgError(str(e)) from e

    def settings_to_dict(self) -> dict[str, Any]:
        """
        Convert settings to plain dictionary.

        Returns
        -------
        dict
            Dictionary of {name: value}.
        """
        return {s.name: s.value for s in self.settings.values()}

    def update_namespace(self, ns: dict[str, Any], keep_defined: bool = True) -> None:
        """
        Insert settings into a namespace (like module locals()).

        Parameters
        ----------
        ns : dict
            Namespace dictionary (e.g., from locals()).
        keep_defined : bool
            Don't overwrite existing variables.
        """
        self.collect_settings()

        if self.user_missing and self.create_missing:
            if self.verbose:
                self.log.warning(f"Creating user config: {self.f_user}")
            self.write_user_settings(errors_only=True)

        d = self.settings_to_dict()

        if keep_defined:
            for k, v in d.items():
                if k not in ns:
                    ns[k] = v
        else:
            ns.update(d)


if __name__ == "__main__":
    # Simple test
    import tempfile

    # Create test config
    cfg_content = """
[GENERAL]
int-count = 42
float-rate = 3.14
bool-enabled = true
list-items = one, two, three
name = test_name  # A comment

[PATHS]
data_dir = /tmp
"""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".cfg", delete=False) as f:
        f.write(cfg_content)
        cfg_path = f.name

    try:
        parser = SettingsParser(cfg_path)
        result = parser.parse()

        print("Parsed settings:")
        for name, setting in result.items():
            print(f"  {name}: {setting.value!r} ({setting.vtype.__name__})")

        print(f"\nParser: {parser}")
    finally:
        os.unlink(cfg_path)
