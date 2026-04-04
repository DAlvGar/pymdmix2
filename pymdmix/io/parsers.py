"""
Configuration File Parsers
==========================

Parse system and MD settings configuration files.

Migrated from the original pyMDMix/Parsers.py.
"""

from __future__ import annotations

import configparser
import logging
from dataclasses import dataclass
from difflib import get_close_matches
from pathlib import Path
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from pymdmix.project.settings import MDSettings

log = logging.getLogger(__name__)


# =============================================================================
# Legacy field-name mapping  (old cfg key → MDSettings dataclass field name)
# =============================================================================

#: Maps legacy/old-style cfg key names to the corresponding ``MDSettings``
#: dataclass field names used in pymdmix2.
_LEGACY_TO_MDSETTINGS: dict[str, str] = {
    # Primary restraint keys (as produced by MDSettingsConfigFileParser.parse)
    "restr_mode": "restraint_mode",
    "restr_force": "restraint_force",
    "restrain_mask": "restraint_mask",
    # Old-style config file names
    "restrmode": "restraint_mode",
    "restrforce": "restraint_force",
    "restrmask": "restraint_mask",
    "alignmask": "align_mask",
    # Temperature aliases
    "temp": "temperature",
    "md_timestep": "timestep",
    # Trajectory / production settings
    "trajfrequency": "trajectory_frequency",
    "prod_steps": "production_steps",
    "npt_eq_steps": "npt_equilibration_steps",
    "minsteps": "minimization_steps",
    "heating_steps": "heating_steps",
    "parm_heating_tempi": "heating_initial_temp",
    "npt_pressure": "npt_pressure",
    # Format flags
    "mdnetcdf": "use_netcdf",
    "iwrap": "wrap_coordinates",
    # Program / ensemble
    "mdprogram": "md_program",
    "production_ensemble": "production_ensemble",
    # Folder names
    "mdfolder": "md_folder",
    "eqfolder": "eq_folder",
    "minfolder": "min_folder",
    "alignfolder": "align_folder",
    "densityfolder": "density_folder",
    "energyfolder": "energy_folder",
    # File templates
    "mdoutfiletemplate": "md_file_template",
    "eqoutfiletemplate": "eq_file_template",
    # Force fields
    "ff": "force_fields",
}

# =============================================================================
# Exceptions
# =============================================================================


class ParserError(Exception):
    """Base exception for parsers."""

    pass


class SystemParserError(ParserError):
    """Error parsing system configuration."""

    pass


class BadFile(ParserError):
    """File not found or invalid."""

    pass


class MDSettingsParserError(ParserError):
    """Error parsing MD settings configuration."""

    pass


class BadSolvent(MDSettingsParserError):
    """Invalid solvent specified."""

    pass


# =============================================================================
# System Configuration Parser
# =============================================================================


class SystemConfigFileParser:
    """
    Parse system configuration files.

    Reads a config file with [SYSTEM] section and returns a dictionary
    of parameters suitable for creating a System object.

    Examples
    --------
    >>> parser = SystemConfigFileParser()
    >>> params = parser.parse('system.cfg')
    >>> params['name']
    'my_system'
    """

    def __init__(self):
        self.log = logging.getLogger("SystemConfigFileParser")

    def parse(self, config_file: str | Path) -> dict[str, Any]:
        """
        Parse system configuration file.

        Parameters
        ----------
        config_file : str | Path
            Path to configuration file.

        Returns
        -------
        dict
            Dictionary with system parameters:
            - name: System name
            - pdb: PDB file path (if provided)
            - off: OFF file path (if provided)
            - top: PRMTOP file path (if provided)
            - crd: PRMCRD file path (if provided)
            - unit_name: Unit name in OFF file
            - extra_residues: List of extra residue names
            - restrain_mask: Restraint mask
            - align_mask: Alignment mask
            - extra_ff: List of extra force field files

        Raises
        ------
        BadFile
            If file not found.
        SystemParserError
            If [SYSTEM] section missing or invalid.
        """
        config_file = Path(config_file)

        if not config_file.exists():
            raise BadFile(f"File not found: {config_file}")

        config = configparser.ConfigParser()
        config.read(config_file)

        # Check for SYSTEM section
        if "SYSTEM" not in config.sections():
            raise SystemParserError("Section [SYSTEM] missing.")

        section = dict(config.items("SYSTEM"))
        params: dict[str, Any] = {}
        config_dir = config_file.parent

        # System name
        params["name"] = section.get("name", "system")

        # Input files (priority: off > pdb > top+crd)
        off = section.get("off")
        pdb = section.get("pdb")
        top = section.get("top")
        crd = section.get("crd")

        if off:
            # OFF file has priority
            if pdb:
                self.log.warning("Ignoring PDB entry, OFF has priority.")

            if off != "test":  # 'test' is special case
                off_path = Path(off)
                if not off_path.is_absolute():
                    off_path = config_dir / off
                if not off_path.exists():
                    raise BadFile(f"OFF file not found: {off}")
                params["off"] = str(off_path.resolve())
            else:
                params["off"] = "test"

            params["unit_name"] = section.get("uname")

        elif pdb:
            pdb_path = Path(pdb)
            if not pdb_path.is_absolute():
                pdb_path = config_dir / pdb
            if not pdb_path.exists():
                raise BadFile(f"PDB file not found: {pdb}")
            params["pdb"] = str(pdb_path.resolve())

        elif top and crd:
            top_path = Path(top)
            crd_path = Path(crd)
            if not top_path.is_absolute():
                top_path = config_dir / top
            if not crd_path.is_absolute():
                crd_path = config_dir / crd

            if not top_path.exists():
                raise BadFile(f"PRMTOP file not found: {top}")
            if not crd_path.exists():
                raise BadFile(f"PRMCRD file not found: {crd}")

            params["top"] = str(top_path.resolve())
            params["crd"] = str(crd_path.resolve())

        # Extra residues
        extra_res = section.get("extrares", "")
        if extra_res:
            params["extra_residues"] = [r.strip() for r in extra_res.split(",")]
            self.log.info(f"Extra residue list: {params['extra_residues']}")
        else:
            params["extra_residues"] = []

        # Restraint and alignment masks
        params["restrain_mask"] = section.get("restrmask", "auto")
        params["align_mask"] = section.get("alignmask", "auto")

        # Extra force field files
        extra_ff = section.get("extraff", "")
        if extra_ff:
            params["extra_ff"] = [f.strip() for f in extra_ff.split(",")]
        else:
            params["extra_ff"] = []

        return params


# =============================================================================
# MD Settings Configuration Parser
# =============================================================================


class MDSettingsConfigFileParser:
    """
    Parse MD settings configuration files.

    Reads config files with [MDSETTINGS] sections and returns a list of
    parameter dictionaries for creating MDSettings objects.

    Default values for extra configuration options are loaded automatically
    from the package defaults file (``data/defaults/md-settings.cfg``), and
    user-supplied overrides in ``~/.mdmix/md-settings.cfg`` are applied on top
    — mirroring the behaviour of the original pyMDMix ``SettingsManager``.

    Supports per-replica customisation using the slash notation::

        RESTR = HA, WAT//FREE, ETA/1/BB

    Examples
    --------
    >>> parser = MDSettingsConfigFileParser()
    >>> settings_list = parser.parse('md_settings.cfg')
    >>> len(settings_list)
    3
    """

    #: Path to the user-level MD settings override file.
    USER_MD_CFG: Path = Path("~/.mdmix/md-settings.cfg").expanduser()

    def __init__(self, available_solvents: list[str] | None = None):
        """
        Initialize parser.

        Parameters
        ----------
        available_solvents : list[str] | None
            List of valid solvent names for validation.
        """
        self.log = logging.getLogger("MDSettingsConfigFileParser")
        self.available_solvents = available_solvents or []
        self.solvents: list[str] = []
        self.solv_nrepl: dict[str, int] = {}
        self._defaults_cache: dict[str, Any] | None = None

    def _make_settings_manager(self, default_cfg: Path, user_cfg: Path) -> Any:
        """Create a :class:`SettingsManager` instance for the given files."""
        from pymdmix.utils.settings_parser import SettingsManager

        return SettingsManager(
            f_default=default_cfg,
            f_user=user_cfg,
            create_missing=False,
            verbose=False,
        )

    def _load_defaults(self) -> dict[str, Any]:
        """
        Load and cache defaults from ``data/defaults/md-settings.cfg``.

        Returns
        -------
        dict
            Mapping of setting name → value, with types applied.
        """
        if self._defaults_cache is not None:
            return self._defaults_cache

        from pymdmix.utils.tools import defaults_root

        default_cfg = defaults_root("md-settings.cfg")
        if not default_cfg.exists():
            self._defaults_cache = {}
            return self._defaults_cache

        manager = self._make_settings_manager(default_cfg, self.USER_MD_CFG)
        manager.collect_settings()
        self._defaults_cache = manager.settings_to_dict()
        return self._defaults_cache

    def _get_default_keys(self) -> list[str]:
        """Return the list of valid key names from the defaults file."""
        return list(self._load_defaults().keys())

    def _get_default_nreplicas(self) -> int:
        """
        Return the default replica count from ``data/defaults/settings.cfg``.

        Reads ``DEF_NREPLICAS`` from the global settings defaults (falling back
        to 1 if the file or key is unavailable), mirroring the old
        ``S.DEF_NREPLICAS`` access pattern.
        """
        from pymdmix.utils.tools import defaults_root

        default_cfg = defaults_root("settings.cfg")
        user_cfg = Path("~/.mdmix/settings.cfg").expanduser()

        if not default_cfg.exists():
            return 1

        try:
            manager = self._make_settings_manager(default_cfg, user_cfg)
            manager.collect_settings()
            return int(manager.settings_to_dict().get("DEF_NREPLICAS", 1))
        except Exception:
            return 1

    @staticmethod
    def _map_legacy_key(key: str) -> str:
        """
        Translate a legacy/old-style cfg key to its ``MDSettings`` field name.

        Parameters
        ----------
        key : str
            Key name from the cfg file (may be camelCase or old-style snake_case).

        Returns
        -------
        str
            Canonical ``MDSettings`` dataclass field name, or the original key
            if no mapping exists.
        """
        return _LEGACY_TO_MDSETTINGS.get(key.lower(), key)

    def _validate_solvents(self, solvents: list[str]) -> None:
        """
        Validate solvent list.

        Raises
        ------
        BadSolvent
            If any solvent is not in available list.
        """
        if not self.available_solvents:
            return  # No validation if no list provided

        missing = set(solvents) - set(self.available_solvents)
        if missing:
            raise BadSolvent(f"Invalid solvents: {missing}")

    def _split_per_replica_slash(
        self, string: str | None, val_control: list[str] | None = None, val_format: type = str
    ) -> dict[str, dict[int, Any]]:
        """
        Parse per-replica values using slash notation.

        Format: "value" or "SOLV//value" or "SOLV/repl/value"

        Parameters
        ----------
        string : str | None
            Input string.
        val_control : list[str] | None
            Allowed values (for validation).
        val_format : type
            Type to convert values to.

        Returns
        -------
        dict
            {solvent: {replica_num: value}}
        """
        # Initialize result for all solvent/replica combos
        result: dict[str, dict[int, Any]] = {}
        for solv, nrepl in self.solv_nrepl.items():
            result[solv] = {i: None for i in range(1, nrepl + 1)}

        if not string:
            return result

        string = string.strip()
        temp: dict[str, Any] = {}

        if "," in string:
            parts = string.split(",")
            for part in parts:
                part = part.strip()
                if "/" in part:
                    pieces = part.split("/")
                    solv = pieces[0].strip()
                    repl_spec = pieces[1].strip() if len(pieces) > 2 else ""
                    value = val_format(pieces[-1].strip())

                    if val_control and value not in val_control:
                        raise MDSettingsParserError(
                            f"Value '{value}' not valid. Options: {val_control}"
                        )

                    if solv not in temp:
                        temp[solv] = {}

                    # Parse replica spec
                    if not repl_spec:
                        replicas = list(range(1, self.solv_nrepl.get(solv, 1) + 1))
                    elif "-" in repl_spec:
                        start, end = map(int, repl_spec.split("-"))
                        replicas = list(range(start, end + 1))
                    else:
                        replicas = [int(repl_spec)]

                    for r in replicas:
                        temp[solv][r] = value
                else:
                    # Common value
                    value = val_format(part)
                    if val_control and value not in val_control:
                        raise MDSettingsParserError(
                            f"Value '{value}' not valid. Options: {val_control}"
                        )
                    temp["COMMON"] = value
        else:
            value = val_format(string)
            if val_control and value not in val_control:
                raise MDSettingsParserError(f"Value '{value}' not valid. Options: {val_control}")
            temp["COMMON"] = value

        # Apply common and specific values
        common = temp.get("COMMON")
        for solv, nrepl in self.solv_nrepl.items():
            for i in range(1, nrepl + 1):
                if solv in temp and i in temp[solv]:
                    result[solv][i] = temp[solv][i]
                elif common is not None:
                    result[solv][i] = common

        return result

    def _split_per_replica_colon(self, string: str, val_format: type = int) -> dict[str, Any]:
        """
        Parse per-solvent values using colon notation.

        Format: "value" or "SOLV:value,SOLV2:value"

        Parameters
        ----------
        string : str
            Input string.
        val_format : type
            Type to convert values to.

        Returns
        -------
        dict
            {solvent: value}
        """
        string = string.strip()
        temp: dict[str, Any] = {}

        if "," in string:
            parts = string.split(",")
            for part in parts:
                part = part.strip()
                if ":" in part:
                    solv, val = part.split(":")
                    temp[solv.strip()] = val_format(val.strip())
                else:
                    temp["COMMON"] = val_format(part)
        else:
            temp["COMMON"] = val_format(string)

        # Assign to all solvents
        result = {}
        common = temp.get("COMMON")
        for solv in self.solvents:
            if solv in temp:
                result[solv] = temp[solv]
            elif common is not None:
                result[solv] = common
            else:
                raise MDSettingsParserError(f"No value for solvent {solv} and no common value")

        return result

    def parse(
        self,
        config_file: str | Path,
        md_settings_keys: list[str] | None = None,
        default_nreplicas: int | None = None,
    ) -> list[dict[str, Any]]:
        """
        Parse MD settings configuration file.

        Parameters
        ----------
        config_file : str | Path
            Path to configuration file.
        md_settings_keys : list[str] | None
            Valid MD settings parameter names (for fuzzy matching).  When
            ``None`` the list is loaded automatically from the package
            defaults file ``data/defaults/md-settings.cfg``.
        default_nreplicas : int | None
            Default number of replicas when ``NREPL`` is absent.  When
            ``None`` the value is read from ``data/defaults/settings.cfg``
            (``DEF_NREPLICAS``), falling back to 1.

        Returns
        -------
        list[dict]
            List of parameter dictionaries, one per solvent/replica combo.

        Raises
        ------
        BadFile
            If file not found.
        MDSettingsParserError
            If parsing fails.
        """
        config_file = Path(config_file)

        if not config_file.exists():
            raise BadFile(f"Config file not found: {config_file}")

        # Resolve default replica count from package settings when not provided
        if default_nreplicas is None:
            default_nreplicas = self._get_default_nreplicas()

        config = configparser.ConfigParser()
        config.read(config_file)

        # Find MDSETTINGS sections
        sections = config.sections()
        md_sections = [s for s in sections if s.startswith("MDSETTINGS")]

        if not md_sections:
            raise MDSettingsParserError("No sections starting with 'MDSETTINGS' found.")

        settings_list: list[dict[str, Any]] = []

        for section in md_sections:
            section_dict = dict(config.items(section))

            # Get solvents
            solvents_str = section_dict.get("solvents") or section_dict.get("solvent")
            if not solvents_str:
                raise MDSettingsParserError("SOLVENT(S) option missing in config file.")

            self.solvents = [s.strip() for s in solvents_str.split(",")]
            self._validate_solvents(self.solvents)

            # Number of replicas (per-solvent supported)
            nrepl_str = section_dict.get("nrepl", str(default_nreplicas))
            self.solv_nrepl = self._split_per_replica_colon(nrepl_str)

            # Per-replica options
            restr = section_dict.get("restr")
            if restr:
                restr = restr.upper().strip()
            restr_split = self._split_per_replica_slash(restr, val_control=["BB", "HA", "FREE"])

            force_str = section_dict.get("force")
            force_split = self._split_per_replica_slash(force_str, val_format=float)

            nanos_str = section_dict.get("nanos")
            nanos_split = self._split_per_replica_slash(nanos_str, val_format=int)

            temp_str = section_dict.get("temp")
            temp_split = self._split_per_replica_slash(temp_str, val_format=float)

            # Common masks
            restrain_mask = section_dict.get("restrmask", "")
            align_mask = section_dict.get("alignmask", "")

            # Parse extra configuration options
            main_opts = {
                "restr",
                "solvents",
                "solvent",
                "nanos",
                "temp",
                "force",
                "nrepl",
                "restrmask",
                "alignmask",
            }

            extra_cfg: dict[str, Any] = {}
            # Use caller-supplied keys, or fall back to defaults from md-settings.cfg
            effective_keys = md_settings_keys if md_settings_keys is not None else self._get_default_keys()
            defaults = self._load_defaults()

            if effective_keys:
                for k, v in section_dict.items():
                    if v is None or k in main_opts:
                        continue

                    # Fuzzy match to known settings (defaults or caller-supplied)
                    matches = get_close_matches(k, effective_keys, n=1, cutoff=0.8)
                    if matches:
                        matched_key = matches[0]
                        # Apply type from defaults if available
                        if matched_key in defaults and isinstance(defaults[matched_key], (int, float, bool)):
                            try:
                                v = type(defaults[matched_key])(v)
                            except (ValueError, TypeError):
                                pass
                        extra_cfg[matched_key] = v
                    else:
                        raise MDSettingsParserError(f"Unknown attribute '{k}' in MD settings")

            # Generate settings for each solvent/replica combo
            for solv, nrepl in self.solv_nrepl.items():
                for i in range(1, nrepl + 1):
                    params = {
                        "solvent": solv,
                        "nanos": nanos_split[solv][i],
                        "restr_mode": restr_split[solv][i],
                        "restr_force": force_split[solv][i],
                        "temp": temp_split[solv][i],
                        "restrain_mask": restrain_mask,
                        "align_mask": align_mask,
                        **extra_cfg,
                    }
                    settings_list.append(params)

        return settings_list

    def parse_to_mdsettings(
        self,
        config_file: str | Path,
        default_nreplicas: int | None = None,
    ) -> list[MDSettings]:
        """
        Parse an MD settings config file and return :class:`MDSettings` objects.

        This is the high-level convenience method that mirrors the old
        ``pyMDMix.MDSettings.parseSettingsConfigFile()`` behaviour.  It:

        1. Calls :meth:`parse` to get raw parameter dicts.
        2. Translates legacy/old-style key names to ``MDSettings`` field names
           using :data:`_LEGACY_TO_MDSETTINGS`.
        3. Instantiates :class:`~pymdmix.project.settings.MDSettings` objects,
           filling in any missing fields from the package defaults.

        Parameters
        ----------
        config_file : str | Path
            Path to configuration file with ``[MDSETTINGS]`` section(s).
        default_nreplicas : int
            Default number of replicas when ``NREPL`` is absent.

        Returns
        -------
        list[MDSettings]
            One ``MDSettings`` object per solvent × replica combination.
        """
        from pymdmix.project.settings import MDSettings

        raw_list = self.parse(config_file, default_nreplicas=default_nreplicas)
        defaults = self._load_defaults()

        result: list[MDSettings] = []
        for raw in raw_list:
            kwargs: dict[str, Any] = {}

            # Merge defaults first (convert legacy names)
            for cfg_key, default_val in defaults.items():
                mds_key = self._map_legacy_key(cfg_key)
                if mds_key in MDSettings.__dataclass_fields__:
                    kwargs[mds_key] = default_val

            # Apply values from the parsed params dict (override defaults)
            for cfg_key, val in raw.items():
                if val is None:
                    continue
                mds_key = self._map_legacy_key(cfg_key)
                if mds_key not in MDSettings.__dataclass_fields__:
                    # Skip keys that don't map to MDSettings fields
                    continue
                # Convert bool-like int fields using the dataclass field annotation
                field_obj = MDSettings.__dataclass_fields__[mds_key]
                field_type = field_obj.type
                try:
                    if field_type is bool or (isinstance(field_type, type) and issubclass(field_type, bool)):
                        kwargs[mds_key] = bool(int(val)) if not isinstance(val, bool) else val
                    else:
                        kwargs[mds_key] = val
                except (TypeError, ValueError):
                    kwargs[mds_key] = val

            # solvent is required
            solvent = raw.get("solvent") or kwargs.get("solvent")
            if not solvent:
                raise MDSettingsParserError("Solvent not found in parsed parameters.")
            kwargs["solvent"] = solvent

            result.append(MDSettings(**{k: v for k, v in kwargs.items() if v is not None}))

        return result

    def parse_no_solvent(
        self, config_file: str | Path, md_settings_keys: list[str] | None = None
    ) -> dict[str, Any]:
        """
        Parse MD settings without solvent information.

        Parameters
        ----------
        config_file : str | Path
            Path to configuration file.
        md_settings_keys : list[str] | None
            Valid MD settings parameter names.

        Returns
        -------
        dict
            Single settings dictionary.
        """
        config_file = Path(config_file)

        if not config_file.exists():
            raise BadFile(f"Config file not found: {config_file}")

        config = configparser.ConfigParser()
        config.read(config_file)

        sections = config.sections()
        md_sections = [s for s in sections if s.startswith("MDSETTINGS")]

        if not md_sections:
            raise MDSettingsParserError("No sections starting with 'MDSETTINGS' found.")

        section_dict = dict(config.items(md_sections[0]))

        restr_raw = section_dict.get("restr")
        restr = restr_raw.upper().strip() if restr_raw else None

        force_raw = section_dict.get("force")
        force = float(force_raw) if force_raw else None

        nanos_raw = section_dict.get("nanos")
        nanos = int(nanos_raw) if nanos_raw else None

        temp_raw = section_dict.get("temp")
        temp = float(temp_raw) if temp_raw else None

        restrain_mask = section_dict.get("restrmask", "")
        align_mask = section_dict.get("alignmask", "")

        main_opts = {"restr", "nanos", "temp", "force", "nrepl", "restrmask", "alignmask"}

        extra_cfg: dict[str, Any] = {}
        effective_keys = md_settings_keys if md_settings_keys is not None else self._get_default_keys()
        defaults = self._load_defaults()

        if effective_keys:
            for k, v in section_dict.items():
                if v is None or k in main_opts:
                    continue
                matches = get_close_matches(k, effective_keys, n=1, cutoff=0.8)
                if matches:
                    matched_key = matches[0]
                    if matched_key in defaults and isinstance(defaults[matched_key], (int, float, bool)):
                        try:
                            v = type(defaults[matched_key])(v)
                        except (ValueError, TypeError):
                            pass
                    extra_cfg[matched_key] = v
                else:
                    raise MDSettingsParserError(f"Unknown attribute '{k}' in MD settings")

        return {
            "nanos": nanos,
            "restr_mode": restr,
            "restr_force": force,
            "temp": temp,
            "restrain_mask": restrain_mask,
            "align_mask": align_mask,
            **extra_cfg,
        }

if __name__ == "__main__":
    import tempfile

    # Test SystemConfigFileParser
    system_cfg = """
[SYSTEM]
name = test_system
pdb = /tmp/test.pdb
restrmask = :1-100@CA
alignmask = auto
"""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".cfg", delete=False) as f:
        f.write(system_cfg)
        system_path = f.name

    try:
        parser = SystemConfigFileParser()
        # Will fail on pdb not existing, but shows structure works
        try:
            params = parser.parse(system_path)
            print(f"System params: {params}")
        except BadFile as e:
            print(f"Expected error (no test.pdb): {e}")
    finally:
        import os

        os.unlink(system_path)

    # Test MDSettingsConfigFileParser
    md_cfg = """
[MDSETTINGS]
solvents = ETA, WAT
nrepl = 2, WAT:1
nanos = 10
restr = HA
temp = 300.0
"""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".cfg", delete=False) as f:
        f.write(md_cfg)
        md_path = f.name

    try:
        md_parser = MDSettingsConfigFileParser()
        settings = md_parser.parse(md_path)
        print(f"\nMD settings ({len(settings)} total):")
        for s in settings:
            print(f"  {s['solvent']}: nanos={s['nanos']}, restr={s['restr_mode']}")
    finally:
        os.unlink(md_path)


# =============================================================================
# Convenience Functions
# =============================================================================


@dataclass
class SystemConfig:
    """Parsed system configuration."""

    name: str
    input_file: Path
    unit_name: str | None = None
    extra_residues: list[str] | None = None
    extra_forcefields: list[str] | None = None


@dataclass
class ReplicaConfig:
    """Parsed replica configuration."""

    system: str
    solvent: str
    nanos: int = 20
    restraint_mode: str = "FREE"
    restraint_mask: str | None = None
    restraint_force: float = 5.0
    align_mask: str | None = None
    md_settings: Path | None = None


def parse_system_config(config_file: str | Path) -> SystemConfig:
    """
    Parse a system configuration file.

    Parameters
    ----------
    config_file : Path
        Path to system.cfg file

    Returns
    -------
    SystemConfig
        Parsed system configuration
    """
    parser = SystemConfigFileParser()
    params = parser.parse(config_file)

    input_path = params.get("off") or params.get("pdb")
    if input_path is None:
        raise SystemParserError("System config requires either 'off' or 'pdb'")

    return SystemConfig(
        name=params["name"],
        input_file=Path(input_path),
        unit_name=params.get("unit_name"),
        extra_residues=params.get("extra_residues"),
        extra_forcefields=params.get("extra_ff"),
    )


def parse_replica_config(config_file: str | Path) -> ReplicaConfig:
    """
    Parse a replica configuration file.

    Parameters
    ----------
    config_file : Path
        Path to replica.cfg file

    Returns
    -------
    ReplicaConfig
        Parsed replica configuration
    """
    config_file = Path(config_file)

    cfg = configparser.ConfigParser()
    cfg.read(config_file)

    if "REPLICA" not in cfg:
        raise ParserError(f"No [REPLICA] section in {config_file}")

    section = cfg["REPLICA"]

    return ReplicaConfig(
        system=section.get("SYSTEM", section.get("system", "")),
        solvent=section.get("SOLVENT", section.get("solvent", "")),
        nanos=int(section.get("NANOS", section.get("nanos", 20))),
        restraint_mode=section.get("RESTRMODE", section.get("restrmode", "FREE")),
        restraint_mask=section.get("RESTRAINMASK", section.get("restrainmask")),
        restraint_force=float(section.get("RESTRFORCE", section.get("restrforce", 5.0))),
        align_mask=section.get("ALIGNMASK", section.get("alignmask")),
        md_settings=Path(section["MDSETTINGS"]) if section.get("MDSETTINGS") else None,
    )


def parse_settings_config_file(
    config_file: str | Path,
    no_solvent: bool = False,
    default_nreplicas: int | None = None,
) -> list[MDSettings] | list[dict[str, Any]]:
    """
    Parse a settings config file and return :class:`~pymdmix.project.settings.MDSettings` objects.

    Convenience wrapper around
    :meth:`MDSettingsConfigFileParser.parse_to_mdsettings` that mirrors
    the old ``pyMDMix.MDSettings.parseSettingsConfigFile()`` API.

    Parameters
    ----------
    config_file : str | Path
        Path to the configuration file containing ``[MDSETTINGS]`` section(s).
    no_solvent : bool
        If ``True``, parse without solvent information (returns a single dict
        instead of a list of ``MDSettings``).  Retained for backward compat.
    default_nreplicas : int | None
        Default replica count when ``NREPL`` is not given.  ``None`` reads the
        value from ``data/defaults/settings.cfg`` (``DEF_NREPLICAS``).

    Returns
    -------
    list[MDSettings]
        One ``MDSettings`` per solvent × replica combination, or a single
        settings dict if ``no_solvent=True``.

    Raises
    ------
    BadFile
        If the config file does not exist.
    MDSettingsParserError
        If parsing fails.
    """
    parser = MDSettingsConfigFileParser()
    if no_solvent:
        return [parser.parse_no_solvent(config_file)]
    return parser.parse_to_mdsettings(config_file, default_nreplicas=default_nreplicas)


@dataclass
class ProjectConfig:
    """
    Parsed full project configuration (combines system + MD settings).

    This is the result of reading a combined project ``.cfg`` file that contains
    both a ``[SYSTEM]`` section and one or more ``[MDSETTINGS]`` sections.
    It represents the complete setup needed to bootstrap a pyMDMix project
    from a single configuration file.

    Attributes
    ----------
    system : SystemConfig
        Parsed system definition.
    settings : list[Any]
        List of :class:`~pymdmix.project.settings.MDSettings` objects, one per
        solvent × replica combination.

    Examples
    --------
    >>> cfg = parse_project_config("project.cfg")
    >>> print(cfg.system.name)
    'MyProtein'
    >>> len(cfg.settings)
    3
    """

    system: SystemConfig
    settings: list[MDSettings]


def parse_project_config(
    config_file: str | Path,
    default_nreplicas: int | None = None,
) -> ProjectConfig:
    """
    Parse a *full* project configuration file.

    A full project config contains both a ``[SYSTEM]`` section and one or more
    ``[MDSETTINGS]`` sections — the canonical format supported by the original
    pyMDMix.  This function is the primary entry-point for project
    initialisation:

    .. code-block:: bash

        pymdmix create project -n myproject -f project.cfg

    The companion CLI command calls this function, then registers the system and
    creates all replicas in a single step.

    Parameters
    ----------
    config_file : str | Path
        Path to the combined ``.cfg`` file with ``[SYSTEM]`` and
        ``[MDSETTINGS]`` sections.
    default_nreplicas : int | None
        Default replica count used when ``NREPL`` is absent.  ``None`` reads
        the value from ``data/defaults/settings.cfg`` (``DEF_NREPLICAS``).

    Returns
    -------
    ProjectConfig
        Parsed project configuration with ``system`` and ``settings`` fields.

    Raises
    ------
    BadFile
        If the config file does not exist.
    SystemParserError
        If the ``[SYSTEM]`` section is missing or invalid.
    MDSettingsParserError
        If no ``[MDSETTINGS]`` section is found or parsing fails.

    Examples
    --------
    >>> cfg = parse_project_config("project.cfg")
    >>> print(cfg.system.name, len(cfg.settings))
    MyProtein 6
    """
    config_file = Path(config_file)
    if not config_file.exists():
        raise BadFile(f"Config file not found: {config_file}")

    system_cfg = parse_system_config(config_file)
    settings_list = parse_settings_config_file(
        config_file, default_nreplicas=default_nreplicas
    )
    return ProjectConfig(system=system_cfg, settings=settings_list)

