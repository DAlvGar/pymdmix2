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
from typing import Any

log = logging.getLogger(__name__)


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

    Supports per-replica customization using the slash notation:
        RESTR = HA, WAT//FREE, ETA/1/BB

    Examples
    --------
    >>> parser = MDSettingsConfigFileParser()
    >>> settings_list = parser.parse('md_settings.cfg')
    >>> len(settings_list)
    3
    """

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
        default_nreplicas: int = 3,
    ) -> list[dict[str, Any]]:
        """
        Parse MD settings configuration file.

        Parameters
        ----------
        config_file : str | Path
            Path to configuration file.
        md_settings_keys : list[str] | None
            Valid MD settings parameter names (for fuzzy matching).
        default_nreplicas : int
            Default number of replicas.

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
            if md_settings_keys:
                for k, v in section_dict.items():
                    if v is None or k in main_opts:
                        continue

                    # Fuzzy match to known settings
                    matches = get_close_matches(k, md_settings_keys, n=1, cutoff=0.8)
                    if matches:
                        extra_cfg[matches[0]] = v
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
        if md_settings_keys:
            for k, v in section_dict.items():
                if v is None or k in main_opts:
                    continue
                matches = get_close_matches(k, md_settings_keys, n=1, cutoff=0.8)
                if matches:
                    extra_cfg[matches[0]] = v
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
        unit_name=params.get("uname"),
        extra_residues=params.get("extrares"),
        extra_forcefields=params.get("extraff"),
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
