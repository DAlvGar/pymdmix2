"""
pyMDMix Utility Functions
=========================

Provides path utilities, mask parsing, logging, and other helper functions.

Migrated from the original pyMDMix/tools.py and Biskit.tools.
"""

from __future__ import annotations

import logging
import os
import shutil
import sys
import traceback
from itertools import groupby
from operator import itemgetter
from pathlib import Path
from types import FrameType
from typing import Any

# =============================================================================
# Exceptions
# =============================================================================


class ToolsError(Exception):
    """Base exception for tools module."""

    pass


class InvalidPath(ToolsError):
    """Raised when a path is invalid or doesn't exist."""

    pass


class InvalidBinary(ToolsError):
    """Raised when a binary path is invalid."""

    pass


# =============================================================================
# Path Utilities
# =============================================================================


def project_root(file: str | None = None) -> Path:
    """
    Get the root directory of the pymdmix package.

    Parameters
    ----------
    file : str | None
        Optional file to join with the root path.

    Returns
    -------
    Path
        Absolute path to package root (or file within it).

    Examples
    --------
    >>> project_root()
    PosixPath('/path/to/pymdmix')
    >>> project_root('data/solvents')
    PosixPath('/path/to/pymdmix/data/solvents')
    """
    import pymdmix

    root = Path(pymdmix.__file__).parent.resolve()
    if file:
        return root / file
    return root


def data_root(subfolder: str = "") -> Path:
    """
    Get the data directory of the pymdmix package.

    Parameters
    ----------
    subfolder : str
        Optional subfolder within data directory.

    Returns
    -------
    Path
        Absolute path to data directory.

    Examples
    --------
    >>> data_root()
    PosixPath('/path/to/pymdmix/data')
    >>> data_root('solvents')
    PosixPath('/path/to/pymdmix/data/solvents')
    """
    return project_root("data") / subfolder if subfolder else project_root("data")


def test_root(*subpath: str) -> Path:
    """
    Get the test data directory.

    Parameters
    ----------
    *subpath : str
        Subpath components to join.

    Returns
    -------
    Path
        Absolute path to test data directory.

    Examples
    --------
    >>> test_root()
    PosixPath('/path/to/pymdmix/data/test')
    >>> test_root('pep', 'system.pdb')
    PosixPath('/path/to/pymdmix/data/test/pep/system.pdb')
    """
    base = data_root("test")
    if subpath:
        return base.joinpath(*subpath)
    return base


def templates_root(file: str | None = None) -> Path:
    """
    Get the templates directory.

    Parameters
    ----------
    file : str | None
        Optional filename within templates directory.

    Returns
    -------
    Path
        Absolute path to templates directory.
    """
    base = data_root("templates")
    if file:
        return base / file
    return base


def solvents_root(file: str | None = None) -> Path:
    """
    Get the solvents library directory.

    Parameters
    ----------
    file : str | None
        Optional filename within solvents directory.

    Returns
    -------
    Path
        Absolute path to solvents directory.
    """
    base = data_root("solvents")
    if file:
        return base / file
    return base


def absfile(filename: str | Path) -> Path:
    """
    Get absolute path with ~ expansion.

    Parameters
    ----------
    filename : str | Path
        Filename or path.

    Returns
    -------
    Path
        Absolute path.
    """
    return Path(filename).expanduser().resolve()


def valid_path(path: str | Path) -> Path:
    """
    Validate that a path exists.

    Parameters
    ----------
    path : str | Path
        Path to validate.

    Returns
    -------
    Path
        Validated absolute path.

    Raises
    ------
    InvalidPath
        If path doesn't exist.

    Examples
    --------
    >>> valid_path('/tmp')
    PosixPath('/tmp')
    >>> valid_path('/nonexistent')  # doctest: +IGNORE_EXCEPTION_DETAIL
    Traceback (most recent call last):
        ...
    InvalidPath: Invalid path: '/nonexistent'
    """
    try:
        path = absfile(path)
        if not path.exists():
            raise InvalidPath(f"Invalid path: '{path}'")
        return path
    except InvalidPath:
        raise
    except Exception as e:
        raise InvalidPath(f"Error validating path: {e}") from e


def valid_binary(name: str) -> Path:
    """
    Validate that a binary exists and is executable.

    Parameters
    ----------
    name : str
        Binary name or path.

    Returns
    -------
    Path
        Absolute path to binary.

    Raises
    ------
    InvalidBinary
        If binary not found.

    Examples
    --------
    >>> valid_binary('python')  # doctest: +SKIP
    PosixPath('/usr/bin/python')
    """
    if not name:
        raise InvalidBinary("Empty binary name")

    # Try as path first
    path = Path(name).expanduser()
    if path.exists() and path.is_file():
        if os.access(path, os.X_OK):
            return path.resolve()

    # Try to find in PATH
    binary = shutil.which(name)
    if binary:
        return Path(binary).resolve()

    raise InvalidBinary(f"Binary not found: '{name}'")


def file_permissions(path: str | Path) -> dict[str, bool]:
    """
    Check file permissions.

    Parameters
    ----------
    path : str | Path
        Path to check.

    Returns
    -------
    dict[str, bool]
        Dictionary with keys: EXISTS, READ, WRITE, EXECUTE.

    Examples
    --------
    >>> perms = file_permissions('/tmp')
    >>> perms['EXISTS']
    True
    """
    path = Path(path)
    return {
        "EXISTS": path.exists(),
        "READ": os.access(path, os.R_OK),
        "WRITE": os.access(path, os.W_OK),
        "EXECUTE": os.access(path, os.X_OK),
    }


def backup(path: str | Path, suffix: str = ".bak") -> Path | None:
    """
    Create backup of a file if it exists.

    Parameters
    ----------
    path : str | Path
        File to backup.
    suffix : str
        Suffix for backup file.

    Returns
    -------
    Path | None
        Backup path if created, None if original didn't exist.
    """
    path = Path(path)
    if path.exists():
        backup_path = path.with_suffix(path.suffix + suffix)
        shutil.copy2(path, backup_path)
        return backup_path
    return None


def try_remove(path: str | Path, tree: bool = False) -> bool:
    """
    Try to remove a file or directory.

    Parameters
    ----------
    path : str | Path
        Path to remove.
    tree : bool
        If True, remove entire directory tree.

    Returns
    -------
    bool
        True if removed, False otherwise.
    """
    path = Path(path)
    try:
        if tree and path.is_dir():
            shutil.rmtree(path)
        elif path.exists():
            path.unlink()
        return True
    except Exception:
        return False


def temp_dir() -> Path:
    """
    Get temporary directory path.

    Returns
    -------
    Path
        Path to temp directory.
    """
    import tempfile

    return Path(tempfile.gettempdir())


# =============================================================================
# List Utilities
# =============================================================================


def simplify_nested_list(nested_list: list) -> list:
    """
    Flatten a nested list into a single list.

    Parameters
    ----------
    nested_list : list
        Nested list of any depth.

    Returns
    -------
    list
        Flattened list.

    Examples
    --------
    >>> simplify_nested_list([[1, 2], [3, [4, 5]]])
    [1, 2, 3, 4, 5]
    >>> simplify_nested_list([1, [2, [3, [4]]]])
    [1, 2, 3, 4]
    """
    result: list[Any] = []

    def _flatten(items):
        for item in items:
            if isinstance(item, list):
                _flatten(item)
            else:
                result.append(item)

    _flatten(nested_list)
    return result


# =============================================================================
# Mask Utilities
# =============================================================================


def amber_mask_to_dict(mask_string: str) -> dict[str, list[str]]:
    """
    Parse AMBER-style residue@atom mask to dictionary.

    Parameters
    ----------
    mask_string : str
        Mask string like "RES@AT1,AT2;RES2@AT3".

    Returns
    -------
    dict[str, list[str]]
        Dictionary mapping residue names to atom name lists.

    Examples
    --------
    >>> amber_mask_to_dict(':ETA@O1,C1;WAT@O')
    {'ETA': ['O1', 'C1'], 'WAT': ['O']}
    >>> amber_mask_to_dict('ETA')
    {'ETA': ['all']}
    """
    # Remove leading colon if present
    if mask_string.startswith(":"):
        mask_string = mask_string[1:]

    result: dict[str, list[str]] = {}
    parts = mask_string.split(";")

    for part in parts:
        if "@" in part:
            residue, atoms = part.split("@", 1)
            atom_list = [a.strip() for a in atoms.split(",")]
            if residue in result:
                result[residue].extend(atom_list)
            else:
                result[residue] = atom_list
        else:
            # Just residue name, select all atoms
            result[part.strip()] = ["all"]

    return result


def parse_num_mask(mask: str) -> list[int]:
    """
    Parse numeric mask to list of integers.

    Parameters
    ----------
    mask : str
        Mask string like "1,2,5:10" or "1-5,7,10-12".

    Returns
    -------
    list[int]
        Sorted list of integers.

    Examples
    --------
    >>> parse_num_mask('1,2,5:10')
    [1, 2, 5, 6, 7, 8, 9, 10]
    >>> parse_num_mask('1-3,7,10-12')
    [1, 2, 3, 7, 10, 11, 12]
    """
    result: list[int] = []
    parts = mask.split(",")

    for part in parts:
        part = part.strip()
        # Support both ':' and '-' as range separators
        if ":" in part:
            a, b = map(int, part.split(":"))
            result.extend(range(a, b + 1))
        elif "-" in part:
            a, b = map(int, part.split("-"))
            result.extend(range(a, b + 1))
        else:
            result.append(int(part))

    result.sort()
    return result


def num_list_to_mask(num_list: list[int]) -> str:
    """
    Convert list of integers to mask string.

    Parameters
    ----------
    num_list : list[int]
        List of integers.

    Returns
    -------
    str
        Mask string with ranges like "1-5,7,8-20".

    Examples
    --------
    >>> num_list_to_mask([1, 2, 3, 7, 10, 11, 12])
    '1-3,7,10-12'
    >>> num_list_to_mask([5])
    '5'
    >>> num_list_to_mask([])
    ''
    """
    if not num_list:
        return ""

    if len(num_list) == 1:
        return str(num_list[0])

    num_list = sorted(set(num_list))
    ranges = []

    for k, g in groupby(enumerate(num_list), lambda x: x[0] - x[1]):
        group = list(map(itemgetter(1), g))
        if len(group) > 1:
            ranges.append(f"{group[0]}-{group[-1]}")
        else:
            ranges.append(str(group[0]))

    return ",".join(ranges)


# =============================================================================
# Logging Utilities
# =============================================================================


class LogFormatter(logging.Formatter):
    """
    Custom log formatter that formats INFO differently from warnings/errors.

    Examples
    --------
    >>> handler = logging.StreamHandler()
    >>> handler.setFormatter(LogFormatter())
    """

    def format(self, record: logging.LogRecord) -> str:
        if record.levelno == logging.INFO:
            self._style._fmt = "%(levelname)s\t%(message)s"
        else:
            self._style._fmt = "###\t%(levelname)s (%(name)s)\t%(message)s"
        return super().format(record)


def traceback_plus() -> str:
    """
    Generate detailed traceback with local variables.

    Returns
    -------
    str
        Formatted traceback string with local variables.

    Examples
    --------
    >>> try:
    ...     x = 1 / 0
    ... except:
    ...     tb = traceback_plus()
    """
    output: list[str] = []

    # Get the exception info
    tb = sys.exc_info()[2]
    if tb is None:
        return "No exception traceback available"

    # Navigate to the last frame
    while tb.tb_next:
        tb = tb.tb_next

    # Build stack from tb_frame upward
    stack: list[FrameType] = []
    f: FrameType | None = tb.tb_frame
    while f:
        stack.append(f)
        f = f.f_back
    stack.reverse()

    # Format standard traceback
    output.append("".join(traceback.format_exc()))
    output.append("\nLocals by frame, innermost last:\n")

    for frame in stack:
        output.append(
            f"\nFrame {frame.f_code.co_name} in {frame.f_code.co_filename} "
            f"at line {frame.f_lineno}\n"
        )
        for key, value in frame.f_locals.items():
            try:
                output.append(f"\t{key:>20} = {value}\n")
            except Exception:
                output.append(f"\t{key:>20} = <ERROR WHILE PRINTING VALUE>\n")

    return "".join(output)


def clip_str(s: str, length: int) -> str:
    """
    Clip a string to maximum length.

    Parameters
    ----------
    s : str
        Input string.
    length : int
        Maximum length.

    Returns
    -------
    str
        Clipped string with '...' if truncated.
    """
    if len(s) <= length:
        return s
    return s[: length - 3] + "..."


# =============================================================================
# Convenience aliases (compatible with old API)
# =============================================================================

projectRoot = project_root
dataRoot = data_root
testRoot = test_root
templatesRoot = templates_root
solventsRoot = solvents_root
validPath = valid_path
validBinary = valid_binary
filePermission = file_permissions
simplifyNestedList = simplify_nested_list
amberMaskToDict = amber_mask_to_dict
parseNumMask = parse_num_mask
numListToMask = num_list_to_mask
tracebackPlus = traceback_plus


if __name__ == "__main__":
    print(f"Project root: {project_root()}")
    print(f"Data root: {data_root()}")
    print(f"Solvents root: {solvents_root()}")

    # Test nested list flattening
    nested = [[1, 2, 3], 4, [[5, 6], [[7], [8], [9, 10, 11]], [[12, 13], 14]]]
    print(f"Flattened: {simplify_nested_list(nested)}")

    # Test mask parsing
    print(f"Parse '1,2,5:10': {parse_num_mask('1,2,5:10')}")
    print(f"Mask from [1,2,3,7,10,11,12]: {num_list_to_mask([1, 2, 3, 7, 10, 11, 12])}")
