"""
pyMDMix Utilities
=================

Utility functions for path handling, mask parsing, logging, and configuration.
"""
from __future__ import annotations

from pymdmix.utils.settings_parser import (
    CaseSensitiveConfigParser,
    InvalidFile,
    InvalidType,
    InvalidValue,
    Setting,
    SettingsError,
    SettingsManager,
    SettingsParser,
    SettingsWarning,
    WriteCfgError,
)
from pymdmix.utils.tools import (
    # Exceptions
    InvalidBinary,
    InvalidPath,
    LogFormatter,
    ToolsError,
    absfile,
    # Mask utilities
    amber_mask_to_dict,
    amberMaskToDict,
    backup,
    clip_str,
    data_root,
    dataRoot,
    file_permissions,
    filePermission,
    num_list_to_mask,
    numListToMask,
    parse_num_mask,
    parseNumMask,
    # Path utilities
    project_root,
    # Legacy aliases
    projectRoot,
    # List utilities
    simplify_nested_list,
    simplifyNestedList,
    solvents_root,
    solventsRoot,
    temp_dir,
    templates_root,
    templatesRoot,
    test_root,
    testRoot,
    traceback_plus,
    tracebackPlus,
    try_remove,
    valid_binary,
    valid_path,
    validBinary,
    validPath,
)

__all__ = [
    # Path utilities
    "project_root",
    "data_root",
    "test_root",
    "templates_root",
    "solvents_root",
    "absfile",
    "valid_path",
    "valid_binary",
    "file_permissions",
    "backup",
    "try_remove",
    "temp_dir",
    # List utilities
    "simplify_nested_list",
    # Mask utilities
    "amber_mask_to_dict",
    "parse_num_mask",
    "num_list_to_mask",
    # Logging utilities
    "LogFormatter",
    "traceback_plus",
    "clip_str",
    # Settings
    "Setting",
    "SettingsParser",
    "SettingsManager",
    "CaseSensitiveConfigParser",
    # Exceptions
    "ToolsError",
    "InvalidPath",
    "InvalidBinary",
    "SettingsError",
    "InvalidType",
    "InvalidValue",
    "InvalidFile",
    "SettingsWarning",
    "WriteCfgError",
    # Legacy aliases (camelCase)
    "projectRoot",
    "dataRoot",
    "testRoot",
    "templatesRoot",
    "solventsRoot",
    "validPath",
    "validBinary",
    "filePermission",
    "simplifyNestedList",
    "amberMaskToDict",
    "parseNumMask",
    "numListToMask",
    "tracebackPlus",
]
