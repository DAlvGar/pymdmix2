"""
Tests for pymdmix.utils module.
"""

import os
import tempfile
from pathlib import Path

import pytest

from pymdmix.utils import (
    InvalidBinary,
    # Exceptions
    InvalidPath,
    # Logging
    LogFormatter,
    absfile,
    # Mask utilities
    amber_mask_to_dict,
    backup,
    clip_str,
    data_root,
    file_permissions,
    num_list_to_mask,
    parse_num_mask,
    # Path utilities
    project_root,
    # List utilities
    simplify_nested_list,
    solvents_root,
    temp_dir,
    templates_root,
    test_root,
    traceback_plus,
    try_remove,
    valid_binary,
    valid_path,
)


class TestPathUtilities:
    """Tests for path utility functions."""

    def test_project_root(self):
        """Test project_root returns package directory."""
        root = project_root()
        assert root.exists()
        assert root.name == "pymdmix"
        assert (root / "__init__.py").exists()

    def test_project_root_with_file(self):
        """Test project_root with file argument."""
        path = project_root("data")
        assert path.name == "data"
        assert "pymdmix" in str(path.parent)

    def test_data_root(self):
        """Test data_root returns data directory."""
        root = data_root()
        assert root.name == "data"

    def test_data_root_with_subfolder(self):
        """Test data_root with subfolder."""
        path = data_root("solvents")
        assert path.name == "solvents"

    def test_test_root(self):
        """Test test_root returns test data directory."""
        root = test_root()
        assert root.name == "test"

    def test_test_root_with_subpath(self):
        """Test test_root with subpath components."""
        path = test_root("pep", "system.pdb")
        assert str(path).endswith("pep/system.pdb")

    def test_templates_root(self):
        """Test templates_root returns templates directory."""
        root = templates_root()
        assert root.name == "templates"

    def test_solvents_root(self):
        """Test solvents_root returns solvents directory."""
        root = solvents_root()
        assert root.name == "solvents"

    def test_absfile(self):
        """Test absfile expands and resolves paths."""
        path = absfile("~")
        assert path.is_absolute()
        assert "~" not in str(path)

    def test_valid_path_exists(self):
        """Test valid_path with existing path."""
        result = valid_path("/tmp")
        assert result.exists()
        assert result.is_absolute()

    def test_valid_path_not_exists(self):
        """Test valid_path raises for non-existent path."""
        with pytest.raises(InvalidPath):
            valid_path("/nonexistent/path/that/does/not/exist")

    def test_valid_binary_python(self):
        """Test valid_binary finds python."""
        # Use python3 as it's more universally available
        result = valid_binary("python3")
        assert result.exists()
        assert os.access(result, os.X_OK)

    def test_valid_binary_not_found(self):
        """Test valid_binary raises for non-existent binary."""
        with pytest.raises(InvalidBinary):
            valid_binary("nonexistent_binary_12345")

    def test_file_permissions(self):
        """Test file_permissions returns dict."""
        perms = file_permissions("/tmp")
        assert perms["EXISTS"] is True
        assert perms["READ"] is True
        assert "WRITE" in perms
        assert "EXECUTE" in perms

    def test_file_permissions_nonexistent(self):
        """Test file_permissions for non-existent file."""
        perms = file_permissions("/nonexistent/path")
        assert perms["EXISTS"] is False

    def test_backup(self):
        """Test backup creates backup file."""
        with tempfile.NamedTemporaryFile(delete=False) as f:
            f.write(b"test content")
            path = f.name

        try:
            backup_path = backup(path)
            assert backup_path is not None
            assert backup_path.exists()
            assert backup_path.suffix == ".bak"
        finally:
            os.unlink(path)
            if backup_path:
                os.unlink(backup_path)

    def test_backup_nonexistent(self):
        """Test backup returns None for non-existent file."""
        result = backup("/nonexistent/file")
        assert result is None

    def test_try_remove(self):
        """Test try_remove removes file."""
        with tempfile.NamedTemporaryFile(delete=False) as f:
            path = f.name

        assert Path(path).exists()
        result = try_remove(path)
        assert result is True
        assert not Path(path).exists()

    def test_try_remove_nonexistent(self):
        """Test try_remove returns True for non-existent."""
        result = try_remove("/nonexistent/file")
        assert result is True

    def test_temp_dir(self):
        """Test temp_dir returns valid temp directory."""
        td = temp_dir()
        assert td.exists()
        assert td.is_dir()


class TestListUtilities:
    """Tests for list utility functions."""

    def test_simplify_nested_list_simple(self):
        """Test simplify_nested_list with simple nesting."""
        result = simplify_nested_list([[1, 2], [3, 4]])
        assert result == [1, 2, 3, 4]

    def test_simplify_nested_list_deep(self):
        """Test simplify_nested_list with deep nesting."""
        nested = [[1, 2, 3], 4, [[5, 6], [[7], [8], [9, 10, 11]], [[12, 13], 14]]]
        result = simplify_nested_list(nested)
        assert result == [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]

    def test_simplify_nested_list_empty(self):
        """Test simplify_nested_list with empty list."""
        result = simplify_nested_list([])
        assert result == []

    def test_simplify_nested_list_flat(self):
        """Test simplify_nested_list with already flat list."""
        result = simplify_nested_list([1, 2, 3])
        assert result == [1, 2, 3]


class TestMaskUtilities:
    """Tests for mask utility functions."""

    def test_amber_mask_to_dict_simple(self):
        """Test amber_mask_to_dict with simple mask."""
        result = amber_mask_to_dict("ETA@O1,C1")
        assert result == {"ETA": ["O1", "C1"]}

    def test_amber_mask_to_dict_with_colon(self):
        """Test amber_mask_to_dict with leading colon."""
        result = amber_mask_to_dict(":ETA@O1,C1;WAT@O")
        assert result == {"ETA": ["O1", "C1"], "WAT": ["O"]}

    def test_amber_mask_to_dict_residue_only(self):
        """Test amber_mask_to_dict with residue only."""
        result = amber_mask_to_dict("ETA")
        assert result == {"ETA": ["all"]}

    def test_parse_num_mask_simple(self):
        """Test parse_num_mask with simple list."""
        result = parse_num_mask("1,2,3")
        assert result == [1, 2, 3]

    def test_parse_num_mask_with_range_colon(self):
        """Test parse_num_mask with colon range."""
        result = parse_num_mask("1,2,5:10")
        assert result == [1, 2, 5, 6, 7, 8, 9, 10]

    def test_parse_num_mask_with_range_dash(self):
        """Test parse_num_mask with dash range."""
        result = parse_num_mask("1-3,7,10-12")
        assert result == [1, 2, 3, 7, 10, 11, 12]

    def test_num_list_to_mask_simple(self):
        """Test num_list_to_mask with simple list."""
        result = num_list_to_mask([1, 2, 3])
        assert result == "1-3"

    def test_num_list_to_mask_with_gaps(self):
        """Test num_list_to_mask with gaps."""
        result = num_list_to_mask([1, 2, 3, 7, 10, 11, 12])
        assert result == "1-3,7,10-12"

    def test_num_list_to_mask_empty(self):
        """Test num_list_to_mask with empty list."""
        result = num_list_to_mask([])
        assert result == ""

    def test_num_list_to_mask_single(self):
        """Test num_list_to_mask with single element."""
        result = num_list_to_mask([5])
        assert result == "5"


class TestLoggingUtilities:
    """Tests for logging utility functions."""

    def test_log_formatter_instantiation(self):
        """Test LogFormatter can be instantiated."""
        formatter = LogFormatter()
        assert formatter is not None

    def test_clip_str_no_truncation(self):
        """Test clip_str without truncation."""
        result = clip_str("hello", 10)
        assert result == "hello"

    def test_clip_str_with_truncation(self):
        """Test clip_str with truncation."""
        result = clip_str("hello world", 8)
        assert result == "hello..."
        assert len(result) == 8

    def test_traceback_plus_no_exception(self):
        """Test traceback_plus with no active exception."""
        result = traceback_plus()
        assert "No exception" in result


class TestLegacyAliases:
    """Test that legacy camelCase aliases work."""

    def test_projectRoot(self):
        from pymdmix.utils import projectRoot

        assert projectRoot() == project_root()

    def test_dataRoot(self):
        from pymdmix.utils import dataRoot

        assert dataRoot() == data_root()

    def test_testRoot(self):
        from pymdmix.utils import testRoot

        assert testRoot() == test_root()

    def test_parseNumMask(self):
        from pymdmix.utils import parseNumMask

        assert parseNumMask("1,2,3") == parse_num_mask("1,2,3")

    def test_numListToMask(self):
        from pymdmix.utils import numListToMask

        assert numListToMask([1, 2, 3]) == num_list_to_mask([1, 2, 3])
