"""
Tests for Project Browser
=========================

Tests for pymdmix.project.browser module.
"""
import pytest
import os
from pathlib import Path
from dataclasses import dataclass

from pymdmix.project.browser import (
    Browser,
    BrowserError,
    DirectoryNotFoundError,
    find_project_root,
    list_projects,
)


@dataclass
class MockReplica:
    """Mock replica for testing."""
    path: Path


@dataclass
class MockProject:
    """Mock project for testing."""
    path: Path


class TestBrowser:
    """Tests for Browser class."""
    
    @pytest.fixture
    def test_dir(self, tmp_path):
        """Create a test directory structure."""
        # Create project structure
        (tmp_path / "MD").mkdir()
        (tmp_path / "MD" / "replica1").mkdir()
        (tmp_path / "MD" / "replica2").mkdir()
        (tmp_path / "analysis").mkdir()
        
        return tmp_path
    
    @pytest.fixture
    def browser(self, test_dir):
        """Create browser with test dir as home."""
        original_cwd = os.getcwd()
        os.chdir(test_dir)
        
        b = Browser()
        b.set_home(test_dir)
        
        yield b
        
        os.chdir(original_cwd)
    
    def test_create_browser(self):
        """Test creating a browser."""
        browser = Browser()
        
        assert browser.cwd == Path.cwd()
        assert browser.home == Path.cwd()
    
    def test_set_home(self, test_dir):
        """Test setting home directory."""
        original_cwd = os.getcwd()
        
        try:
            browser = Browser()
            browser.set_home(test_dir)
            
            assert browser.home == test_dir
            assert browser.cwd == test_dir
        finally:
            os.chdir(original_cwd)
    
    def test_set_home_nonexistent(self, tmp_path):
        """Test setting nonexistent home."""
        browser = Browser()
        
        with pytest.raises(DirectoryNotFoundError):
            browser.set_home(tmp_path / "nonexistent")
    
    def test_chdir(self, browser, test_dir):
        """Test changing directory."""
        browser.chdir(test_dir / "MD")
        
        assert browser.cwd == test_dir / "MD"
        assert browser.prev_dir == test_dir
    
    def test_chdir_nonexistent(self, browser, test_dir):
        """Test chdir to nonexistent directory."""
        with pytest.raises(DirectoryNotFoundError):
            browser.chdir(test_dir / "nonexistent")
    
    def test_go_home(self, browser, test_dir):
        """Test returning to home."""
        browser.chdir(test_dir / "MD")
        browser.go_home()
        
        assert browser.cwd == test_dir
    
    def test_go_back(self, browser, test_dir):
        """Test going back to previous directory."""
        browser.chdir(test_dir / "MD")
        browser.chdir(test_dir / "analysis")
        browser.go_back()
        
        assert browser.cwd == test_dir / "MD"
    
    def test_go_up(self, browser, test_dir):
        """Test going to parent directory."""
        browser.chdir(test_dir / "MD" / "replica1")
        browser.go_up()
        
        assert browser.cwd == test_dir / "MD"
    
    def test_go_md(self, browser, test_dir):
        """Test navigating to MD directory."""
        browser.go_md()
        
        assert browser.cwd == test_dir / "MD"
    
    def test_go_replica(self, browser, test_dir):
        """Test navigating to replica."""
        replica = MockReplica(path=test_dir / "MD" / "replica1")
        browser.go_replica(replica)
        
        assert browser.cwd == test_dir / "MD" / "replica1"
    
    def test_go_replica_with_subfolder(self, browser, test_dir):
        """Test navigating to replica subfolder."""
        # Create subfolder
        (test_dir / "MD" / "replica1" / "output").mkdir()
        
        replica = MockReplica(path=test_dir / "MD" / "replica1")
        browser.go_replica(replica, subfolder="output")
        
        assert browser.cwd == test_dir / "MD" / "replica1" / "output"
    
    def test_getcwd(self, browser, test_dir):
        """Test getcwd updates internal state."""
        os.chdir(test_dir / "MD")
        cwd = browser.getcwd()
        
        assert cwd == test_dir / "MD"
        assert browser.cwd == test_dir / "MD"
    
    def test_list_dir(self, browser, test_dir):
        """Test listing directory contents."""
        browser.go_md()
        contents = browser.list_dir()
        
        names = [p.name for p in contents]
        assert "replica1" in names
        assert "replica2" in names
    
    def test_list_dir_pattern(self, browser, test_dir):
        """Test listing with pattern."""
        browser.go_md()
        contents = browser.list_dir("replica*")
        
        assert len(contents) == 2
    
    def test_list_replicas(self, browser, test_dir):
        """Test listing replica directories."""
        replicas = browser.list_replicas()
        
        names = [p.name for p in replicas]
        assert "replica1" in names
        assert "replica2" in names
    
    def test_list_replicas_no_md(self, tmp_path):
        """Test listing replicas when no MD dir."""
        original_cwd = os.getcwd()
        
        try:
            browser = Browser()
            browser.set_home(tmp_path)
            
            replicas = browser.list_replicas()
            assert replicas == []
        finally:
            os.chdir(original_cwd)
    
    def test_set_project(self, test_dir):
        """Test setting project."""
        original_cwd = os.getcwd()
        
        try:
            project = MockProject(path=test_dir)
            browser = Browser()
            browser.set_project(project)
            
            assert browser.project == project
            assert browser.home == test_dir
        finally:
            os.chdir(original_cwd)
    
    def test_repr(self, browser, test_dir):
        """Test string representation."""
        repr_str = repr(browser)
        
        assert "Browser" in repr_str
        assert "cwd=" in repr_str
        assert "home=" in repr_str


class TestFindProjectRoot:
    """Tests for find_project_root function."""
    
    def test_find_with_config(self, tmp_path):
        """Test finding root via mdmix.cfg."""
        (tmp_path / "mdmix.cfg").touch()
        (tmp_path / "MD").mkdir()
        
        original_cwd = os.getcwd()
        os.chdir(tmp_path / "MD")
        
        try:
            root = find_project_root()
            assert root == tmp_path
        finally:
            os.chdir(original_cwd)
    
    def test_find_with_dotdir(self, tmp_path):
        """Test finding root via .mdmix directory."""
        (tmp_path / ".mdmix").mkdir()
        (tmp_path / "nested" / "deep").mkdir(parents=True)
        
        root = find_project_root(tmp_path / "nested" / "deep")
        assert root == tmp_path
    
    def test_not_found(self, tmp_path):
        """Test when no project root exists."""
        (tmp_path / "random").mkdir()
        
        root = find_project_root(tmp_path / "random")
        assert root is None
    
    def test_find_from_start(self, tmp_path):
        """Test finding from explicit start path."""
        (tmp_path / "project").mkdir(parents=True, exist_ok=True)
        (tmp_path / "project" / "mdmix.cfg").touch()
        
        root = find_project_root(tmp_path / "project")
        assert root == tmp_path / "project"


class TestListProjects:
    """Tests for list_projects function."""
    
    def test_list_projects(self, tmp_path):
        """Test listing projects."""
        # Create some projects
        (tmp_path / "proj1").mkdir(parents=True, exist_ok=True)
        (tmp_path / "proj1" / "mdmix.cfg").touch()
        (tmp_path / "proj2" / ".mdmix").mkdir(parents=True)
        (tmp_path / "not_a_project").mkdir()
        
        projects = list_projects(tmp_path)
        
        names = [p.name for p in projects]
        assert "proj1" in names
        assert "proj2" in names
        assert "not_a_project" not in names
    
    def test_list_empty(self, tmp_path):
        """Test listing with no projects."""
        projects = list_projects(tmp_path)
        assert projects == []
    
    def test_list_sorted(self, tmp_path):
        """Test projects are sorted."""
        (tmp_path / "z_proj").mkdir(parents=True, exist_ok=True)
        (tmp_path / "z_proj" / "mdmix.cfg").touch()
        (tmp_path / "a_proj").mkdir(parents=True, exist_ok=True)
        (tmp_path / "a_proj" / "mdmix.cfg").touch()
        
        projects = list_projects(tmp_path)
        
        assert projects[0].name == "a_proj"
        assert projects[1].name == "z_proj"
