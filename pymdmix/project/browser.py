"""
Project Browser
===============

Directory navigation and browsing utilities for MDMix projects.
Tracks current directory state and provides navigation shortcuts.

Examples
--------
>>> from pymdmix.project.browser import Browser
>>> browser = Browser()
>>> browser.set_home("/path/to/project")
>>> browser.go_md()  # Navigate to MD folder
>>> browser.go_replica(replica)  # Navigate to replica
>>> browser.go_back()  # Return to previous directory
"""

from __future__ import annotations

import logging
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from pymdmix.project.project import Project
    from pymdmix.project.replica import Replica

log = logging.getLogger(__name__)


class BrowserError(Exception):
    """Base exception for browser errors."""

    pass


class DirectoryNotFoundError(BrowserError):
    """Directory does not exist."""

    pass


@dataclass
class Browser:
    """
    Project directory browser.

    Wraps os directory functions to work with MDMix projects.
    Tracks navigation history and provides shortcuts for
    common project locations.

    Parameters
    ----------
    project : Project, optional
        MDMix project to browse

    Attributes
    ----------
    cwd : Path
        Current working directory
    home : Path
        Project home directory
    prev_dir : Path
        Previous directory (for go_back)

    Examples
    --------
    >>> browser = Browser()
    >>> browser.set_home("/path/to/project")
    >>>
    >>> # Navigate around
    >>> browser.chdir("MD/replica1")
    >>> browser.go_home()
    >>> browser.go_back()  # Returns to MD/replica1
    """

    project: Project | None = None
    cwd: Path = field(default_factory=Path.cwd)
    home: Path = field(default_factory=Path.cwd)
    prev_dir: Path = field(default_factory=Path.cwd)

    def __post_init__(self) -> None:
        """Initialize with project home if provided."""
        if self.project is not None:
            self.set_home(self.project.path)

    def set_project(self, project: Project) -> None:
        """
        Set the project to browse.

        Parameters
        ----------
        project : Project
            MDMix project
        """
        self.project = project
        self.set_home(project.path)

    def set_home(self, home_dir: str | Path | None = None) -> Path:
        """
        Set the project home directory.

        Parameters
        ----------
        home_dir : str or Path, optional
            Home directory path. If None, uses current directory.

        Returns
        -------
        Path
            The home directory
        """
        if home_dir is not None:
            home_path = Path(home_dir)
            if not home_path.exists():
                raise DirectoryNotFoundError(f"Directory does not exist: {home_path}")
            os.chdir(home_path)

        self.home = self._update()
        return self.home

    def _update(self) -> Path:
        """Update directory state after navigation."""
        self.prev_dir = self.cwd
        self.cwd = Path.cwd()
        log.debug(f"Directory change: {self.prev_dir} -> {self.cwd}")
        return self.cwd

    def chdir(self, new_dir: str | Path) -> Path:
        """
        Change to a new directory.

        Parameters
        ----------
        new_dir : str or Path
            Target directory

        Returns
        -------
        Path
            New current directory
        """
        new_path = Path(new_dir)
        if not new_path.exists():
            raise DirectoryNotFoundError(f"Directory does not exist: {new_path}")

        os.chdir(new_path)
        return self._update()

    def go_home(self) -> Path:
        """
        Navigate to project home directory.

        Returns
        -------
        Path
            Home directory path
        """
        return self.chdir(self.home)

    def go_back(self) -> Path:
        """
        Navigate to previous directory (like cd -).

        Returns
        -------
        Path
            Previous directory path
        """
        return self.chdir(self.prev_dir)

    def go_up(self) -> Path:
        """
        Navigate to parent directory (like cd ..).

        Returns
        -------
        Path
            Parent directory path
        """
        return self.chdir(self.cwd.parent)

    def go_md(self) -> Path:
        """
        Navigate to project MD directory.

        Returns
        -------
        Path
            MD directory path
        """
        self.go_home()
        return self.chdir("MD")

    def go_replica(
        self,
        replica: Replica,
        subfolder: str = "",
    ) -> Path:
        """
        Navigate to a replica directory.

        Parameters
        ----------
        replica : Replica
            Replica to navigate to
        subfolder : str, optional
            Subfolder within replica directory

        Returns
        -------
        Path
            Replica directory path
        """
        path = replica.path
        if subfolder:
            path = path / subfolder
        return self.chdir(path)

    def getcwd(self) -> Path:
        """
        Get current working directory.

        Also updates internal state if directory changed externally.

        Returns
        -------
        Path
            Current working directory
        """
        actual_cwd = Path.cwd()
        if self.cwd != actual_cwd:
            self._update()
        return actual_cwd

    def list_dir(self, pattern: str = "*") -> list[Path]:
        """
        List contents of current directory.

        Parameters
        ----------
        pattern : str
            Glob pattern for filtering

        Returns
        -------
        list[Path]
            Matching paths
        """
        return sorted(self.cwd.glob(pattern))

    def list_replicas(self) -> list[Path]:
        """
        List replica directories in MD folder.

        Returns
        -------
        list[Path]
            Paths to replica directories
        """
        md_dir = self.home / "MD"
        if not md_dir.exists():
            return []

        return sorted([p for p in md_dir.iterdir() if p.is_dir() and not p.name.startswith(".")])

    def __repr__(self) -> str:
        return f"Browser(cwd='{self.cwd}', home='{self.home}')"


def find_project_root(start: str | Path | None = None) -> Path | None:
    """
    Find MDMix project root by searching for markers.

    Searches upward from start directory for mdmix.cfg or
    .mdmix directory.

    Parameters
    ----------
    start : str or Path, optional
        Starting directory (default: current directory)

    Returns
    -------
    Path or None
        Project root path, or None if not found
    """
    current = Path(start) if start else Path.cwd()

    # Search upward
    for parent in [current] + list(current.parents):
        if (parent / "mdmix.cfg").exists():
            return parent
        if (parent / ".mdmix").is_dir():
            return parent

    return None


def list_projects(root: str | Path | None = None) -> list[Path]:
    """
    List MDMix projects in a directory.

    Parameters
    ----------
    root : str or Path, optional
        Root directory to search (default: current directory)

    Returns
    -------
    list[Path]
        Paths to project directories
    """
    root_path = Path(root) if root else Path.cwd()
    projects = []

    for path in root_path.iterdir():
        if not path.is_dir():
            continue
        if (path / "mdmix.cfg").exists() or (path / ".mdmix").is_dir():
            projects.append(path)

    return sorted(projects)
