"""
Shared pytest fixtures for pyMDMix tests.
"""

import os
import shutil
import tempfile
from pathlib import Path

import numpy as np
import pytest

# ---------------------------------------------------------------------------
# Test data directory
# ---------------------------------------------------------------------------

# All test data lives under tests/data/ — NOT inside the pymdmix package.
_TEST_DATA_DIR = Path(__file__).parent / "data"
_PEP_DATA_DIR = _TEST_DATA_DIR / "pep"
_AMBER_DATA_DIR = _TEST_DATA_DIR / "amber"
_GRIDS_DATA_DIR = _TEST_DATA_DIR / "grids"


# ---------------------------------------------------------------------------
# Custom marker auto-skip hooks
# ---------------------------------------------------------------------------


def pytest_configure(config):
    config.addinivalue_line(
        "markers",
        "ambertools: marks tests that require AmberTools (tleap, LEaP) to be installed",
    )
    config.addinivalue_line(
        "markers",
        "cpptraj: marks tests that require the cpptraj binary from AmberTools",
    )
    config.addinivalue_line(
        "markers",
        "real_data: marks tests that require large test data files (traj.nc, DX grids, etc.) in tests/data/",
    )


def _tleap_available() -> bool:
    """Return True if tleap can be found via AMBERHOME or PATH."""
    amber_home = os.environ.get("AMBERHOME")
    if amber_home:
        for sub in ("bin", "exe"):
            candidate = Path(amber_home) / sub / "tleap"
            if candidate.exists():
                return True
    return shutil.which("tleap") is not None


def _cpptraj_available() -> bool:
    """Return True if cpptraj can be found via AMBER_PTRAJ or PATH."""
    return bool(os.environ.get("AMBER_PTRAJ")) or shutil.which("cpptraj") is not None


def _real_data_available() -> bool:
    """Return True if the test trajectory file is present (LFS objects checked out)."""
    return (_AMBER_DATA_DIR / "traj.nc").exists()


def pytest_collection_modifyitems(config, items):
    """Auto-skip ambertools/cpptraj/real_data tests when the required tools or data are absent."""
    amber_ok = _tleap_available()
    cpptraj_ok = _cpptraj_available()
    real_data_ok = _real_data_available()

    skip_amber = pytest.mark.skip(
        reason=(
            "AmberTools not available. "
            "Set AMBERHOME or run inside Docker: ./scripts/run_ambertools_tests.sh"
        )
    )
    skip_cpptraj = pytest.mark.skip(
        reason=(
            "cpptraj not available. "
            "Set AMBER_PTRAJ or add cpptraj to PATH; run inside Docker: "
            "./scripts/run_ambertools_tests.sh"
        )
    )
    skip_real_data = pytest.mark.skip(
        reason=(
            "Real test data not available (Git LFS objects not checked out). "
            "Run: git lfs pull"
        )
    )

    for item in items:
        if "ambertools" in item.keywords and not amber_ok:
            item.add_marker(skip_amber)
        if "cpptraj" in item.keywords and not cpptraj_ok:
            item.add_marker(skip_cpptraj)
        if "real_data" in item.keywords and not real_data_ok:
            item.add_marker(skip_real_data)


# ---------------------------------------------------------------------------
# Test data path fixtures (tests/data/)
# ---------------------------------------------------------------------------


@pytest.fixture(scope="session")
def test_data_dir() -> Path:
    """Root of the tests/data/ directory."""
    return _TEST_DATA_DIR


@pytest.fixture(scope="session")
def pep_data_dir() -> Path:
    """Path to tests/data/pep/ containing the bundled test peptide."""
    return _PEP_DATA_DIR


@pytest.fixture(scope="session")
def pep_pdb_path() -> Path:
    """Path to the test peptide PDB file (tests/data/pep/pep.pdb)."""
    p = _PEP_DATA_DIR / "pep.pdb"
    assert p.exists(), f"Test PDB not found: {p}"
    return p


@pytest.fixture(scope="session")
def pep_prmtop_path() -> Path:
    """Path to the test peptide Amber topology (tests/data/pep/pep.prmtop)."""
    p = _PEP_DATA_DIR / "pep.prmtop"
    assert p.exists(), f"Test prmtop not found: {p}"
    return p


@pytest.fixture(scope="session")
def pep_prmcrd_path() -> Path:
    """Path to the test peptide Amber coordinates (tests/data/pep/pep.prmcrd)."""
    p = _PEP_DATA_DIR / "pep.prmcrd"
    assert p.exists(), f"Test prmcrd not found: {p}"
    return p


@pytest.fixture(scope="session")
def pep_off_path() -> Path:
    """Path to the test peptide LEaP object file (tests/data/pep/pep.off)."""
    p = _PEP_DATA_DIR / "pep.off"
    assert p.exists(), f"Test OFF not found: {p}"
    return p


@pytest.fixture(scope="session")
def solvated_pep_pdb_path() -> Path:
    """
    Path to the pre-solvated peptide structure (tests/data/amber/pep_WAT_WAT_1.pdb).

    This is a WAT-solvated system from the legacy test data.
    Marked ``real_data`` — skip if file absent (LFS not checked out).
    """
    return _AMBER_DATA_DIR / "pep_WAT_WAT_1.pdb"


@pytest.fixture(scope="session")
def amber_traj_nc_path() -> Path:
    """
    Path to the Amber NetCDF trajectory (tests/data/amber/traj.nc).

    This is a multi-frame WAT-solvated trajectory of the test peptide from the
    legacy test data.  Marked ``real_data`` — skip if file absent (LFS not
    checked out).
    """
    return _AMBER_DATA_DIR / "traj.nc"


@pytest.fixture(scope="session")
def eta_ct_dx_path() -> Path:
    """
    Path to the pre-computed ETA CT probe density grid (tests/data/grids/ETA_CT.dx).

    This is a full-size DX density grid computed from the legacy test trajectory.
    Marked ``real_data`` — skip if file absent (LFS not checked out).
    """
    return _GRIDS_DATA_DIR / "ETA_CT.dx"


# ---------------------------------------------------------------------------
# General-purpose fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def tmp_output_dir():
    """Temporary directory for test outputs."""
    path = Path(tempfile.mkdtemp(prefix="pymdmix_test_"))
    yield path
    shutil.rmtree(path, ignore_errors=True)


@pytest.fixture
def sample_coordinates():
    """Sample atomic coordinates for testing."""
    # Small protein-like coordinates (100 atoms)
    np.random.seed(42)
    coords = np.random.randn(100, 3) * 10  # ~10 Angstrom spread
    coords += np.array([25, 25, 25])  # Center around (25, 25, 25)
    return coords.astype(np.float64)


@pytest.fixture
def sample_trajectory_coords():
    """Sample trajectory coordinates (10 frames, 100 atoms)."""
    np.random.seed(42)
    n_frames = 10
    n_atoms = 100

    # Base coordinates
    base = np.random.randn(n_atoms, 3) * 10
    base += np.array([25, 25, 25])

    # Add small perturbations for each frame
    frames = []
    for i in range(n_frames):
        perturbation = np.random.randn(n_atoms, 3) * 0.5
        frames.append((base + perturbation).astype(np.float64))

    return frames


@pytest.fixture
def sample_pdb_content():
    """Sample PDB file content."""
    return """ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00  0.00           C
ATOM      3  C   ALA A   1       2.009   1.420   0.000  1.00  0.00           C
ATOM      4  O   ALA A   1       1.246   2.390   0.000  1.00  0.00           O
ATOM      5  CB  ALA A   1       1.986  -0.760  -1.216  1.00  0.00           C
ATOM      6  N   GLY A   2       3.326   1.544   0.000  1.00  0.00           N
ATOM      7  CA  GLY A   2       3.941   2.861   0.000  1.00  0.00           C
ATOM      8  C   GLY A   2       5.457   2.779   0.000  1.00  0.00           C
ATOM      9  O   GLY A   2       6.030   1.689   0.000  1.00  0.00           O
ATOM     10  OXT GLY A   2       6.073   3.848   0.000  1.00  0.00           O
TER
HETATM   11  O   WAT A   3      10.000  10.000  10.000  1.00  0.00           O
HETATM   12  H1  WAT A   3      10.757  10.000  10.586  1.00  0.00           H
HETATM   13  H2  WAT A   3       9.243  10.000  10.586  1.00  0.00           H
END
"""


@pytest.fixture
def sample_pdb_file(tmp_output_dir, sample_pdb_content):
    """Sample PDB file."""
    pdb_path = tmp_output_dir / "sample.pdb"
    pdb_path.write_text(sample_pdb_content)
    return pdb_path


@pytest.fixture
def sample_dx_content():
    """Sample DX file content (3x3x3 grid)."""
    return """# OpenDX density grid
# Generated by pyMDMix
object 1 class gridpositions counts 3 3 3
origin 0.000000 0.000000 0.000000
delta 1.000000 0.000000 0.000000
delta 0.000000 1.000000 0.000000
delta 0.000000 0.000000 1.000000
object 2 class gridconnections counts 3 3 3
object 3 class array type double rank 0 items 27 data follows
1.000000e+00 2.000000e+00 3.000000e+00
4.000000e+00 5.000000e+00 6.000000e+00
7.000000e+00 8.000000e+00 9.000000e+00
1.000000e+01 1.100000e+01 1.200000e+01
1.300000e+01 1.400000e+01 1.500000e+01
1.600000e+01 1.700000e+01 1.800000e+01
1.900000e+01 2.000000e+01 2.100000e+01
2.200000e+01 2.300000e+01 2.400000e+01
2.500000e+01 2.600000e+01 2.700000e+01
attribute "dep" string "positions"
object "density" class field
component "positions" value 1
component "connections" value 2
component "data" value 3
"""


@pytest.fixture
def sample_dx_file(tmp_output_dir, sample_dx_content):
    """Sample DX file."""
    dx_path = tmp_output_dir / "sample.dx"
    dx_path.write_text(sample_dx_content)
    return dx_path


class MockTrajectoryReader:
    """Mock trajectory reader for testing without real files."""

    def __init__(self, frames: list[np.ndarray]):
        self._frames = frames
        self._n_atoms = frames[0].shape[0] if frames else 0

    @property
    def n_frames(self) -> int:
        return len(self._frames)

    @property
    def n_atoms(self) -> int:
        return self._n_atoms

    def __len__(self) -> int:
        return self.n_frames

    def __iter__(self):
        from pymdmix.core.trajectory import Frame

        for coords in self._frames:
            yield Frame(coordinates=coords)


@pytest.fixture
def mock_trajectory(sample_trajectory_coords):
    """Mock trajectory reader with sample coordinates."""
    return MockTrajectoryReader(sample_trajectory_coords)
