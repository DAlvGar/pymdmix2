"""
Tests for DCD Parser
====================

Tests for pymdmix.io.dcd_parser module.
"""

import struct
from pathlib import Path

import numpy as np
import pytest

from pymdmix.io.dcd_parser import (
    DCDFormatError,
    DCDFrame,
    DCDHeader,
    DCDReader,
    read_dcd,
)


def create_test_dcd(path: Path, n_frames: int = 5, n_atoms: int = 10, has_box: bool = True) -> None:
    """Create a minimal valid DCD file for testing."""
    with open(path, "wb") as f:
        # Header block (92 bytes total)
        # Format: >I 4s 9I f 11I = 4 + 4 + 36 + 4 + 44 = 92 bytes, 23 items
        # Items 0-10: blocksize, magic, NSET, ISTART, NSAVC, NTOT, 5 unused (11 items)
        # Item 11: DELTA (float)
        # Items 12-22: extrablock, 4dims, 8 unused, blocksize (11 items)
        f.write(struct.pack(">I", 84))  # [0] Block size
        f.write(b"CORD")  # [1] Magic
        f.write(struct.pack(">I", n_frames))  # [2] NSET
        f.write(struct.pack(">I", 0))  # [3] ISTART
        f.write(struct.pack(">I", 1))  # [4] NSAVC
        f.write(struct.pack(">I", n_frames))  # [5] NTOT
        f.write(struct.pack(">5I", 0, 0, 0, 0, 0))  # [6-10] 5 unused
        f.write(struct.pack(">f", 0.002))  # [11] DELTA
        f.write(struct.pack(">I", 1 if has_box else 0))  # [12] Extra block flag
        f.write(struct.pack(">I", 1))  # [13] Charmm version (non-zero!)
        f.write(struct.pack(">8I", 0, 0, 0, 0, 0, 0, 0, 0))  # [14-21] 8 unused
        f.write(struct.pack(">I", 84))  # [22] Block size

        # Title block
        # Single title line
        f.write(struct.pack(">I", 84))  # Block size = 4 + 80
        f.write(struct.pack(">I", 1))  # Number of title lines
        f.write(struct.pack(">80s", b"Test DCD file"))
        f.write(struct.pack(">I", 84))  # Closing block

        # Atom count block
        f.write(struct.pack(">3I", 4, n_atoms, 4))

        # Frames
        for frame_idx in range(n_frames):
            # Extra block (unit cell) if present
            if has_box:
                f.write(struct.pack(">I", 48))  # Block size
                # Unit cell: a, gamma, b, beta, alpha, c
                unitcell = np.array([50.0, 90.0, 50.0, 90.0, 90.0, 50.0], dtype=">f8")
                f.write(unitcell.tobytes())
                f.write(struct.pack(">I", 48))

            # Coordinates (x, y, z stored separately)
            float_block_size = 4 * n_atoms

            # X coordinates
            f.write(struct.pack(">I", float_block_size))
            x = np.arange(n_atoms, dtype=">f4") + frame_idx
            f.write(x.tobytes())
            f.write(struct.pack(">I", float_block_size))

            # Y coordinates
            f.write(struct.pack(">I", float_block_size))
            y = np.arange(n_atoms, dtype=">f4") + frame_idx + 10
            f.write(y.tobytes())
            f.write(struct.pack(">I", float_block_size))

            # Z coordinates
            f.write(struct.pack(">I", float_block_size))
            z = np.arange(n_atoms, dtype=">f4") + frame_idx + 20
            f.write(z.tobytes())
            f.write(struct.pack(">I", float_block_size))


class TestDCDHeader:
    """Tests for DCDHeader dataclass."""

    def test_header_creation(self):
        """Test creating a DCDHeader."""
        header = DCDHeader(
            n_frames=100,
            istart=0,
            nsavc=500,
            n_total_steps=50000,
            delta=0.002,
            has_extra_block=True,
            has_4dims=False,
            n_atoms=1000,
            title=["Test title"],
        )

        assert header.n_frames == 100
        assert header.n_atoms == 1000
        assert header.has_extra_block is True

    def test_header_size(self):
        """Test header size calculation."""
        header = DCDHeader(
            n_frames=10,
            istart=0,
            nsavc=1,
            n_total_steps=10,
            delta=0.002,
            has_extra_block=False,
            has_4dims=False,
            n_atoms=100,
            title=["Title 1"],  # 1 title line
        )

        # 116 + 80 * 1 = 196
        assert header.header_size == 196

    def test_frame_size(self):
        """Test frame size calculation."""
        header = DCDHeader(
            n_frames=10,
            istart=0,
            nsavc=1,
            n_total_steps=10,
            delta=0.002,
            has_extra_block=False,
            has_4dims=False,
            n_atoms=100,
            title=["Title"],
        )

        # 3 * 4 * 100 + 24 = 1224
        assert header.frame_size == 1224

    def test_frame_size_with_box(self):
        """Test frame size with extra block."""
        header = DCDHeader(
            n_frames=10,
            istart=0,
            nsavc=1,
            n_total_steps=10,
            delta=0.002,
            has_extra_block=True,
            has_4dims=False,
            n_atoms=100,
            title=["Title"],
        )

        # 3 * 4 * 100 + 24 + 56 = 1280
        assert header.frame_size == 1280


class TestDCDFrame:
    """Tests for DCDFrame dataclass."""

    def test_frame_creation(self):
        """Test creating a DCDFrame."""
        coords = np.random.rand(100, 3).astype(np.float32)
        frame = DCDFrame(coordinates=coords)

        assert frame.n_atoms == 100
        assert frame.coordinates.shape == (100, 3)
        assert frame.unitcell is None

    def test_frame_with_unitcell(self):
        """Test frame with unit cell info."""
        coords = np.random.rand(50, 3).astype(np.float32)
        unitcell = np.array([50.0, 90.0, 50.0, 90.0, 90.0, 50.0], dtype=np.float64)
        frame = DCDFrame(coordinates=coords, unitcell=unitcell)

        assert frame.n_atoms == 50
        assert frame.unitcell is not None
        assert len(frame.unitcell) == 6


class TestDCDReader:
    """Tests for DCDReader class."""

    @pytest.fixture
    def dcd_file(self, tmp_path):
        """Create a temporary DCD file."""
        path = tmp_path / "test.dcd"
        create_test_dcd(path, n_frames=5, n_atoms=10, has_box=True)
        return path

    @pytest.fixture
    def dcd_no_box(self, tmp_path):
        """Create a DCD file without box info."""
        path = tmp_path / "test_nobox.dcd"
        create_test_dcd(path, n_frames=3, n_atoms=20, has_box=False)
        return path

    def test_open_dcd(self, dcd_file):
        """Test opening a DCD file."""
        with DCDReader(dcd_file) as reader:
            assert reader.n_frames == 5
            assert reader.n_atoms == 10

    def test_header_properties(self, dcd_file):
        """Test header property access."""
        with DCDReader(dcd_file) as reader:
            header = reader.header
            assert header.n_frames == 5
            assert header.n_atoms == 10
            assert header.has_extra_block is True

    def test_iteration(self, dcd_file):
        """Test iterating through frames."""
        frames = []
        with DCDReader(dcd_file) as reader:
            for frame in reader:
                frames.append(frame)

        assert len(frames) == 5
        for frame in frames:
            assert frame.n_atoms == 10
            assert frame.unitcell is not None

    def test_random_access(self, dcd_file):
        """Test random frame access."""
        with DCDReader(dcd_file) as reader:
            frame0 = reader[0]
            frame2 = reader[2]
            frame4 = reader[4]

            assert frame0.n_atoms == 10
            # Coordinates increase with frame index
            assert frame2.coordinates[0, 0] > frame0.coordinates[0, 0]
            assert frame4.coordinates[0, 0] > frame2.coordinates[0, 0]

    def test_read_all(self, dcd_file):
        """Test reading all frames at once."""
        with DCDReader(dcd_file) as reader:
            coords = reader.read_all()

        assert coords.shape == (5, 10, 3)

    def test_len(self, dcd_file):
        """Test length method."""
        with DCDReader(dcd_file) as reader:
            assert len(reader) == 5

    def test_without_box(self, dcd_no_box):
        """Test reading DCD without box info."""
        with DCDReader(dcd_no_box) as reader:
            assert reader.n_frames == 3
            assert reader.n_atoms == 20

            for frame in reader:
                assert frame.unitcell is None

    def test_repr(self, dcd_file):
        """Test string representation."""
        with DCDReader(dcd_file) as reader:
            repr_str = repr(reader)
            assert "DCDReader" in repr_str
            assert "n_frames=5" in repr_str
            assert "n_atoms=10" in repr_str

    def test_file_not_found(self, tmp_path):
        """Test error on missing file."""
        with pytest.raises(FileNotFoundError):
            DCDReader(tmp_path / "nonexistent.dcd")

    def test_invalid_format(self, tmp_path):
        """Test error on invalid format."""
        bad_file = tmp_path / "bad.dcd"
        bad_file.write_bytes(b"This is not a DCD file")

        with pytest.raises(DCDFormatError):
            DCDReader(bad_file)

    def test_index_out_of_range(self, dcd_file):
        """Test index bounds checking."""
        with DCDReader(dcd_file) as reader:
            with pytest.raises(IndexError):
                reader[10]
            with pytest.raises(IndexError):
                reader[-1]


class TestReadDCDFunction:
    """Tests for read_dcd convenience function."""

    def test_read_dcd(self, tmp_path):
        """Test the convenience function."""
        dcd_path = tmp_path / "test.dcd"
        create_test_dcd(dcd_path, n_frames=3, n_atoms=5)

        coords = read_dcd(dcd_path)

        assert coords.shape == (3, 5, 3)
