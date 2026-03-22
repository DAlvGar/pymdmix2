"""
Tests for pymdmix.core.solvent module.
"""

from pymdmix.core.solvent import (
    Probe,
    Solvent,
    SolventLibrary,
    SolventResidue,
    create_standard_solvents,
)


class TestProbe:
    """Test Probe dataclass."""

    def test_probe_creation(self):
        """Test creating a probe."""
        probe = Probe(
            name="OH",
            residue="ETA",
            atoms=["O1"],
            description="Ethanol hydroxyl",
        )

        assert probe.name == "OH"
        assert probe.residue == "ETA"
        assert probe.atoms == ["O1"]

    def test_probe_with_types(self):
        """Test creating a probe with chemical types."""
        probe = Probe(
            name="OH",
            residue="ETA",
            atoms=["O1"],
            description="Ethanol hydroxyl",
            probe_types=["Don", "Acc"],
            probability=0.000266,
        )

        assert probe.probe_types == ["Don", "Acc"]
        assert probe.probability == 0.000266

    def test_probe_selection(self):
        """Test probe MDAnalysis selection string."""
        probe = Probe(name="OH", residue="ETA", atoms=["O1", "H1"])

        selection = probe.selection
        assert "resname ETA" in selection
        assert "name O1 H1" in selection

    def test_probe_repr(self):
        """Test probe string representation."""
        probe = Probe(name="OH", residue="ETA", atoms=["O1"])

        repr_str = repr(probe)
        assert "Probe" in repr_str
        assert "OH" in repr_str


class TestSolventResidue:
    """Test SolventResidue dataclass."""

    def test_residue_creation(self):
        """Test creating a solvent residue."""
        res = SolventResidue(name="ETA", count=200, description="Ethanol")

        assert res.name == "ETA"
        assert res.count == 200

    def test_residue_defaults(self):
        """Test default values."""
        res = SolventResidue(name="WAT")

        assert res.count == 0
        assert res.description == ""


class TestSolvent:
    """Test Solvent dataclass."""

    def test_solvent_creation(self):
        """Test creating a solvent."""
        solvent = Solvent(
            name="ETA",
            description="20% Ethanol",
            residues=[
                SolventResidue("ETA", 200),
                SolventResidue("WAT", 800),
            ],
            probes=[
                Probe("OH", "ETA", ["O1"]),
            ],
        )

        assert solvent.name == "ETA"
        assert len(solvent.residues) == 2
        assert len(solvent.probes) == 1

    def test_solvent_with_extra_fields(self):
        """Test creating a solvent with volume, is_ionic, etc."""
        solvent = Solvent(
            name="ETA",
            description="20% Ethanol",
            residues=[
                SolventResidue("ETA", 200),
                SolventResidue("WAT", 800),
            ],
            probes=[
                Probe("OH", "ETA", ["O1"], "Hydroxyl", ["Don", "Acc"], 0.000266),
                Probe("CT", "ETA", ["C1"], "Hydrophobic", ["Hyd"], 0.000266),
            ],
            volume=7988.43,
            is_ionic=False,
            total_charge=0.0,
            corrections={},
        )

        assert solvent.volume == 7988.43
        assert solvent.is_ionic is False
        assert solvent.total_charge == 0.0
        assert solvent.corrections == {}

    def test_get_types_map(self):
        """Test getting types map from solvent."""
        solvent = Solvent(
            name="ETA",
            probes=[
                Probe("OH", "ETA", ["O1"], "Hydroxyl", ["Don", "Acc"]),
                Probe("CT", "ETA", ["C1"], "Hydrophobic", ["Hyd"]),
                Probe("WAT", "WAT", ["O"], "Water", ["Wat"]),
            ],
        )

        types_map = solvent.get_types_map()

        assert types_map["OH"] == "Don,Acc"
        assert types_map["CT"] == "Hyd"
        assert types_map["WAT"] == "Wat"

    def test_calculate_probability(self):
        """Test probability calculation."""
        solvent = Solvent(
            name="ETA",
            probes=[
                Probe("OH", "ETA", ["O1"], probability=0.000266),
            ],
            volume=7988.43,
        )

        # Should return stored probability
        prob = solvent.calculate_probability("OH")
        assert prob == 0.000266

        # Missing probe should return None
        assert solvent.calculate_probability("MISSING") is None

    def test_get_probe(self):
        """Test getting probe by name."""
        solvent = Solvent(
            name="ETA",
            probes=[
                Probe("OH", "ETA", ["O1"]),
                Probe("CT", "ETA", ["C2"]),
            ],
        )

        probe = solvent.get_probe("OH")
        assert probe is not None
        assert probe.name == "OH"

        missing = solvent.get_probe("nonexistent")
        assert missing is None

    def test_get_residue_names(self):
        """Test getting residue names."""
        solvent = Solvent(
            name="ETA",
            residues=[
                SolventResidue("ETA", 200),
                SolventResidue("WAT", 800),
            ],
        )

        names = solvent.get_residue_names()
        assert names == ["ETA", "WAT"]

    def test_solvent_to_dict(self):
        """Test converting solvent to dictionary."""
        solvent = Solvent(
            name="ETA",
            description="20% Ethanol",
            residues=[SolventResidue("ETA", 200)],
            probes=[Probe("OH", "ETA", ["O1"])],
        )

        data = solvent.to_dict()

        assert data["name"] == "ETA"
        assert len(data["residues"]) == 1
        assert len(data["probes"]) == 1

    def test_solvent_from_dict(self):
        """Test creating solvent from dictionary."""
        data = {
            "name": "ETA",
            "description": "20% Ethanol",
            "residues": [{"name": "ETA", "count": 200}],
            "probes": [{"name": "OH", "residue": "ETA", "atoms": ["O1"]}],
        }

        solvent = Solvent.from_dict(data)

        assert solvent.name == "ETA"
        assert len(solvent.residues) == 1
        assert len(solvent.probes) == 1

    def test_solvent_json_roundtrip(self, tmp_output_dir):
        """Test saving and loading solvent from JSON."""
        solvent = Solvent(
            name="ETA",
            description="20% Ethanol",
            residues=[SolventResidue("ETA", 200)],
            probes=[Probe("OH", "ETA", ["O1"])],
        )

        json_path = tmp_output_dir / "eta.json"
        solvent.to_json(json_path)

        loaded = Solvent.from_json(json_path)

        assert loaded.name == solvent.name
        assert loaded.description == solvent.description
        assert len(loaded.residues) == len(solvent.residues)
        assert len(loaded.probes) == len(solvent.probes)

    def test_write_off(self, tmp_output_dir):
        """write_off() copies the configured OFF file content."""
        source_off = tmp_output_dir / "source.off"
        source_off.write_text("!mock off content\n")

        solvent = Solvent(name="ETA", off_file=source_off)
        output_off = tmp_output_dir / "output.off"

        chars_written = solvent.write_off(output_off)

        assert chars_written == len("!mock off content\n")
        assert output_off.read_text() == "!mock off content\n"

    def test_from_file_json(self, tmp_output_dir):
        """from_file() with .json extension delegates to from_json()."""
        solvent = Solvent(
            name="ETA",
            description="Ethanol",
            residues=[SolventResidue("ETA", 200)],
            probes=[Probe("OH", "ETA", ["O1"])],
        )
        json_path = tmp_output_dir / "eta.json"
        solvent.to_json(json_path)

        loaded = Solvent.from_file(json_path)

        assert loaded.name == "ETA"
        assert len(loaded.residues) == 1
        assert len(loaded.probes) == 1

    def test_from_file_fallback_json(self, tmp_output_dir):
        """from_file() with unknown extension falls back to JSON parsing."""
        solvent = Solvent(name="ETA")
        # Save as .json then rename to .dat — should still parse as JSON
        json_path = tmp_output_dir / "eta.json"
        solvent.to_json(json_path)
        dat_path = tmp_output_dir / "eta.dat"
        json_path.rename(dat_path)

        loaded = Solvent.from_file(dat_path)

        assert loaded.name == "ETA"

    def test_from_file_missing_raises(self):
        """from_file() raises FileNotFoundError for missing path."""
        import pytest

        with pytest.raises(FileNotFoundError):
            Solvent.from_file("/nonexistent/solvent.json")


class TestSolventLibrary:
    """Test SolventLibrary class."""

    def test_library_creation(self, tmp_output_dir):
        """Test creating a library."""
        library = SolventLibrary(library_path=tmp_output_dir)

        assert library.library_path == tmp_output_dir
        assert len(library) == 0  # Empty directory

    def test_library_add_and_get(self, tmp_output_dir):
        """Test adding and retrieving solvents."""
        library = SolventLibrary(library_path=tmp_output_dir)

        solvent = Solvent(name="TEST", description="Test solvent")
        library.add(solvent)

        retrieved = library.get("TEST")
        assert retrieved is not None
        assert retrieved.name == "TEST"

    def test_library_list_solvents(self, tmp_output_dir):
        """Test listing solvents."""
        library = SolventLibrary(library_path=tmp_output_dir)

        library.add(Solvent(name="AAA"))
        library.add(Solvent(name="ZZZ"))

        names = library.list_solvents()

        assert names == ["AAA", "ZZZ"]  # Should be sorted

    def test_library_save_all(self, tmp_output_dir):
        """Test saving all solvents."""
        library = SolventLibrary(library_path=tmp_output_dir)

        library.add(Solvent(name="TEST1"))
        library.add(Solvent(name="TEST2"))
        library.save_all()

        # Check files exist
        assert (tmp_output_dir / "test1.json").exists()
        assert (tmp_output_dir / "test2.json").exists()

    def test_library_load_from_directory(self, tmp_output_dir):
        """Test loading library from directory with JSON files."""
        # Create JSON files
        solvent1 = Solvent(name="SOL1", description="Solvent 1")
        solvent2 = Solvent(name="SOL2", description="Solvent 2")

        solvent1.to_json(tmp_output_dir / "sol1.json")
        solvent2.to_json(tmp_output_dir / "sol2.json")

        # Load library
        library = SolventLibrary(library_path=tmp_output_dir)

        assert len(library) == 2
        assert library.get("SOL1") is not None
        assert library.get("SOL2") is not None


class TestStandardSolvents:
    """Test standard solvent definitions."""

    def test_create_standard_solvents(self):
        """Test creating standard solvents."""
        solvents = create_standard_solvents()

        assert len(solvents) >= 3  # WAT, ETA, MAM at minimum

        names = [s.name for s in solvents]
        assert "WAT" in names
        assert "ETA" in names

    def test_standard_solvent_probes(self):
        """Test standard solvents have probes defined."""
        solvents = create_standard_solvents()

        for solvent in solvents:
            assert len(solvent.probes) > 0, f"{solvent.name} has no probes"

    def test_ethanol_solvent(self):
        """Test ethanol solvent definition."""
        solvents = create_standard_solvents()
        eta = next(s for s in solvents if s.name == "ETA")

        # Should have ETA and WAT residues
        residue_names = eta.get_residue_names()
        assert "ETA" in residue_names
        assert "WAT" in residue_names

        # Should have OH and CT probes
        probe_names = eta.get_probe_names()
        assert "OH" in probe_names
        assert "CT" in probe_names

        # Should have volume and probe types
        assert eta.volume is not None
        assert eta.volume > 0

        oh_probe = eta.get_probe("OH")
        assert "Don" in oh_probe.probe_types
        assert "Acc" in oh_probe.probe_types


class TestSolventLibraryIntegration:
    """Integration tests for SolventLibrary with actual JSON files."""

    def test_load_all_solvents_with_extra_fields(self):
        """Test loading all solvents from the library with extra fields."""
        library = SolventLibrary()

        # Should load at least 9 solvents
        assert len(library) >= 9

        # Check each solvent has the extra fields
        for solvent in library:
            # Volume should be set
            assert solvent.volume is not None, f"{solvent.name} missing volume"
            assert solvent.volume > 0, f"{solvent.name} has invalid volume"

            # is_ionic should be a bool
            assert isinstance(solvent.is_ionic, bool), f"{solvent.name} is_ionic not bool"

            # total_charge should be set
            assert isinstance(solvent.total_charge, (int, float)), (
                f"{solvent.name} total_charge invalid"
            )

            # Each probe should have types
            for probe in solvent.probes:
                assert len(probe.probe_types) > 0, (
                    f"{solvent.name}/{probe.name} missing probe_types"
                )
                assert probe.probability is not None, (
                    f"{solvent.name}/{probe.name} missing probability"
                )

    def test_ionic_solvent_correctly_marked(self):
        """Test that ION solvent is marked as ionic."""
        library = SolventLibrary()
        ion = library.get("ION")

        assert ion is not None
        assert ion.is_ionic is True

        # Check non-ionic solvents
        eta = library.get("ETA")
        assert eta.is_ionic is False

    def test_probe_types_mapping(self):
        """Test that probe types are correctly loaded."""
        library = SolventLibrary()

        # ETA should have Don,Acc for OH and Hyd for CT
        eta = library.get("ETA")
        types_map = eta.get_types_map()

        assert "Don" in types_map.get("OH", "")
        assert "Acc" in types_map.get("OH", "")
        assert "Hyd" in types_map.get("CT", "")
        assert "Wat" in types_map.get("WAT", "")

    def test_probability_values_realistic(self):
        """Test that probability values are in realistic range."""
        library = SolventLibrary()

        for solvent in library:
            for probe in solvent.probes:
                # Probability should be small positive number
                assert probe.probability > 0, f"{solvent.name}/{probe.name} probability <= 0"
                assert probe.probability < 1, f"{solvent.name}/{probe.name} probability >= 1"
                # Typical range is 1e-6 to 1e-2
                assert probe.probability < 0.01, (
                    f"{solvent.name}/{probe.name} probability too large"
                )

    def test_eta_probability_matches_expected(self):
        """Test that ETA probabilities match expected values from old DB."""
        library = SolventLibrary()
        eta = library.get("ETA")

        oh_probe = eta.get_probe("OH")
        ct_probe = eta.get_probe("CT")
        wat_probe = eta.get_probe("WAT")

        # Values from the old SOLVENTS.db
        assert abs(oh_probe.probability - 0.00026600970383814976) < 1e-10
        assert abs(ct_probe.probability - 0.00026600970383814976) < 1e-10
        assert abs(wat_probe.probability - 0.003270354594245488) < 1e-10

    def test_roundtrip_json(self, tmp_output_dir):
        """Test that solvent can be saved and loaded with all fields."""
        library = SolventLibrary()
        eta = library.get("ETA")

        # Save to JSON
        json_path = tmp_output_dir / "eta_test.json"
        eta.to_json(json_path)

        # Load back
        loaded = Solvent.from_json(json_path)

        # Verify all fields preserved
        assert loaded.volume == eta.volume
        assert loaded.is_ionic == eta.is_ionic
        assert loaded.total_charge == eta.total_charge

        for orig_probe, loaded_probe in zip(eta.probes, loaded.probes):
            assert loaded_probe.probe_types == orig_probe.probe_types
            assert loaded_probe.probability == orig_probe.probability
