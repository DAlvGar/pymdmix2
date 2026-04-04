"""
Tests for pymdmix.project module.
"""

from pymdmix.project.config import Config, MDSettings
from pymdmix.project.project import Project
from pymdmix.project.replica import Replica, ReplicaState, create_replica


class TestMDSettings:
    """Test MDSettings dataclass."""

    def test_default_settings(self):
        """Test default MD settings."""
        settings = MDSettings()

        assert settings.nsteps == 500000
        assert settings.timestep == 2.0
        assert settings.temperature == 300.0

    def test_custom_settings(self):
        """Test custom MD settings."""
        settings = MDSettings(nsteps=1000000, temperature=310.0)

        assert settings.nsteps == 1000000
        assert settings.temperature == 310.0

    def test_to_dict(self):
        """Test serialization to dict."""
        settings = MDSettings(nsteps=100000)
        data = settings.to_dict()

        assert data["nsteps"] == 100000
        assert "temperature" in data

    def test_from_dict(self):
        """Test deserialization from dict."""
        data = {"nsteps": 200000, "temperature": 298.0}
        settings = MDSettings.from_dict(data)

        assert settings.nsteps == 200000
        assert settings.temperature == 298.0

    def test_to_amber_mdin(self):
        """Test Amber mdin generation."""
        settings = MDSettings(nsteps=100000)
        mdin = settings.to_amber_mdin()

        assert "nstlim=100000" in mdin
        assert "&cntrl" in mdin


class TestConfig:
    """Test Config dataclass."""

    def test_default_config(self):
        """Test default configuration."""
        config = Config()

        assert config.default_solvent == "ETA"
        assert config.n_replicas == 3
        assert isinstance(config.md_settings, MDSettings)

    def test_config_with_amber_home(self, tmp_output_dir):
        """Test config with AMBER_HOME set."""
        config = Config(amber_home=tmp_output_dir)

        assert config.amber_home == tmp_output_dir

    def test_config_validation(self):
        """Test config validation."""
        config = Config(amber_home=None)
        errors = config.validate()

        assert len(errors) > 0
        assert any("AMBERHOME" in e for e in errors)

    def test_config_json_roundtrip(self, tmp_output_dir):
        """Test saving and loading config as JSON."""
        config = Config(
            default_solvent="MAM",
            n_replicas=5,
            md_settings=MDSettings(nsteps=100000),
        )

        json_path = tmp_output_dir / "config.json"
        config.to_json(json_path)

        loaded = Config.from_json(json_path)

        assert loaded.default_solvent == "MAM"
        assert loaded.n_replicas == 5
        assert loaded.md_settings.nsteps == 100000


class TestReplicaState:
    """Test ReplicaState enum."""

    def test_state_values(self):
        """Test state enum values exist."""
        assert ReplicaState.CREATED
        assert ReplicaState.SETUP
        assert ReplicaState.READY
        assert ReplicaState.COMPLETE
        assert ReplicaState.ANALYZED


class TestReplica:
    """Test Replica dataclass."""

    def test_replica_creation(self):
        """Test creating a replica."""
        replica = Replica(name="ETA_1", solvent="ETA")

        assert replica.name == "ETA_1"
        assert replica.solvent == "ETA"
        assert replica.state == ReplicaState.CREATED

    def test_replica_with_path(self, tmp_output_dir):
        """Test replica with path."""
        replica = Replica(
            name="ETA_1",
            solvent="ETA",
            path=tmp_output_dir / "ETA_1",
        )

        assert replica.path == tmp_output_dir / "ETA_1"
        assert not replica.exists()

    def test_create_directory(self, tmp_output_dir):
        """Test creating replica directory."""
        replica = Replica(
            name="ETA_1",
            solvent="ETA",
            path=tmp_output_dir / "ETA_1",
        )

        replica.create_directory()

        assert replica.exists()
        assert (replica.path / "min").exists()
        assert (replica.path / "eq").exists()

    def test_set_state(self, tmp_output_dir):
        """Test changing replica state."""
        replica = Replica(name="ETA_1", solvent="ETA")

        replica.set_state(ReplicaState.SETUP)
        assert replica.state == ReplicaState.SETUP

        replica.set_state(ReplicaState.READY)
        assert replica.state == ReplicaState.READY

    def test_replica_json_roundtrip(self, tmp_output_dir):
        """Test saving and loading replica."""
        replica = Replica(
            name="ETA_1",
            solvent="ETA",
            path=tmp_output_dir / "ETA_1",
            topology="system.prmtop",
            seed=12345,
        )
        replica.create_directory()
        replica.save()

        loaded = Replica.load(tmp_output_dir / "ETA_1" / "replica.json")

        assert loaded.name == "ETA_1"
        assert loaded.solvent == "ETA"
        assert loaded.topology == "system.prmtop"
        assert loaded.seed == 12345

    def test_create_replica_function(self, tmp_output_dir):
        """Test create_replica convenience function."""
        replica = create_replica(
            name="MAM_1",
            solvent="MAM",
            base_path=tmp_output_dir,
            seed=99999,
        )

        assert replica.name == "MAM_1"
        assert replica.exists()
        assert replica.seed == 99999


class TestProject:
    """Test Project class."""

    def test_project_creation(self, tmp_output_dir):
        """Test creating a project."""
        project = Project(
            name="test_project",
            path=tmp_output_dir / "test_project",
        )

        assert project.name == "test_project"
        assert project.n_replicas == 0

    def test_create_project(self, tmp_output_dir):
        """Test Project.create factory method."""
        project = Project.create(
            name="my_project",
            path=tmp_output_dir / "my_project",
        )

        assert project.path.exists()
        assert project.replicas_path.exists()
        assert project.input_path.exists()

    def test_add_replica(self, tmp_output_dir):
        """Test adding a replica."""
        project = Project.create(
            name="test",
            path=tmp_output_dir / "test",
        )

        replica = Replica(
            name="ETA_1",
            solvent="ETA",
            path=project.replicas_path / "ETA_1",
        )
        replica.create_directory()

        project.add_replica(replica)

        assert project.n_replicas == 1
        assert project.get_replica("ETA_1") is replica

    def test_add_replicas(self, tmp_output_dir):
        """Test adding multiple replicas."""
        project = Project.create(
            name="test",
            path=tmp_output_dir / "test",
        )

        replicas = project.add_replicas("ETA", n_replicas=3)

        assert len(replicas) == 3
        assert project.n_replicas == 3
        assert project.solvents == ["ETA"]

    def test_get_replicas_by_solvent(self, tmp_output_dir):
        """Test filtering replicas by solvent."""
        project = Project.create(
            name="test",
            path=tmp_output_dir / "test",
        )

        project.add_replicas("ETA", n_replicas=2)
        project.add_replicas("MAM", n_replicas=3)

        eta_replicas = project.get_replicas_by_solvent("ETA")
        mam_replicas = project.get_replicas_by_solvent("MAM")

        assert len(eta_replicas) == 2
        assert len(mam_replicas) == 3

    def test_project_status(self, tmp_output_dir):
        """Test project status."""
        project = Project.create(
            name="test",
            path=tmp_output_dir / "test",
        )
        project.add_replicas("ETA", n_replicas=2)

        status = project.status()

        assert status["n_replicas"] == 2
        assert "ETA" in status["solvents"]

    def test_project_json_roundtrip(self, tmp_output_dir):
        """Test saving and loading project."""
        project = Project.create(
            name="test",
            path=tmp_output_dir / "test",
        )
        project.add_replicas("ETA", n_replicas=2)
        project.save()

        loaded = Project.load(tmp_output_dir / "test")

        assert loaded.name == "test"
        assert loaded.n_replicas == 2
        assert len(loaded.solvents) == 1

    def test_project_iteration(self, tmp_output_dir):
        """Test iterating over project replicas."""
        project = Project.create(
            name="test",
            path=tmp_output_dir / "test",
        )
        project.add_replicas("ETA", n_replicas=3)

        names = [r.name for r in project]

        assert len(names) == 3
        assert "ETA_1" in names


# =============================================================================
# Workflow Integration Tests
# =============================================================================


class TestReplicaSettingsConversion:
    """Test Replica handles different MDSettings types in __post_init__."""

    def test_replica_mdsettings_from_settings_module(self, tmp_output_dir):
        """Replica converts settings.py MDSettings to replica.py MDSettings."""
        from pymdmix.project.replica import MDSettings as ReplicaMDSettings
        from pymdmix.project.settings import MDSettings as ProjMDSettings

        proj_settings = ProjMDSettings(
            solvent="ETA",
            nanos=10,
            temperature=310.0,
            production_steps=250000,
            trajectory_frequency=500,
            restraint_mode="BB",
            restraint_force=2.5,
        )

        replica = Replica(
            name="prot_ETA_1",
            solvent="ETA",
            path=tmp_output_dir / "prot_ETA_1",
            settings=proj_settings,
        )

        # After __post_init__, settings must be the local MDSettings
        assert isinstance(replica.settings, ReplicaMDSettings)
        assert replica.settings.nanos == 10
        assert replica.settings.temperature == 310.0
        assert replica.settings.prod_steps == 250000
        assert replica.settings.traj_frequency == 500
        assert replica.settings.restraint_mode == "BB"
        assert replica.settings.restraint_force == 2.5

    def test_replica_native_mdsettings_unchanged(self, tmp_output_dir):
        """Replica leaves native replica.MDSettings untouched."""
        from pymdmix.project.replica import MDSettings as ReplicaMDSettings

        native_settings = ReplicaMDSettings(nanos=30, temperature=298.0)
        replica = Replica(
            name="prot_ETA_1",
            solvent="ETA",
            path=tmp_output_dir / "prot_ETA_1",
            settings=native_settings,
        )

        assert replica.settings is native_settings


class TestProjectSystemsPath:
    """Test the new systems_path on Project."""

    def test_systems_path_property(self, tmp_output_dir):
        """Project.systems_path points to <project>/systems."""
        project = Project(
            name="test",
            path=tmp_output_dir / "test",
        )
        assert project.systems_path == tmp_output_dir / "test" / "systems"

    def test_create_directories_includes_systems(self, tmp_output_dir):
        """create_directories() creates the systems sub-folder."""
        project = Project(
            name="test",
            path=tmp_output_dir / "test",
        )
        project.create_directories()

        assert project.systems_path.exists()
        assert project.systems_path.is_dir()


class TestAmberWriterWorkflow:
    """Test AmberWriter produces correct files when topology/coordinates are set."""

    def test_write_replica_input_creates_files(self, tmp_output_dir):
        """AmberWriter creates min.in, eq1.in, eq2.in, prod.in."""
        from pymdmix.engines.amber import AmberWriter
        from pymdmix.project.replica import MDSettings as ReplicaMDSettings

        settings = ReplicaMDSettings(
            nanos=1, temperature=300.0, prod_steps=500000, traj_frequency=500
        )
        replica = Replica(
            name="prot_ETA_1",
            solvent="ETA",
            path=tmp_output_dir / "prot_ETA_1",
            settings=settings,
        )
        replica.create_directory()
        # Provide fake topology/coordinates so write_commands can reference them
        replica.topology = "system.prmtop"
        replica.coordinates = "system.inpcrd"

        writer = AmberWriter(replica)
        writer.write_replica_input()
        writer.write_commands()

        assert (replica.min_path / "min.in").exists()
        assert (replica.eq_path / "eq1.in").exists()
        assert (replica.eq_path / "eq2.in").exists()
        assert (replica.md_path / "prod.in").exists()
        assert (replica.path / "COMMANDS.sh").exists()

    def test_commands_sh_contains_expected_sections(self, tmp_output_dir):
        """COMMANDS.sh includes Minimization, Equilibration, Production blocks."""
        from pymdmix.engines.amber import AmberWriter
        from pymdmix.project.replica import MDSettings as ReplicaMDSettings

        settings = ReplicaMDSettings(
            nanos=1, temperature=300.0, prod_steps=500000, traj_frequency=500
        )
        replica = Replica(
            name="prot_ETA_1",
            solvent="ETA",
            path=tmp_output_dir / "prot_ETA_1",
            settings=settings,
        )
        replica.create_directory()
        replica.topology = "system.prmtop"
        replica.coordinates = "system.inpcrd"

        writer = AmberWriter(replica)
        writer.write_commands()

        commands = (replica.path / "COMMANDS.sh").read_text()
        assert "#!/bin/bash" in commands
        assert "Minimization" in commands
        assert "Equilibration" in commands
        assert "Production" in commands
        assert "system.prmtop" in commands
