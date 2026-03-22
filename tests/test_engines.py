"""
Tests for pymdmix.engines module.
"""

import pytest

from pymdmix.engines.amber import AmberEngine, AmberInput
from pymdmix.engines.queue import QueueConfig, generate_queue_script


class TestAmberInput:
    """Test AmberInput dataclass."""

    def test_basic_input(self):
        """Test creating basic Amber input."""
        inp = AmberInput(
            title="Test",
            cntrl={"imin": 1, "maxcyc": 1000},
        )

        assert inp.title == "Test"
        assert inp.cntrl["imin"] == 1


class TestAmberEngine:
    """Test AmberEngine class."""

    def test_default_engine(self):
        """Test default engine settings."""
        engine = AmberEngine()
        assert engine.exe == "pmemd.cuda"

    def test_minimization_input(self):
        """Test minimization input generation."""
        engine = AmberEngine()
        mdin = engine.minimization_input(maxcyc=1000, ncyc=500)

        assert "imin=1" in mdin
        assert "maxcyc=1000" in mdin
        assert "ncyc=500" in mdin

    def test_minimization_with_restraints(self):
        """Test minimization with restraints."""
        engine = AmberEngine()
        mdin = engine.minimization_input(
            restraint_wt=5.0,
            restraint_mask="@CA,C,N",
        )

        assert "ntr=1" in mdin
        assert "restraint_wt=5.0" in mdin

    def test_heating_input(self):
        """Test heating input generation."""
        engine = AmberEngine()
        mdin = engine.heating_input(
            nsteps=25000,
            temp_start=0.0,
            temp_end=300.0,
        )

        assert "nstlim=25000" in mdin
        assert "tempi=0.0" in mdin
        assert "temp0=300.0" in mdin
        assert "ntp=0" in mdin  # NVT

    def test_equilibration_input(self):
        """Test equilibration input generation."""
        engine = AmberEngine()
        mdin = engine.equilibration_input(
            nsteps=50000,
            temp=300.0,
            pressure=1.0,
        )

        assert "nstlim=50000" in mdin
        assert "ntp=1" in mdin  # NPT
        assert "barostat=2" in mdin

    def test_production_input(self):
        """Test production input generation."""
        engine = AmberEngine()
        mdin = engine.production_input(
            nsteps=500000,
            temp=300.0,
            ntwx=500,
        )

        assert "nstlim=500000" in mdin
        assert "ntwx=500" in mdin
        assert "ntr=0" in mdin  # No restraints

    def test_run_command(self, tmp_output_dir):
        """Test run command generation."""
        engine = AmberEngine(exe="pmemd.cuda")

        cmd = engine.run_command(
            topology=tmp_output_dir / "system.prmtop",
            coordinates=tmp_output_dir / "system.inpcrd",
            input_file=tmp_output_dir / "prod.in",
            output_prefix="md1",
        )

        assert "pmemd.cuda" in cmd
        assert "-O" in cmd
        assert "-p" in cmd
        assert "-c" in cmd
        assert "md1.nc" in cmd


class TestQueueConfig:
    """Test QueueConfig dataclass."""

    def test_default_config(self):
        """Test default queue configuration."""
        config = QueueConfig()

        assert config.system == "slurm"
        assert config.n_gpus == 1
        assert config.time_hours == 24

    def test_custom_config(self):
        """Test custom queue configuration."""
        config = QueueConfig(
            system="pbs",
            partition="batch",
            n_cpus=8,
            memory_gb=32,
        )

        assert config.system == "pbs"
        assert config.n_cpus == 8


class TestGenerateQueueScript:
    """Test queue script generation."""

    def test_slurm_script(self):
        """Test SLURM script generation."""
        config = QueueConfig(
            system="slurm",
            partition="gpu",
            n_gpus=1,
            time_hours=12,
        )

        script = generate_queue_script(
            config=config,
            job_name="test_job",
            commands=["echo 'Hello'", "sleep 10"],
        )

        assert "#!/bin/bash" in script
        assert "#SBATCH --job-name=test_job" in script
        assert "#SBATCH --partition=gpu" in script
        assert "#SBATCH --gres=gpu:1" in script
        assert "echo 'Hello'" in script

    def test_sge_script(self):
        """Test SGE script generation."""
        config = QueueConfig(
            system="sge",
            partition="all.q",
            n_cpus=4,
        )

        script = generate_queue_script(
            config=config,
            job_name="sge_test",
            commands=["hostname"],
        )

        assert "#$$ -N sge_test" in script
        assert "#$$ -q all.q" in script

    def test_pbs_script(self):
        """Test PBS script generation."""
        config = QueueConfig(
            system="pbs",
            partition="batch",
            n_nodes=2,
        )

        script = generate_queue_script(
            config=config,
            job_name="pbs_test",
            commands=["mpirun -np 16 ./program"],
        )

        assert "#PBS -N pbs_test" in script
        assert "#PBS -q batch" in script

    def test_local_script(self):
        """Test local script generation."""
        config = QueueConfig(system="local")

        script = generate_queue_script(
            config=config,
            job_name="local_test",
            commands=["./run.sh"],
        )

        assert "#!/bin/bash" in script
        assert "./run.sh" in script

    def test_script_with_modules(self):
        """Test script with module loading."""
        config = QueueConfig(
            system="slurm",
            modules=["amber/22", "cuda/11.8"],
        )

        script = generate_queue_script(
            config=config,
            job_name="test",
            commands=["pmemd.cuda"],
        )

        assert "module load amber/22" in script
        assert "module load cuda/11.8" in script

    def test_script_with_environment(self):
        """Test script with environment variables."""
        config = QueueConfig(
            system="slurm",
            environment={"CUDA_VISIBLE_DEVICES": "0"},
        )

        script = generate_queue_script(
            config=config,
            job_name="test",
            commands=["nvidia-smi"],
        )

        assert "export CUDA_VISIBLE_DEVICES=0" in script

    def test_unknown_system(self):
        """Test unknown queue system raises error."""
        config = QueueConfig(system="unknown")

        with pytest.raises(ValueError, match="Unknown queue system"):
            generate_queue_script(config, "test", ["echo"])


class TestQueueConfigFromFile:
    """Test QueueConfig.from_file() factory method."""

    def test_from_json(self, tmp_path):
        """Test loading QueueConfig from a JSON file."""
        import json

        data = {
            "system": "slurm",
            "partition": "gpu",
            "n_gpus": 2,
            "time_hours": 48,
        }
        cfg_file = tmp_path / "queue.json"
        cfg_file.write_text(json.dumps(data))

        config = QueueConfig.from_file(cfg_file)

        assert config.system == "slurm"
        assert config.partition == "gpu"
        assert config.n_gpus == 2
        assert config.time_hours == 48

    def test_from_toml(self, tmp_path):
        """Test loading QueueConfig from a TOML file."""
        toml_content = """
system = "pbs"
partition = "batch"
n_cpus = 8
memory_gb = 64
"""
        cfg_file = tmp_path / "queue.toml"
        cfg_file.write_text(toml_content)

        config = QueueConfig.from_file(cfg_file)

        assert config.system == "pbs"
        assert config.partition == "batch"
        assert config.n_cpus == 8
        assert config.memory_gb == 64

    def test_from_file_unknown_fields_ignored(self, tmp_path):
        """Unknown fields in file are silently ignored."""
        import json

        data = {"system": "slurm", "unknown_future_field": "ignored"}
        cfg_file = tmp_path / "queue.json"
        cfg_file.write_text(json.dumps(data))

        config = QueueConfig.from_file(cfg_file)  # should not raise

        assert config.system == "slurm"

    def test_from_file_missing_raises(self):
        """Missing file raises FileNotFoundError."""
        with pytest.raises(FileNotFoundError):
            QueueConfig.from_file("/nonexistent/queue.toml")

    def test_from_file_bad_extension_raises(self, tmp_path):
        """Unsupported extension raises ValueError."""
        cfg_file = tmp_path / "queue.cfg"
        cfg_file.write_text("system=slurm")

        with pytest.raises(ValueError, match="Unsupported"):
            QueueConfig.from_file(cfg_file)
