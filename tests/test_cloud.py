"""Tests for the pymdmix cloud module."""

from __future__ import annotations

from unittest.mock import MagicMock

import pytest

# =============================================================================
# AWSConfig tests
# =============================================================================


class TestAWSConfig:
    def test_defaults(self):
        from pymdmix.cloud.config import AWSConfig

        cfg = AWSConfig()
        assert cfg.region == "us-east-1"
        assert cfg.instance_type == "g4dn.xlarge"
        assert cfg.use_spot is True
        assert cfg.ssh_user == "ubuntu"
        assert cfg.terminate_on_completion is True
        assert cfg.s3_prefix == "pymdmix/"

    def test_prefix_trailing_slash_normalized(self):
        from pymdmix.cloud.config import AWSConfig

        cfg = AWSConfig(s3_prefix="myprefix")
        assert cfg.s3_prefix == "myprefix/"

    def test_s3_replica_prefix(self):
        from pymdmix.cloud.config import AWSConfig

        cfg = AWSConfig(s3_prefix="mdmix/", s3_bucket="bucket")
        assert cfg.s3_replica_prefix("ETA_1") == "mdmix/ETA_1/"

    def test_s3_uri(self):
        from pymdmix.cloud.config import AWSConfig

        cfg = AWSConfig(s3_bucket="my-bucket", s3_prefix="pfx/")
        assert cfg.s3_uri("REP_1") == "s3://my-bucket/pfx/REP_1/"

    def test_validate_missing_required(self):
        from pymdmix.cloud.config import AWSConfig

        cfg = AWSConfig(region="", s3_bucket="", key_pair_name="")
        errors = cfg.validate()
        assert any("region" in e for e in errors)
        assert any("s3_bucket" in e for e in errors)
        assert any("key_pair_name" in e for e in errors)

    def test_validate_ok(self):
        from pymdmix.cloud.config import AWSConfig

        cfg = AWSConfig(region="us-east-1", s3_bucket="bucket", key_pair_name="my-key")
        assert cfg.validate() == []

    def test_to_dict_from_dict_roundtrip(self):
        from pymdmix.cloud.config import AWSConfig

        cfg = AWSConfig(
            region="eu-west-1",
            instance_type="g5.xlarge",
            s3_bucket="test-bucket",
            key_pair_name="kp",
            use_spot=False,
        )
        d = cfg.to_dict()
        assert d["region"] == "eu-west-1"
        assert d["use_spot"] is False

        restored = AWSConfig.from_dict(d)
        assert restored.region == cfg.region
        assert restored.instance_type == cfg.instance_type
        assert restored.use_spot is False

    def test_ssh_key_path_from_env(self, monkeypatch, tmp_path):
        from pymdmix.cloud.config import AWSConfig

        key_file = tmp_path / "key.pem"
        key_file.touch()
        monkeypatch.setenv("PYMDMIX_AWS_KEY_PATH", str(key_file))

        cfg = AWSConfig()
        assert cfg.ssh_key_path == key_file

    def test_ssh_key_path_none_when_not_set(self, monkeypatch):
        from pymdmix.cloud.config import AWSConfig

        monkeypatch.delenv("PYMDMIX_AWS_KEY_PATH", raising=False)
        cfg = AWSConfig()
        assert cfg.ssh_key_path is None

    def test_unknown_fields_ignored_in_from_dict(self):
        from pymdmix.cloud.config import AWSConfig

        cfg = AWSConfig.from_dict({"region": "ap-northeast-1", "future_unknown_field": 42})
        assert cfg.region == "ap-northeast-1"


# =============================================================================
# AMI catalog tests
# =============================================================================


class TestAMICatalog:
    def test_load_ami_catalog(self):
        from pymdmix.cloud.config import load_ami_catalog

        catalog = load_ami_catalog()
        assert "regions" in catalog
        assert "instance_types" in catalog
        assert "ambertools_version" in catalog

    def test_catalog_has_expected_regions(self):
        from pymdmix.cloud.config import load_ami_catalog

        catalog = load_ami_catalog()
        regions = set(catalog["regions"].keys())
        assert "us-east-1" in regions
        assert "eu-west-1" in regions

    def test_catalog_instance_types(self):
        from pymdmix.cloud.config import load_ami_catalog

        catalog = load_ami_catalog()
        assert "g4dn.xlarge" in catalog["instance_types"]


# =============================================================================
# Config integration tests
# =============================================================================


class TestConfigAWSIntegration:
    def test_config_has_aws_config_field(self):
        from pymdmix.project.config import Config

        cfg = Config()
        assert cfg.aws_config is None

    def test_config_to_dict_includes_aws(self):
        from pymdmix.cloud.config import AWSConfig
        from pymdmix.project.config import Config

        aws = AWSConfig(region="us-west-2", s3_bucket="b", key_pair_name="kp")
        cfg = Config(aws_config=aws)
        d = cfg.to_dict()
        assert d["aws_config"] is not None
        assert d["aws_config"]["region"] == "us-west-2"

    def test_config_to_dict_no_aws(self):
        from pymdmix.project.config import Config

        cfg = Config()
        d = cfg.to_dict()
        assert d["aws_config"] is None

    def test_config_from_dict_with_aws(self):
        from pymdmix.project.config import Config

        data = {
            "aws_config": {
                "region": "us-east-2",
                "s3_bucket": "mybucket",
                "key_pair_name": "mykey",
            }
        }
        cfg = Config.from_dict(data)
        assert cfg.aws_config is not None
        assert cfg.aws_config.region == "us-east-2"

    def test_config_from_dict_without_aws(self):
        from pymdmix.project.config import Config

        cfg = Config.from_dict({})
        assert cfg.aws_config is None

    def test_config_validate_delegates_to_aws(self):
        from pymdmix.cloud.config import AWSConfig
        from pymdmix.project.config import Config

        bad_aws = AWSConfig(region="", s3_bucket="", key_pair_name="")
        cfg = Config(aws_config=bad_aws)
        errors = cfg.validate()
        # Should include AWS errors
        aws_errors = bad_aws.validate()
        for err in aws_errors:
            assert err in errors

    def test_config_roundtrip_json(self, tmp_path):
        from pymdmix.cloud.config import AWSConfig
        from pymdmix.project.config import Config

        aws = AWSConfig(region="us-west-1", s3_bucket="bkt", key_pair_name="kp", use_spot=False)
        cfg = Config(aws_config=aws)
        path = tmp_path / "config.json"
        cfg.to_json(path)
        loaded = Config.from_json(path)
        assert loaded.aws_config is not None
        assert loaded.aws_config.region == "us-west-1"
        assert loaded.aws_config.use_spot is False


# =============================================================================
# RunJob / JobRegistry tests
# =============================================================================


class TestRunJob:
    def test_defaults(self):
        from pymdmix.cloud.monitor import JOB_LAUNCHING, RunJob

        job = RunJob(replica_name="ETA_1", s3_prefix="pymdmix/ETA_1/")
        assert job.state == JOB_LAUNCHING
        assert job.is_active is True
        assert job.instance_id is None

    def test_to_dict_from_dict_roundtrip(self):
        from pymdmix.cloud.monitor import RunJob

        job = RunJob(
            replica_name="ETA_1",
            s3_prefix="pfx/ETA_1/",
            instance_id="i-abc123",
            public_ip="1.2.3.4",
        )
        d = job.to_dict()
        restored = RunJob.from_dict(d)
        assert restored.replica_name == "ETA_1"
        assert restored.instance_id == "i-abc123"
        assert restored.public_ip == "1.2.3.4"

    def test_unknown_fields_ignored(self):
        from pymdmix.cloud.monitor import RunJob

        job = RunJob.from_dict({"replica_name": "X", "s3_prefix": "p/", "new_future_field": 99})
        assert job.replica_name == "X"

    def test_elapsed_seconds_returns_float(self):
        from pymdmix.cloud.monitor import RunJob

        job = RunJob(replica_name="R", s3_prefix="p/")
        elapsed = job.elapsed_seconds
        assert elapsed is not None
        assert elapsed >= 0.0


class TestJobRegistry:
    def test_empty_registry(self, tmp_path):
        from pymdmix.cloud.monitor import JobRegistry

        reg = JobRegistry(tmp_path)
        assert reg.jobs == []
        assert reg.active_jobs == []

    def test_add_and_get(self, tmp_path):
        from pymdmix.cloud.monitor import JobRegistry, RunJob

        reg = JobRegistry(tmp_path)
        job = RunJob(replica_name="ETA_1", s3_prefix="pfx/ETA_1/", instance_id="i-abc")
        reg.add(job)

        found = reg.get("ETA_1")
        assert found is not None
        assert found.instance_id == "i-abc"

    def test_persist_and_reload(self, tmp_path):
        from pymdmix.cloud.monitor import JobRegistry, RunJob

        reg = JobRegistry(tmp_path)
        job = RunJob(replica_name="MAM_2", s3_prefix="pfx/MAM_2/", instance_id="i-xyz")
        reg.add(job)

        # Reload from disk
        reg2 = JobRegistry(tmp_path)
        assert len(reg2.jobs) == 1
        assert reg2.jobs[0].replica_name == "MAM_2"
        assert reg2.jobs[0].instance_id == "i-xyz"

    def test_update_job(self, tmp_path):
        from pymdmix.cloud.monitor import JOB_DONE, JOB_LAUNCHING, JobRegistry, RunJob

        reg = JobRegistry(tmp_path)
        job = RunJob(replica_name="ETA_1", s3_prefix="pfx/ETA_1/", state=JOB_LAUNCHING)
        reg.add(job)

        job.state = JOB_DONE
        reg.update(job)

        reloaded = JobRegistry(tmp_path)
        assert reloaded.get("ETA_1").state == JOB_DONE  # type: ignore[union-attr]

    def test_remove_job(self, tmp_path):
        from pymdmix.cloud.monitor import JobRegistry, RunJob

        reg = JobRegistry(tmp_path)
        reg.add(RunJob(replica_name="ETA_1", s3_prefix="pfx/ETA_1/"))
        reg.add(RunJob(replica_name="ETA_2", s3_prefix="pfx/ETA_2/"))

        removed = reg.remove("ETA_1")
        assert removed is True
        assert reg.get("ETA_1") is None
        assert reg.get("ETA_2") is not None

    def test_active_jobs_filter(self, tmp_path):
        from pymdmix.cloud.monitor import JOB_DONE, JOB_RUNNING, JobRegistry, RunJob

        reg = JobRegistry(tmp_path)
        reg.add(RunJob(replica_name="ETA_1", s3_prefix="pfx/ETA_1/", state=JOB_RUNNING))
        reg.add(RunJob(replica_name="ETA_2", s3_prefix="pfx/ETA_2/", state=JOB_DONE))
        active = reg.active_jobs
        assert len(active) == 1
        assert active[0].replica_name == "ETA_1"


# =============================================================================
# print_status_table smoke test
# =============================================================================


class TestPrintStatusTable:
    def test_prints_no_jobs(self, capsys):
        from pymdmix.cloud.monitor import print_status_table

        print_status_table([])
        captured = capsys.readouterr()
        assert "No cloud jobs" in captured.out

    def test_prints_jobs(self, capsys):
        from pymdmix.cloud.monitor import JOB_RUNNING, RunJob, print_status_table

        jobs = [
            RunJob(
                replica_name="ETA_1",
                s3_prefix="pfx/",
                state=JOB_RUNNING,
                instance_id="i-abc123",
                public_ip="1.2.3.4",
            )
        ]
        print_status_table(jobs)
        captured = capsys.readouterr()
        assert "ETA_1" in captured.out
        assert "running" in captured.out


# =============================================================================
# EC2Manager user-data builder (no AWS call needed)
# =============================================================================


class TestEC2ManagerUserData:
    def _make_replica(self, name: str) -> MagicMock:
        replica = MagicMock()
        replica.name = name
        return replica

    def test_user_data_contains_replica_name(self):
        from pymdmix.cloud.config import AWSConfig
        from pymdmix.cloud.ec2 import EC2Manager

        cfg = AWSConfig(region="us-east-1", s3_bucket="b", key_pair_name="k")
        replica = self._make_replica("ETA_1")
        script = EC2Manager.build_user_data(replica, cfg)
        assert "ETA_1" in script

    def test_user_data_contains_s3_bucket(self):
        from pymdmix.cloud.config import AWSConfig
        from pymdmix.cloud.ec2 import EC2Manager

        cfg = AWSConfig(region="us-east-1", s3_bucket="my-bucket", key_pair_name="k")
        replica = self._make_replica("MAM_1")
        script = EC2Manager.build_user_data(replica, cfg)
        assert "my-bucket" in script

    def test_user_data_no_termination_when_disabled(self):
        from pymdmix.cloud.config import AWSConfig
        from pymdmix.cloud.ec2 import EC2Manager

        cfg = AWSConfig(
            region="us-east-1", s3_bucket="b", key_pair_name="k", terminate_on_completion=False
        )
        replica = self._make_replica("ETA_1")
        script = EC2Manager.build_user_data(replica, cfg)
        assert "terminate-instances" not in script
        assert "skipping self-termination" in script

    def test_user_data_has_termination_when_enabled(self):
        from pymdmix.cloud.config import AWSConfig
        from pymdmix.cloud.ec2 import EC2Manager

        cfg = AWSConfig(
            region="us-east-1", s3_bucket="b", key_pair_name="k", terminate_on_completion=True
        )
        replica = self._make_replica("ETA_1")
        script = EC2Manager.build_user_data(replica, cfg)
        assert "terminate-instances" in script

    def test_user_data_shebang(self):
        from pymdmix.cloud.config import AWSConfig
        from pymdmix.cloud.ec2 import EC2Manager

        cfg = AWSConfig(region="us-east-1", s3_bucket="b", key_pair_name="k")
        script = EC2Manager.build_user_data(self._make_replica("R"), cfg)
        assert script.startswith("#!/bin/bash")


# =============================================================================
# CLI cloud commands (without real AWS)
# =============================================================================


class TestCloudCLI:
    @pytest.fixture
    def runner(self):
        from click.testing import CliRunner

        return CliRunner()

    def test_cloud_help(self, runner):
        from pymdmix.cli import cli

        result = runner.invoke(cli, ["cloud", "--help"])
        assert result.exit_code == 0
        assert "configure" in result.output
        assert "test-connection" in result.output
        assert "list-amis" in result.output
        assert "estimated-cost" in result.output

    def test_run_help(self, runner):
        from pymdmix.cli import cli

        result = runner.invoke(cli, ["run", "--help"])
        assert result.exit_code == 0
        assert "replica" in result.output
        assert "all" in result.output
        assert "status" in result.output
        assert "fetch" in result.output
        assert "cancel" in result.output
        assert "logs" in result.output

    def test_list_amis(self, runner):
        from pymdmix.cli import cli

        result = runner.invoke(cli, ["cloud", "list-amis"])
        assert result.exit_code == 0
        assert "us-east-1" in result.output
        assert "g4dn.xlarge" in result.output

    def test_test_connection_no_boto3(self, runner, monkeypatch):
        """Test-connection should fail gracefully when boto3 is absent."""
        import pymdmix.cloud as cloud_module

        monkeypatch.setattr(cloud_module, "HAS_BOTO3", False)

        from pymdmix.cli import cli

        result = runner.invoke(cli, ["cloud", "test-connection"])
        assert result.exit_code != 0
        assert "boto3" in result.output

    def test_cloud_configure_dry_run(self, runner, tmp_path):
        """Test configure command with all options provided (no prompts)."""
        from pymdmix.cli import cli
        from pymdmix.project import Project
        from pymdmix.project.config import Config

        # Create a minimal project
        proj = Project(name="test", config=Config(), path=tmp_path)
        proj.save()

        result = runner.invoke(
            cli,
            [
                "cloud",
                "configure",
                "-p",
                str(tmp_path),
                "--region",
                "eu-west-1",
                "--instance-type",
                "g4dn.xlarge",
                "--key-pair",
                "my-kp",
                "--s3-bucket",
                "my-bucket",
                "--s3-prefix",
                "pfx/",
                "--no-spot",
            ],
            catch_exceptions=False,
        )
        assert result.exit_code == 0
        assert "saved" in result.output.lower()

        # Verify it was persisted
        proj2 = Project.load(tmp_path)
        assert proj2.config.aws_config is not None
        assert proj2.config.aws_config.region == "eu-west-1"
        assert proj2.config.aws_config.s3_bucket == "my-bucket"

    def test_run_status_no_jobs(self, runner, tmp_path):
        """run status with no jobs shows helpful message."""
        from pymdmix.cli import cli
        from pymdmix.project import Project
        from pymdmix.project.config import Config

        proj = Project(name="test", config=Config(), path=tmp_path)
        proj.save()

        result = runner.invoke(cli, ["run", "status", "-p", str(tmp_path)])
        assert result.exit_code == 0
        assert "No cloud jobs" in result.output
