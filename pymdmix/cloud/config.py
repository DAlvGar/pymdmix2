"""
AWS Cloud Configuration
========================

AWSConfig dataclass for configuring AWS EC2 GPU instances used to run
MD replicas remotely.

AWS credentials are NEVER stored here — they are resolved via the standard
boto3 credential chain:
  1. Environment variables: AWS_ACCESS_KEY_ID / AWS_SECRET_ACCESS_KEY
  2. ~/.aws/credentials
  3. IAM instance profile (if running on EC2)

The SSH private key path is read from the environment variable
PYMDMIX_AWS_KEY_PATH and is never stored in config files.

Examples
--------
>>> from pymdmix.cloud.config import AWSConfig
>>> cfg = AWSConfig(
...     region="us-east-1",
...     s3_bucket="my-pymdmix-bucket",
...     key_pair_name="my-keypair",
...     instance_type="g4dn.xlarge",
... )
>>> cfg.to_dict()
"""

from __future__ import annotations

import os
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any


@dataclass
class AWSConfig:
    """
    AWS configuration for cloud-based MD replica execution.

    Attributes
    ----------
    region : str
        AWS region (e.g., 'us-east-1')
    instance_type : str
        EC2 instance type. Defaults to g4dn.xlarge (1× NVIDIA T4).
    ami_id : str | None
        AMI ID with AmberTools + CUDA. If None, a bootstrap script installs
        AmberTools on first boot (adds ~10 min overhead).
    key_pair_name : str
        Name of the EC2 key pair for SSH access.
    security_group_id : str | None
        Security group ID. Must allow inbound SSH (port 22) from your IP.
    subnet_id : str | None
        Subnet ID for VPC placement.
    iam_instance_profile : str | None
        IAM instance profile name or ARN for S3 access without static keys.
    s3_bucket : str
        S3 bucket for staging replica data and results.
    s3_prefix : str
        Key prefix in the S3 bucket (default: 'pymdmix/').
    use_spot : bool
        Use Spot instances (~70% cost reduction). Default: True.
    ssh_user : str
        SSH username on the instance (default: 'ubuntu').
    terminate_on_completion : bool
        Terminate the instance when MD finishes. Default: True.
    poll_interval_seconds : int
        Polling interval for job status checks.
    """

    region: str = "us-east-1"
    instance_type: str = "g4dn.xlarge"
    ami_id: str | None = None
    key_pair_name: str = ""
    security_group_id: str | None = None
    subnet_id: str | None = None
    iam_instance_profile: str | None = None
    s3_bucket: str = ""
    s3_prefix: str = "pymdmix/"
    use_spot: bool = True
    ssh_user: str = "ubuntu"
    terminate_on_completion: bool = True
    poll_interval_seconds: int = 60

    # Additional EBS storage (GB) for large trajectory data
    ebs_volume_gb: int = 100

    # Tags applied to all EC2 resources
    extra_tags: dict[str, str] = field(default_factory=dict)

    def __post_init__(self) -> None:
        # Normalise trailing slash on prefix
        if self.s3_prefix and not self.s3_prefix.endswith("/"):
            self.s3_prefix = self.s3_prefix + "/"

    # ------------------------------------------------------------------
    # Derived helpers
    # ------------------------------------------------------------------

    @property
    def ssh_key_path(self) -> Path | None:
        """
        Path to the local SSH private key (.pem).

        Read from the PYMDMIX_AWS_KEY_PATH environment variable.
        Never stored in config files.
        """
        env_val = os.environ.get("PYMDMIX_AWS_KEY_PATH")
        if env_val:
            return Path(env_val)
        return None

    def s3_replica_prefix(self, replica_name: str) -> str:
        """S3 key prefix for a specific replica."""
        return f"{self.s3_prefix}{replica_name}/"

    def s3_uri(self, replica_name: str) -> str:
        """Full S3 URI for a replica's staging area."""
        return f"s3://{self.s3_bucket}/{self.s3_replica_prefix(replica_name)}"

    def validate(self) -> list[str]:
        """
        Validate the configuration.

        Returns
        -------
        list[str]
            List of validation errors (empty if valid).
        """
        errors: list[str] = []
        if not self.region:
            errors.append("aws_config.region must not be empty")
        if not self.s3_bucket:
            errors.append("aws_config.s3_bucket must not be empty")
        if not self.key_pair_name:
            errors.append("aws_config.key_pair_name must not be empty")
        if self.ebs_volume_gb < 30:
            errors.append("aws_config.ebs_volume_gb should be at least 30 GB")
        return errors

    # ------------------------------------------------------------------
    # Serialisation
    # ------------------------------------------------------------------

    def to_dict(self) -> dict[str, Any]:
        """Serialise to a plain dictionary (safe for JSON/YAML)."""
        d = asdict(self)
        # ssh_key_path is intentionally excluded (env-var only)
        return d

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> AWSConfig:
        """Deserialise from a dictionary."""
        known = {f.name for f in cls.__dataclass_fields__.values()}
        filtered = {k: v for k, v in data.items() if k in known}
        return cls(**filtered)


def load_ami_catalog() -> dict[str, Any]:
    """
    Load the bundled AMI catalog from pymdmix/data/cloud/amis.json.

    Returns
    -------
    dict
        AMI catalog with regions, instance types, and metadata.
    """
    import json
    from importlib.resources import files

    catalog_path = files("pymdmix.data.cloud").joinpath("amis.json")
    return json.loads(catalog_path.read_text())  # type: ignore[arg-type]
