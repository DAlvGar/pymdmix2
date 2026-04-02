"""
Cloud Infrastructure Support
=============================

AWS EC2 GPU support for running MD replicas on cloud instances.

Requires optional dependencies:
    pip install pymdmix[cloud]
    # or: boto3>=1.26 paramiko>=3.0

Usage
-----
>>> from pymdmix.cloud import AWSConfig, EC2Manager, S3Transfer
>>> cfg = AWSConfig(region="us-east-1", s3_bucket="my-bucket", key_pair_name="my-key")
>>> ec2 = EC2Manager(cfg)
>>> transfer = S3Transfer(cfg)
"""

from __future__ import annotations

try:
    import boto3  # noqa: F401

    HAS_BOTO3 = True
except ImportError:
    HAS_BOTO3 = False

try:
    import paramiko  # noqa: F401

    HAS_PARAMIKO = True
except ImportError:
    HAS_PARAMIKO = False


def _require_boto3() -> None:
    """Raise ImportError with install hint if boto3 is not available."""
    if not HAS_BOTO3:
        raise ImportError(
            "boto3 is required for cloud features. "
            "Install with: pip install pymdmix[cloud]  or  pip install boto3"
        )


def _require_paramiko() -> None:
    """Raise ImportError with install hint if paramiko is not available."""
    if not HAS_PARAMIKO:
        raise ImportError(
            "paramiko is required for SSH features. "
            "Install with: pip install pymdmix[cloud]  or  pip install paramiko"
        )


from pymdmix.cloud.config import AWSConfig  # noqa: E402
from pymdmix.cloud.ec2 import EC2Manager  # noqa: E402
from pymdmix.cloud.monitor import JobRegistry, RunJob  # noqa: E402
from pymdmix.cloud.ssh import run_remote_command, tail_remote_log  # noqa: E402
from pymdmix.cloud.transfer import S3Transfer  # noqa: E402

__all__ = [
    "HAS_BOTO3",
    "HAS_PARAMIKO",
    "AWSConfig",
    "EC2Manager",
    "S3Transfer",
    "RunJob",
    "JobRegistry",
    "tail_remote_log",
    "run_remote_command",
]
