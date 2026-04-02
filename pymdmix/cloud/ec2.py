"""
EC2 Instance Lifecycle Management
===================================

Manages the full lifecycle of AWS EC2 GPU instances used to run MD replicas:
launch → wait for SSH → run → terminate.

Requires boto3:
    pip install pymdmix[cloud]

Examples
--------
>>> from pymdmix.cloud.config import AWSConfig
>>> from pymdmix.cloud.ec2 import EC2Manager
>>>
>>> cfg = AWSConfig(region="us-east-1", s3_bucket="my-bucket",
...                 key_pair_name="my-key", ami_id="ami-0abc123")
>>> ec2 = EC2Manager(cfg)
>>> instance_id, public_ip = ec2.launch_instance("ETA_1", user_data_script)
"""

from __future__ import annotations

import base64
import logging
import socket
import time
from typing import TYPE_CHECKING, Any

from pymdmix.cloud import _require_boto3

if TYPE_CHECKING:
    from pymdmix.cloud.config import AWSConfig
    from pymdmix.project.replica import Replica

log = logging.getLogger(__name__)

# Tag applied to all instances launched by pymdmix
MANAGED_TAG_KEY = "ManagedBy"
MANAGED_TAG_VALUE = "pymdmix"


class EC2Error(Exception):
    """Error during EC2 operations."""


class EC2Manager:
    """
    Manage AWS EC2 instance lifecycle for running MD replicas.

    Parameters
    ----------
    aws_config : AWSConfig
        Cloud configuration.
    """

    def __init__(self, aws_config: AWSConfig) -> None:
        _require_boto3()
        import boto3

        self.config = aws_config
        self._ec2 = boto3.client("ec2", region_name=aws_config.region)
        self._ec2_resource = boto3.resource("ec2", region_name=aws_config.region)

    # ------------------------------------------------------------------
    # Instance launch
    # ------------------------------------------------------------------

    def launch_instance(
        self,
        job_name: str,
        user_data_script: str,
        extra_tags: dict[str, str] | None = None,
    ) -> tuple[str, str | None]:
        """
        Launch an EC2 instance for a single MD job.

        Parameters
        ----------
        job_name : str
            Human-readable job name used as the Name tag and for logging.
        user_data_script : str
            Bash script embedded as EC2 user-data. Executed as root on boot.
        extra_tags : dict, optional
            Additional EC2 tags.

        Returns
        -------
        tuple[str, str | None]
            (instance_id, public_ip) — public_ip may be None for private subnets.
        """
        cfg = self.config
        tags = {
            MANAGED_TAG_KEY: MANAGED_TAG_VALUE,
            "Name": f"pymdmix-{job_name}",
            "PyMDMixJob": job_name,
        }
        tags.update(cfg.extra_tags)
        if extra_tags:
            tags.update(extra_tags)

        tag_specs = [
            {
                "ResourceType": "instance",
                "Tags": [{"Key": k, "Value": v} for k, v in tags.items()],
            }
        ]

        encoded_ud = base64.b64encode(user_data_script.encode()).decode()

        launch_kwargs: dict[str, Any] = {
            "ImageId": cfg.ami_id or self._get_default_ami(),
            "InstanceType": cfg.instance_type,
            "KeyName": cfg.key_pair_name,
            "MinCount": 1,
            "MaxCount": 1,
            "UserData": encoded_ud,
            "TagSpecifications": tag_specs,
            "BlockDeviceMappings": [
                {
                    "DeviceName": "/dev/sda1",
                    "Ebs": {
                        "VolumeSize": cfg.ebs_volume_gb,
                        "VolumeType": "gp3",
                        "DeleteOnTermination": True,
                    },
                }
            ],
        }

        if cfg.security_group_id:
            launch_kwargs["SecurityGroupIds"] = [cfg.security_group_id]
        if cfg.subnet_id:
            launch_kwargs["SubnetId"] = cfg.subnet_id
        if cfg.iam_instance_profile:
            launch_kwargs["IamInstanceProfile"] = {"Name": cfg.iam_instance_profile}

        if cfg.use_spot:
            launch_kwargs["InstanceMarketOptions"] = {
                "MarketType": "spot",
                "SpotOptions": {"SpotInstanceType": "one-time"},
            }

        log.info("Launching EC2 instance for job '%s' (type=%s)", job_name, cfg.instance_type)
        response = self._ec2.run_instances(**launch_kwargs)
        instance = response["Instances"][0]
        instance_id: str = instance["InstanceId"]
        public_ip: str | None = instance.get("PublicIpAddress")

        log.info("Launched instance %s for job '%s'", instance_id, job_name)
        return instance_id, public_ip

    # ------------------------------------------------------------------
    # Instance control
    # ------------------------------------------------------------------

    def terminate_instance(self, instance_id: str) -> None:
        """Terminate an EC2 instance permanently."""
        log.info("Terminating instance %s", instance_id)
        self._ec2.terminate_instances(InstanceIds=[instance_id])

    def stop_instance(self, instance_id: str) -> None:
        """Stop (hibernate) an EC2 instance without terminating."""
        log.info("Stopping instance %s", instance_id)
        self._ec2.stop_instances(InstanceIds=[instance_id])

    def get_instance_state(self, instance_id: str) -> str:
        """
        Return the current state of an EC2 instance.

        Returns
        -------
        str
            One of: pending, running, shutting-down, terminated, stopping, stopped
        """
        response = self._ec2.describe_instances(InstanceIds=[instance_id])
        reservations = response.get("Reservations", [])
        if not reservations:
            return "terminated"
        instance = reservations[0]["Instances"][0]
        return str(instance["State"]["Name"])

    def get_public_ip(self, instance_id: str) -> str | None:
        """Return the current public IP of a running instance."""
        response = self._ec2.describe_instances(InstanceIds=[instance_id])
        reservations = response.get("Reservations", [])
        if not reservations:
            return None
        instance = reservations[0]["Instances"][0]
        return instance.get("PublicIpAddress")

    def wait_for_ssh(
        self,
        instance_id: str,
        timeout_seconds: int = 300,
        poll_interval: int = 15,
    ) -> str:
        """
        Wait until the instance is running and SSH port 22 is open.

        Parameters
        ----------
        instance_id : str
            EC2 instance ID.
        timeout_seconds : int
            Maximum wait time before raising EC2Error.
        poll_interval : int
            Seconds between polls.

        Returns
        -------
        str
            Public IP address once SSH is available.

        Raises
        ------
        EC2Error
            If SSH is not reachable within timeout_seconds.
        """
        deadline = time.monotonic() + timeout_seconds
        log.info("Waiting for SSH on instance %s (timeout=%ds)", instance_id, timeout_seconds)

        while time.monotonic() < deadline:
            state = self.get_instance_state(instance_id)
            if state in ("terminated", "shutting-down"):
                raise EC2Error(f"Instance {instance_id} is {state} — cannot connect via SSH")

            if state == "running":
                ip = self.get_public_ip(instance_id)
                if ip and _is_port_open(ip, 22, timeout=5):
                    log.info("SSH available on %s (%s)", instance_id, ip)
                    return ip

            time.sleep(poll_interval)

        raise EC2Error(
            f"Timed out waiting for SSH on instance {instance_id} after {timeout_seconds}s"
        )

    # ------------------------------------------------------------------
    # Discovery
    # ------------------------------------------------------------------

    def list_running_instances(
        self, tag_filter: dict[str, str] | None = None
    ) -> list[dict[str, Any]]:
        """
        List EC2 instances managed by pymdmix.

        Parameters
        ----------
        tag_filter : dict, optional
            Additional tag filters (key → value). The ManagedBy=pymdmix
            filter is always applied.

        Returns
        -------
        list[dict]
            List of instance dicts with keys: instance_id, state, public_ip,
            job_name, launch_time.
        """
        filters = [
            {"Name": f"tag:{MANAGED_TAG_KEY}", "Values": [MANAGED_TAG_VALUE]},
            {"Name": "instance-state-name", "Values": ["pending", "running", "stopping"]},
        ]
        if tag_filter:
            for k, v in tag_filter.items():
                filters.append({"Name": f"tag:{k}", "Values": [v]})

        response = self._ec2.describe_instances(Filters=filters)
        results: list[dict[str, Any]] = []
        for reservation in response.get("Reservations", []):
            for inst in reservation["Instances"]:
                tags = {t["Key"]: t["Value"] for t in inst.get("Tags", [])}
                results.append(
                    {
                        "instance_id": inst["InstanceId"],
                        "state": inst["State"]["Name"],
                        "public_ip": inst.get("PublicIpAddress"),
                        "job_name": tags.get("PyMDMixJob", ""),
                        "launch_time": inst.get("LaunchTime"),
                    }
                )
        return results

    # ------------------------------------------------------------------
    # User-data script builder
    # ------------------------------------------------------------------

    @staticmethod
    def build_user_data(replica: Replica, aws_config: AWSConfig) -> str:
        """
        Build the EC2 user-data bootstrap script for a replica.

        The script:
        1. Configures the instance environment (AMBERHOME, PATH)
        2. Syncs staging data from S3
        3. Executes the pre-generated COMMANDS.sh
        4. Syncs results back to S3
        5. Writes a DONE / ERROR sentinel object
        6. Optionally self-terminates the instance

        Parameters
        ----------
        replica : Replica
            The replica to run. Must have topology, coordinates, and COMMANDS.sh
            already written (i.e., state >= READY).
        aws_config : AWSConfig
            Cloud configuration.

        Returns
        -------
        str
            Bash script suitable for EC2 user-data (cloud-init).
        """
        s3_uri = aws_config.s3_uri(replica.name)
        work_dir = f"/home/{aws_config.ssh_user}/pymdmix/{replica.name}"
        terminate_cmd = (
            'INSTANCE_ID=$(TOKEN=$(curl -s -X PUT -H "X-aws-ec2-metadata-token-ttl-seconds: 21600" '
            'http://169.254.169.254/latest/api/token) && '
            'curl -s -H "X-aws-ec2-metadata-token: $TOKEN" '
            'http://169.254.169.254/latest/meta-data/instance-id) && '
            'aws ec2 terminate-instances --region "$AWS_REGION" --instance-ids "$INSTANCE_ID"'
        )
        maybe_terminate = terminate_cmd if aws_config.terminate_on_completion else "echo 'skipping self-termination'"
        log_file = f"{work_dir}/bootstrap.log"

        script = f"""\
#!/bin/bash
# pymdmix bootstrap script — auto-generated, do not edit
set -euo pipefail
exec > >(tee -a {log_file}) 2>&1

echo "=== pymdmix bootstrap started $(date) ==="

# --- Environment ---
export AWS_REGION="{aws_config.region}"
export AMBERHOME="${{AMBERHOME:-/opt/amber}}"
export PATH="$AMBERHOME/bin:$PATH"
export LD_LIBRARY_PATH="${{LD_LIBRARY_PATH:-}}:$AMBERHOME/lib"

# --- Download staging data ---
mkdir -p {work_dir}
cd {work_dir}
echo "Syncing input from {s3_uri}"
aws s3 sync "{s3_uri}input/" "{work_dir}/" --only-show-errors

# --- Run MD ---
echo "Starting MD simulation..."
if bash COMMANDS.sh; then
    EXIT_CODE=0
    echo "MD simulation completed successfully"
else
    EXIT_CODE=$?
    echo "MD simulation FAILED with exit code $EXIT_CODE"
fi

# --- Upload results ---
echo "Syncing results to {s3_uri}output/"
aws s3 sync "{work_dir}/" "{s3_uri}output/" \\
    --exclude "*.in" --exclude "*.sh" --exclude "*.prmtop" --exclude "*.inpcrd" \\
    --only-show-errors

# --- Write sentinel ---
if [ "$EXIT_CODE" -eq 0 ]; then
    echo "$(date)" | aws s3 cp - "{s3_uri}DONE"
    echo "Sentinel DONE written"
else
    echo "EXIT_CODE=$EXIT_CODE  $(date)" | aws s3 cp - "{s3_uri}ERROR"
    echo "Sentinel ERROR written"
fi

echo "=== pymdmix bootstrap finished $(date) ==="

# --- Self-terminate ---
{maybe_terminate}
"""
        return script

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _get_default_ami(self) -> str:
        """Look up the default AMI for the configured region from the catalog."""
        from pymdmix.cloud.config import load_ami_catalog

        catalog = load_ami_catalog()
        region_info = catalog.get("regions", {}).get(self.config.region, {})
        ami_id = region_info.get("ami_id")
        if ami_id:
            return str(ami_id)

        # Fall back to the latest Ubuntu 22.04 LTS (non-GPU) AMI in the region
        # The user should provide a proper CUDA AMI via AWSConfig.ami_id
        log.warning(
            "No pre-baked pymdmix AMI found for region %s. "
            "Falling back to Ubuntu 22.04 LTS — AmberTools will be installed on first boot. "
            "Set AWSConfig.ami_id to a custom AMI to avoid this.",
            self.config.region,
        )
        return self._find_ubuntu_ami()

    def _find_ubuntu_ami(self) -> str:
        """Find the latest Ubuntu 22.04 LTS AMI in the configured region."""
        response = self._ec2.describe_images(
            Owners=["099720109477"],  # Canonical
            Filters=[
                {"Name": "name", "Values": ["ubuntu/images/hvm-ssd/ubuntu-jammy-22.04-amd64-server-*"]},
                {"Name": "state", "Values": ["available"]},
                {"Name": "architecture", "Values": ["x86_64"]},
            ],
        )
        images = sorted(
            response.get("Images", []),
            key=lambda img: img.get("CreationDate", ""),
            reverse=True,
        )
        if not images:
            raise EC2Error(f"Could not find Ubuntu 22.04 AMI in region {self.config.region}")
        return str(images[0]["ImageId"])


# ---------------------------------------------------------------------------
# Utility
# ---------------------------------------------------------------------------


def _is_port_open(host: str, port: int, timeout: float = 5.0) -> bool:
    """Return True if the TCP port is reachable."""
    try:
        with socket.create_connection((host, port), timeout=timeout):
            return True
    except (TimeoutError, OSError):
        return False
