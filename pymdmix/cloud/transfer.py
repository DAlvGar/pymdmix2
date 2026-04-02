"""
S3 Data Transfer for MD Replicas
==================================

Uploads replica input data to S3 staging, downloads results, and manages
completion sentinels.

Requires boto3:
    pip install pymdmix[cloud]

Examples
--------
>>> from pymdmix.cloud.config import AWSConfig
>>> from pymdmix.cloud.transfer import S3Transfer
>>>
>>> cfg = AWSConfig(region="us-east-1", s3_bucket="my-bucket", ...)
>>> s3 = S3Transfer(cfg)
>>> s3.upload_replica(replica)
>>> # ... wait for job ...
>>> status = s3.check_completion_marker(replica)
>>> if status == "done":
...     s3.download_results(replica)
...     s3.cleanup_staging(replica)
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING

from pymdmix.cloud import _require_boto3

if TYPE_CHECKING:
    from pymdmix.cloud.config import AWSConfig
    from pymdmix.project.replica import Replica

log = logging.getLogger(__name__)

# Files to upload as MD input (relative to replica.path)
_INPUT_GLOBS = [
    "*.prmtop",
    "*.inpcrd",
    "*.rst7",
    "COMMANDS.sh",
    "min/*.in",
    "eq/*.in",
    "md/*.in",
]

# Output directories to download back
_OUTPUT_FOLDERS = ["min", "eq", "md"]


class S3TransferError(Exception):
    """Error during S3 transfer operations."""


class S3Transfer:
    """
    Transfer replica data to/from S3 staging.

    Parameters
    ----------
    aws_config : AWSConfig
        Cloud configuration containing bucket and region.
    """

    def __init__(self, aws_config: AWSConfig) -> None:
        _require_boto3()
        import boto3

        self.config = aws_config
        self._s3 = boto3.client("s3", region_name=aws_config.region)

    # ------------------------------------------------------------------
    # Upload
    # ------------------------------------------------------------------

    def upload_replica(self, replica: Replica) -> list[str]:
        """
        Upload all input files for a replica to S3.

        Uploads: topology, coordinates, .in files, COMMANDS.sh.

        Parameters
        ----------
        replica : Replica
            The replica whose files should be uploaded.

        Returns
        -------
        list[str]
            S3 keys of uploaded objects.

        Raises
        ------
        S3TransferError
            If the replica path is not set or required files are missing.
        """
        if replica.path is None:
            raise S3TransferError(f"Replica '{replica.name}' has no path set")
        replica_path = replica.path
        if not replica_path.exists():
            raise S3TransferError(f"Replica path does not exist: {replica_path}")

        prefix = self.config.s3_replica_prefix(replica.name) + "input/"
        uploaded: list[str] = []

        # Collect files matching input globs
        files_to_upload: list[Path] = []
        for pattern in _INPUT_GLOBS:
            files_to_upload.extend(replica_path.glob(pattern))

        if not files_to_upload:
            raise S3TransferError(
                f"No input files found for replica '{replica.name}' in {replica_path}. "
                "Did you run 'pymdmix queue generate' first?"
            )

        for file_path in files_to_upload:
            rel = file_path.relative_to(replica_path)
            s3_key = prefix + str(rel)
            log.debug("Uploading %s → s3://%s/%s", file_path, self.config.s3_bucket, s3_key)
            self._s3.upload_file(str(file_path), self.config.s3_bucket, s3_key)
            uploaded.append(s3_key)

        log.info(
            "Uploaded %d files for replica '%s' to s3://%s/%s",
            len(uploaded),
            replica.name,
            self.config.s3_bucket,
            prefix,
        )
        return uploaded

    # ------------------------------------------------------------------
    # Download
    # ------------------------------------------------------------------

    def download_results(self, replica: Replica) -> list[Path]:
        """
        Download MD output from S3 back to the local replica directory.

        Downloads the contents of min/, eq/, and md/ from the S3 output prefix
        into the local replica path.

        Parameters
        ----------
        replica : Replica
            The replica to download results for.

        Returns
        -------
        list[Path]
            Local paths of downloaded files.
        """
        if replica.path is None:
            raise S3TransferError(f"Replica '{replica.name}' has no path set")
        replica_path = replica.path

        output_prefix = self.config.s3_replica_prefix(replica.name) + "output/"
        downloaded: list[Path] = []

        paginator = self._s3.get_paginator("list_objects_v2")
        for page in paginator.paginate(Bucket=self.config.s3_bucket, Prefix=output_prefix):
            for obj in page.get("Contents", []):
                key: str = obj["Key"]
                # Only download output folder contents
                rel_key = key[len(output_prefix):]
                if not any(rel_key.startswith(f"{folder}/") for folder in _OUTPUT_FOLDERS):
                    continue

                local_path = replica_path / rel_key
                local_path.parent.mkdir(parents=True, exist_ok=True)
                log.debug("Downloading s3://%s/%s → %s", self.config.s3_bucket, key, local_path)
                self._s3.download_file(self.config.s3_bucket, key, str(local_path))
                downloaded.append(local_path)

        log.info(
            "Downloaded %d files for replica '%s'",
            len(downloaded),
            replica.name,
        )
        return downloaded

    # ------------------------------------------------------------------
    # Completion markers
    # ------------------------------------------------------------------

    def check_completion_marker(self, replica: Replica) -> str:
        """
        Check whether the remote job wrote a completion sentinel to S3.

        The bootstrap script writes either:
          - ``<prefix>DONE``  — successful completion
          - ``<prefix>ERROR`` — job failed

        Returns
        -------
        str
            ``"done"``, ``"error"``, or ``"pending"`` (neither sentinel found).
        """
        prefix = self.config.s3_replica_prefix(replica.name)
        for sentinel, result in [("DONE", "done"), ("ERROR", "error")]:
            key = prefix + sentinel
            try:
                self._s3.head_object(Bucket=self.config.s3_bucket, Key=key)
                log.debug("Found sentinel %s for replica '%s'", sentinel, replica.name)
                return result
            except self._s3.exceptions.ClientError as exc:
                if exc.response["Error"]["Code"] in ("404", "NoSuchKey"):
                    continue
                raise
        return "pending"

    def get_bootstrap_log(self, replica: Replica) -> str | None:
        """
        Fetch the bootstrap.log from the running/completed instance via S3.

        Returns
        -------
        str | None
            Log content, or None if not yet available.
        """
        prefix = self.config.s3_replica_prefix(replica.name)
        key = prefix + "output/bootstrap.log"
        try:
            response = self._s3.get_object(Bucket=self.config.s3_bucket, Key=key)
            return response["Body"].read().decode(errors="replace")
        except Exception as exc:
            log.debug("Bootstrap log not yet available for replica '%s': %s", replica.name, exc)
            return None

    # ------------------------------------------------------------------
    # Cleanup
    # ------------------------------------------------------------------

    def cleanup_staging(self, replica: Replica) -> int:
        """
        Delete all S3 staging objects for a replica.

        Call this after successfully downloading results.

        Returns
        -------
        int
            Number of deleted objects.
        """
        prefix = self.config.s3_replica_prefix(replica.name)
        paginator = self._s3.get_paginator("list_objects_v2")
        to_delete: list[dict] = []

        for page in paginator.paginate(Bucket=self.config.s3_bucket, Prefix=prefix):
            for obj in page.get("Contents", []):
                to_delete.append({"Key": obj["Key"]})

        if not to_delete:
            return 0

        # S3 delete_objects supports up to 1000 keys per call
        batch_size = 1000
        deleted = 0
        for i in range(0, len(to_delete), batch_size):
            batch = to_delete[i : i + batch_size]
            self._s3.delete_objects(
                Bucket=self.config.s3_bucket,
                Delete={"Objects": batch, "Quiet": True},
            )
            deleted += len(batch)

        log.info("Deleted %d staging objects for replica '%s'", deleted, replica.name)
        return deleted
