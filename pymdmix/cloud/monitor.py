"""
Job Monitor
============

Tracks the state of cloud MD jobs and polls for completion.

Jobs are persisted to ``{project_path}/.cloud_jobs.json`` so status survives
across pymdmix invocations.

Examples
--------
>>> from pymdmix.cloud.monitor import JobRegistry
>>> registry = JobRegistry(project_path)
>>> registry.add(job)
>>> for job in registry.jobs:
...     print(job.replica_name, job.state)
"""

from __future__ import annotations

import json
import logging
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from pymdmix.cloud.config import AWSConfig
    from pymdmix.project.replica import Replica

log = logging.getLogger(__name__)

# Job states
JOB_LAUNCHING = "launching"
JOB_RUNNING = "running"
JOB_DONE = "done"
JOB_FAILED = "failed"
JOB_CANCELLED = "cancelled"

_REGISTRY_FILENAME = ".cloud_jobs.json"


@dataclass
class RunJob:
    """
    Represents a single cloud MD run job.

    Attributes
    ----------
    replica_name : str
        Name of the replica being run.
    instance_id : str | None
        EC2 instance ID, set after launch.
    s3_prefix : str
        S3 key prefix for this job's staging data.
    launched_at : str
        ISO-format timestamp of when the job was launched.
    state : str
        Current job state: launching, running, done, failed, cancelled.
    public_ip : str | None
        Public IP of the EC2 instance (may be None for private subnets).
    error_message : str | None
        Error description if state == failed.
    """

    replica_name: str
    s3_prefix: str
    state: str = JOB_LAUNCHING
    instance_id: str | None = None
    public_ip: str | None = None
    launched_at: str = field(default_factory=lambda: datetime.now(timezone.utc).isoformat())
    completed_at: str | None = None
    error_message: str | None = None

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> RunJob:
        known = {f.name for f in cls.__dataclass_fields__.values()}
        return cls(**{k: v for k, v in data.items() if k in known})

    @property
    def is_active(self) -> bool:
        """True if the job is still in progress."""
        return self.state in (JOB_LAUNCHING, JOB_RUNNING)

    @property
    def elapsed_seconds(self) -> float | None:
        """Wall-clock seconds since launch, or total duration if completed."""
        try:
            start = datetime.fromisoformat(self.launched_at)
            if self.completed_at:
                end = datetime.fromisoformat(self.completed_at)
            else:
                end = datetime.now(timezone.utc)
            return (end - start).total_seconds()
        except ValueError as exc:
            log.debug("Could not parse timestamps for elapsed time: %s", exc)
            return None


class JobRegistry:
    """
    Persist and manage RunJob records for a project.

    Parameters
    ----------
    project_path : Path
        Root directory of the pymdmix project.
    """

    def __init__(self, project_path: Path) -> None:
        self.project_path = project_path
        self._registry_path = project_path / _REGISTRY_FILENAME
        self._jobs: list[RunJob] = []
        self._load()

    @property
    def jobs(self) -> list[RunJob]:
        """All jobs (active and historical)."""
        return list(self._jobs)

    @property
    def active_jobs(self) -> list[RunJob]:
        """Jobs that are still launching or running."""
        return [j for j in self._jobs if j.is_active]

    def get(self, replica_name: str) -> RunJob | None:
        """Find a job by replica name (most recent if multiple)."""
        matches = [j for j in self._jobs if j.replica_name == replica_name]
        return matches[-1] if matches else None

    def add(self, job: RunJob) -> None:
        """Add a new job to the registry and persist."""
        self._jobs.append(job)
        self._save()

    def update(self, job: RunJob) -> None:
        """Update an existing job record and persist."""
        for i, existing in enumerate(self._jobs):
            if (
                existing.replica_name == job.replica_name
                and existing.launched_at == job.launched_at
            ):
                self._jobs[i] = job
                self._save()
                return
        # If not found, add it
        self.add(job)

    def remove(self, replica_name: str) -> bool:
        """Remove all jobs for a replica. Returns True if any were removed."""
        before = len(self._jobs)
        self._jobs = [j for j in self._jobs if j.replica_name != replica_name]
        if len(self._jobs) < before:
            self._save()
            return True
        return False

    # ------------------------------------------------------------------
    # Persistence
    # ------------------------------------------------------------------

    def _load(self) -> None:
        if not self._registry_path.exists():
            self._jobs = []
            return
        try:
            data = json.loads(self._registry_path.read_text())
            self._jobs = [RunJob.from_dict(j) for j in data.get("jobs", [])]
        except Exception as exc:
            log.warning("Could not load cloud job registry: %s", exc)
            self._jobs = []

    def _save(self) -> None:
        data = {"version": 1, "jobs": [j.to_dict() for j in self._jobs]}
        self._registry_path.write_text(json.dumps(data, indent=2, default=str))


# ---------------------------------------------------------------------------
# Polling
# ---------------------------------------------------------------------------


def poll_all_jobs(
    registry: JobRegistry,
    aws_config: AWSConfig,
    replicas: list[Replica] | None = None,
) -> list[RunJob]:
    """
    Poll EC2 and S3 to update the state of all active jobs.

    For each active job:
    1. Checks EC2 instance state (running / terminated / stopped)
    2. Checks S3 for DONE / ERROR sentinel
    3. If done: downloads results and marks job complete
    4. If instance terminated but no sentinel: marks job failed

    Parameters
    ----------
    registry : JobRegistry
        Persisted job registry.
    aws_config : AWSConfig
        Cloud configuration.
    replicas : list[Replica] | None
        Replica objects (needed for result download). If None, download is skipped.

    Returns
    -------
    list[RunJob]
        Updated job list.
    """
    from pymdmix.cloud.ec2 import EC2Manager
    from pymdmix.cloud.transfer import S3Transfer

    ec2 = EC2Manager(aws_config)
    s3 = S3Transfer(aws_config)

    replica_map: dict[str, Replica] = {}
    if replicas:
        replica_map = {r.name: r for r in replicas}

    for job in list(registry.active_jobs):
        if job.instance_id is None:
            continue

        try:
            instance_state = ec2.get_instance_state(job.instance_id)
            public_ip = ec2.get_public_ip(job.instance_id)
            if public_ip:
                job.public_ip = public_ip

            sentinel = s3.check_completion_marker(
                _FakeReplica(job.replica_name, aws_config)  # type: ignore[arg-type]
            )

            if sentinel == "done":
                job.state = JOB_DONE
                job.completed_at = datetime.now(timezone.utc).isoformat()
                replica = replica_map.get(job.replica_name)
                if replica:
                    try:
                        s3.download_results(replica)
                        log.info("Downloaded results for replica '%s'", job.replica_name)
                    except Exception as exc:
                        log.warning("Could not download results for '%s': %s", job.replica_name, exc)
            elif sentinel == "error":
                job.state = JOB_FAILED
                job.completed_at = datetime.now(timezone.utc).isoformat()
                job.error_message = "Remote MD job reported ERROR sentinel"
            elif instance_state in ("terminated", "shutting-down", "stopped"):
                job.state = JOB_FAILED
                job.completed_at = datetime.now(timezone.utc).isoformat()
                job.error_message = f"Instance {job.instance_id} is {instance_state} with no completion sentinel"
            else:
                # Still running
                job.state = JOB_RUNNING

            registry.update(job)

        except Exception as exc:
            log.warning("Error polling job for replica '%s': %s", job.replica_name, exc)

    return registry.jobs


# ---------------------------------------------------------------------------
# Display
# ---------------------------------------------------------------------------


def format_duration(seconds: float | None) -> str:
    """Format duration as HH:MM:SS."""
    if seconds is None:
        return "—"
    h = int(seconds // 3600)
    m = int((seconds % 3600) // 60)
    s = int(seconds % 60)
    return f"{h:02d}:{m:02d}:{s:02d}"


def print_status_table(jobs: list[RunJob]) -> None:
    """
    Print a formatted status table of cloud jobs to stdout.

    Uses click.echo for consistent output with the rest of the CLI.
    """
    import click

    if not jobs:
        click.echo("No cloud jobs found.")
        return

    # Header
    header = f"{'REPLICA':<30} {'STATE':<12} {'INSTANCE':<20} {'IP':<16} {'ELAPSED':<10}"
    click.echo(header)
    click.echo("-" * len(header))

    state_colors = {
        JOB_LAUNCHING: "yellow",
        JOB_RUNNING: "cyan",
        JOB_DONE: "green",
        JOB_FAILED: "red",
        JOB_CANCELLED: "magenta",
    }

    for job in sorted(jobs, key=lambda j: j.launched_at):
        color = state_colors.get(job.state, "white")
        state_str = click.style(f"{job.state:<12}", fg=color)
        instance_str = (job.instance_id or "—")[:20]
        ip_str = (job.public_ip or "—")[:16]
        elapsed_str = format_duration(job.elapsed_seconds)
        click.echo(f"{job.replica_name:<30} {state_str} {instance_str:<20} {ip_str:<16} {elapsed_str:<10}")


# ---------------------------------------------------------------------------
# Internal helper
# ---------------------------------------------------------------------------


class _FakeReplica:
    """Minimal duck-type Replica used for S3 key resolution in poll_all_jobs."""

    def __init__(self, name: str, aws_config: AWSConfig) -> None:
        self.name = name
        self.path = None
