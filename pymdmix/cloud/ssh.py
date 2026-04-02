"""
SSH Utilities for Remote Instance Access
==========================================

Optional helpers for streaming logs and executing commands on running
EC2 instances via SSH.

Requires paramiko:
    pip install pymdmix[cloud]

Security note
-------------
SSH host key validation is enforced by default (RejectPolicy). The caller
must either:

1. Provide a ``known_hosts_path`` pointing to a file containing the host key, or
2. Ensure the system's ``~/.ssh/known_hosts`` already contains an entry for the
   instance (e.g., populated by running ``ssh-keyscan <ip> >> ~/.ssh/known_hosts``
   after the instance boots).

Passing ``known_hosts_path=None`` (the default) falls back to the system
``~/.ssh/known_hosts``.

Examples
--------
>>> from pymdmix.cloud.ssh import tail_remote_log, run_remote_command
>>> tail_remote_log("1.2.3.4", "/path/to/key.pem", "/home/ubuntu/pymdmix/ETA_1/bootstrap.log")
"""

from __future__ import annotations

import logging
import time
from collections.abc import Generator
from pathlib import Path

from pymdmix.cloud import _require_paramiko

log = logging.getLogger(__name__)

_DEFAULT_SSH_PORT = 22
_DEFAULT_CONNECT_TIMEOUT = 30
_TAIL_CHUNK_SIZE = 4096


def _make_ssh_client(
    hostname: str,
    port: int,
    username: str,
    key_filename: str,
    known_hosts_path: str | Path | None,
):
    """Create a paramiko SSHClient with host key validation."""
    _require_paramiko()
    import paramiko

    client = paramiko.SSHClient()
    client.set_missing_host_key_policy(paramiko.RejectPolicy())

    # Load host keys: explicit path takes priority, then system known_hosts
    if known_hosts_path is not None:
        client.load_host_keys(str(known_hosts_path))
    else:
        client.load_system_host_keys()

    client.connect(
        hostname=hostname,
        port=port,
        username=username,
        key_filename=key_filename,
        timeout=_DEFAULT_CONNECT_TIMEOUT,
    )
    return client


def tail_remote_log(
    instance_ip: str,
    key_path: str | Path,
    log_path: str,
    ssh_user: str = "ubuntu",
    poll_interval: float = 2.0,
    ssh_port: int = _DEFAULT_SSH_PORT,
    known_hosts_path: str | Path | None = None,
) -> Generator[str, None, None]:
    """
    Stream lines from a remote log file via SSH (like ``tail -f``).

    Yields new lines as they appear. Stops when the remote process closes
    the channel or the caller breaks iteration.

    Parameters
    ----------
    instance_ip : str
        Public IP address of the EC2 instance.
    key_path : str | Path
        Path to the local SSH private key (.pem file).
    log_path : str
        Absolute path to the log file on the remote instance.
    ssh_user : str
        SSH username (default: 'ubuntu').
    poll_interval : float
        Seconds between reads when no new data is available.
    ssh_port : int
        SSH port (default: 22).
    known_hosts_path : str | Path | None
        Path to a known_hosts file containing the instance's host key.
        If None, the system ``~/.ssh/known_hosts`` is used. The instance's
        host key must be present before connecting (use ``ssh-keyscan``).

    Yields
    ------
    str
        Lines from the remote log.

    Examples
    --------
    >>> for line in tail_remote_log("1.2.3.4", "key.pem", "/tmp/md.log"):
    ...     print(line, end="")
    """
    client = _make_ssh_client(
        hostname=instance_ip,
        port=ssh_port,
        username=ssh_user,
        key_filename=str(key_path),
        known_hosts_path=known_hosts_path,
    )
    try:
        channel = client.get_transport().open_session()  # type: ignore[union-attr]
        channel.exec_command(f"tail -f {log_path}")

        buffer = ""
        while True:
            if channel.recv_ready():
                data = channel.recv(_TAIL_CHUNK_SIZE).decode(errors="replace")
                buffer += data
                while "\n" in buffer:
                    line, buffer = buffer.split("\n", 1)
                    yield line + "\n"
            elif channel.exit_status_ready():
                # Drain remaining output
                if buffer:
                    yield buffer
                break
            else:
                time.sleep(poll_interval)
    finally:
        client.close()


def run_remote_command(
    instance_ip: str,
    key_path: str | Path,
    command: str,
    ssh_user: str = "ubuntu",
    ssh_port: int = _DEFAULT_SSH_PORT,
    timeout: int = 60,
    known_hosts_path: str | Path | None = None,
) -> tuple[int, str, str]:
    """
    Execute a command on a remote EC2 instance via SSH.

    Parameters
    ----------
    instance_ip : str
        Public IP address of the EC2 instance.
    key_path : str | Path
        Path to the local SSH private key (.pem file).
    command : str
        Shell command to execute.
    ssh_user : str
        SSH username (default: 'ubuntu').
    ssh_port : int
        SSH port (default: 22).
    timeout : int
        Seconds to wait for the command to complete.
    known_hosts_path : str | Path | None
        Path to a known_hosts file containing the instance's host key.
        If None, the system ``~/.ssh/known_hosts`` is used.

    Returns
    -------
    tuple[int, str, str]
        (exit_code, stdout, stderr)

    Examples
    --------
    >>> rc, out, err = run_remote_command("1.2.3.4", "key.pem", "nvidia-smi")
    >>> print(out)
    """
    client = _make_ssh_client(
        hostname=instance_ip,
        port=ssh_port,
        username=ssh_user,
        key_filename=str(key_path),
        known_hosts_path=known_hosts_path,
    )
    try:
        stdin, stdout_obj, stderr_obj = client.exec_command(command, timeout=timeout)
        exit_code = stdout_obj.channel.recv_exit_status()
        stdout_text = stdout_obj.read().decode(errors="replace")
        stderr_text = stderr_obj.read().decode(errors="replace")
        return exit_code, stdout_text, stderr_text
    finally:
        client.close()
