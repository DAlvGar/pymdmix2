"""
Job Executor
============

Multi-threaded job execution system for running MD simulations
and analysis tasks in parallel.

Examples
--------
>>> from pymdmix.engines.executor import Executor
>>> executor = Executor(n_workers=4)
>>> executor.start()
>>> executor.submit_command("pmemd.cuda -i prod.in -o prod.out", path="/md/replica1")
>>> executor.submit_command("pmemd.cuda -i prod.in -o prod.out", path="/md/replica2")
>>> executor.wait()
>>> executor.terminate()
"""

from __future__ import annotations

import logging
import subprocess
import time
from collections.abc import Callable
from concurrent.futures import Future, ThreadPoolExecutor
from dataclasses import dataclass, field
from pathlib import Path
from queue import Empty, Queue
from threading import Event, Thread
from typing import Any

log = logging.getLogger(__name__)


@dataclass
class JobResult:
    """
    Result of a job execution.

    Attributes
    ----------
    success : bool
        Whether job completed successfully
    return_code : int | None
        Process return code (for commands)
    output : str | None
        Standard output
    error : str | None
        Standard error
    exception : Exception | None
        Exception if job raised one
    """

    success: bool
    return_code: int | None = None
    output: str | None = None
    error: str | None = None
    exception: Exception | None = None


@dataclass
class Job:
    """
    Job to be executed.

    Attributes
    ----------
    command : str | None
        Shell command to execute
    function : Callable | None
        Python function to call
    args : tuple
        Arguments for function
    kwargs : dict
        Keyword arguments for function
    cwd : Path | None
        Working directory for command
    """

    command: str | None = None
    function: Callable | None = None
    args: tuple = field(default_factory=tuple)
    kwargs: dict = field(default_factory=dict)
    cwd: Path | None = None

    @property
    def is_command(self) -> bool:
        """Whether this is a shell command job."""
        return self.command is not None

    def execute(self) -> JobResult:
        """Execute the job."""
        if self.is_command:
            return self._execute_command()
        else:
            return self._execute_function()

    def _execute_command(self) -> JobResult:
        """Execute shell command."""
        try:
            result = subprocess.run(
                self.command,
                shell=True,
                cwd=self.cwd,
                capture_output=True,
                text=True,
            )

            if result.returncode != 0:
                log.warning(f"Command exited with code {result.returncode}: {self.command}")

            return JobResult(
                success=(result.returncode == 0),
                return_code=result.returncode,
                output=result.stdout,
                error=result.stderr,
            )
        except Exception as e:
            log.error(f"Command failed: {self.command}: {e}")
            return JobResult(success=False, exception=e)

    def _execute_function(self) -> JobResult:
        """Execute Python function."""
        try:
            result = self.function(*self.args, **self.kwargs)
            return JobResult(success=True, output=str(result))
        except Exception as e:
            log.error(f"Function {self.function.__name__} failed: {e}")
            return JobResult(success=False, exception=e)


class Executor:
    """
    Multi-threaded job executor.

    Manages a pool of worker threads that execute shell commands
    or Python functions. Supports queue-based job submission
    and graceful shutdown.

    Parameters
    ----------
    n_workers : int
        Number of parallel workers (default: 1)
    wait_interval : float
        Seconds between queue polling (default: 0.5)

    Examples
    --------
    >>> executor = Executor(n_workers=4)
    >>> executor.start()
    >>>
    >>> # Submit shell commands
    >>> executor.submit_command("sleep 1", path="/tmp")
    >>> executor.submit_command("echo hello", path="/tmp")
    >>>
    >>> # Submit Python functions
    >>> executor.submit_function(print, "Hello", "World")
    >>>
    >>> # Wait for completion
    >>> executor.wait()
    >>> executor.terminate()

    >>> # Or use as context manager
    >>> with Executor(n_workers=2) as executor:
    ...     executor.submit_command("echo test")
    """

    def __init__(self, n_workers: int = 1, wait_interval: float = 0.5):
        self.n_workers = max(1, n_workers)
        self.wait_interval = wait_interval

        self._queue: Queue[Job | None] = Queue()
        self._workers: list[Thread] = []
        self._results: list[JobResult] = []
        self._running = Event()
        self._started = False

        log.debug(f"Executor initialized with {self.n_workers} workers")

    def start(self) -> None:
        """Start worker threads."""
        if self._started:
            log.warning("Executor already started")
            return

        self._running.set()
        self._started = True

        for i in range(self.n_workers):
            worker = Thread(
                target=self._worker_loop,
                name=f"Executor-Worker-{i}",
                daemon=True,
            )
            worker.start()
            self._workers.append(worker)

        log.debug(f"Started {self.n_workers} worker threads")

    def _worker_loop(self) -> None:
        """Worker thread main loop."""
        while self._running.is_set():
            try:
                job = self._queue.get(timeout=self.wait_interval)

                # None signals shutdown
                if job is None:
                    self._queue.task_done()
                    break

                result = job.execute()
                self._results.append(result)
                self._queue.task_done()

            except Empty:
                continue
            except Exception as e:
                log.error(f"Worker error: {e}")

    def submit_command(
        self,
        command: str,
        path: str | Path | None = None,
    ) -> None:
        """
        Submit a shell command for execution.

        Parameters
        ----------
        command : str
            Shell command to execute
        path : str or Path, optional
            Working directory for command
        """
        if not self._started:
            self.start()

        cwd = Path(path) if path else None
        job = Job(command=command, cwd=cwd)
        self._queue.put(job)

        log.debug(f"Submitted command: {command} (cwd={cwd or 'current'})")

    def submit_function(
        self,
        func: Callable,
        *args: Any,
        **kwargs: Any,
    ) -> None:
        """
        Submit a Python function for execution.

        Parameters
        ----------
        func : Callable
            Function to call
        *args
            Positional arguments for function
        **kwargs
            Keyword arguments for function
        """
        if not self._started:
            self.start()

        job = Job(function=func, args=args, kwargs=kwargs)
        self._queue.put(job)

        log.debug(f"Submitted function: {func.__name__}")

    def submit(
        self,
        command: str | None = None,
        func: Callable | None = None,
        path: str | Path | None = None,
        *args: Any,
        **kwargs: Any,
    ) -> None:
        """
        Submit a command or function (backwards compatible interface).

        Parameters
        ----------
        command : str, optional
            Shell command to execute
        func : Callable, optional
            Python function to call
        path : str or Path, optional
            Working directory for command
        *args, **kwargs
            Arguments for function
        """
        if command:
            self.submit_command(command, path=path)
        elif func:
            self.submit_function(func, *args, **kwargs)
        else:
            raise ValueError("Must provide either command or func")

    def wait(self, timeout: float | None = None) -> bool:
        """
        Wait for all jobs to complete.

        Parameters
        ----------
        timeout : float, optional
            Maximum seconds to wait (None = forever)

        Returns
        -------
        bool
            True if all jobs completed, False if timeout
        """
        start_time = time.time()

        while not self._queue.empty():
            if timeout and (time.time() - start_time) > timeout:
                return False
            time.sleep(self.wait_interval)

        # Wait for queue to be fully processed
        self._queue.join()
        return True

    def terminate(self) -> None:
        """Terminate all workers."""
        self._running.clear()

        # Send shutdown signals
        for _ in range(self.n_workers):
            self._queue.put(None)

        # Wait for workers to finish
        for worker in self._workers:
            worker.join(timeout=5.0)

        self._workers.clear()
        self._started = False

        log.debug("Executor terminated")

    def set_n_workers(self, n_workers: int) -> None:
        """
        Change number of workers (requires restart).

        Parameters
        ----------
        n_workers : int
            New number of workers
        """
        was_running = self._started
        if was_running:
            self.terminate()

        self.n_workers = max(1, n_workers)
        log.debug(f"Changed to {self.n_workers} workers")

        if was_running:
            self.start()

    @property
    def results(self) -> list[JobResult]:
        """Get all job results."""
        return list(self._results)

    @property
    def pending(self) -> int:
        """Number of pending jobs."""
        return self._queue.qsize()

    @property
    def is_running(self) -> bool:
        """Whether executor is running."""
        return self._started and self._running.is_set()

    def __enter__(self) -> Executor:
        """Context manager entry."""
        self.start()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        """Context manager exit."""
        self.wait()
        self.terminate()

    def __del__(self) -> None:
        """Cleanup on deletion."""
        if self._started:
            self.terminate()


class AsyncExecutor:
    """
    Async executor using ThreadPoolExecutor.

    Simpler alternative to Executor using concurrent.futures.
    Returns Future objects for tracking job completion.

    Parameters
    ----------
    n_workers : int
        Maximum number of parallel workers

    Examples
    --------
    >>> with AsyncExecutor(n_workers=4) as executor:
    ...     future1 = executor.submit_command("echo hello")
    ...     future2 = executor.submit_function(print, "world")
    ...     print(future1.result())  # Wait for result
    """

    def __init__(self, n_workers: int = 4):
        self.n_workers = n_workers
        self._pool: ThreadPoolExecutor | None = None

    def start(self) -> None:
        """Start the executor pool."""
        if self._pool is None:
            self._pool = ThreadPoolExecutor(max_workers=self.n_workers)

    def submit_command(
        self,
        command: str,
        path: str | Path | None = None,
    ) -> Future:
        """Submit a shell command."""
        if self._pool is None:
            self.start()

        job = Job(command=command, cwd=Path(path) if path else None)
        return self._pool.submit(job.execute)

    def submit_function(
        self,
        func: Callable,
        *args: Any,
        **kwargs: Any,
    ) -> Future:
        """Submit a function call."""
        if self._pool is None:
            self.start()

        return self._pool.submit(func, *args, **kwargs)

    def shutdown(self, wait: bool = True) -> None:
        """Shutdown the executor."""
        if self._pool is not None:
            self._pool.shutdown(wait=wait)
            self._pool = None

    def __enter__(self) -> AsyncExecutor:
        self.start()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        self.shutdown()


def run_command(
    command: str,
    cwd: str | Path | None = None,
    capture: bool = True,
    check: bool = True,
) -> subprocess.CompletedProcess:
    """
    Run a shell command.

    Convenience function for single command execution.

    Parameters
    ----------
    command : str
        Shell command to execute
    cwd : str or Path, optional
        Working directory
    capture : bool
        Whether to capture stdout/stderr
    check : bool
        Whether to raise on non-zero exit

    Returns
    -------
    subprocess.CompletedProcess
        Command result
    """
    return subprocess.run(
        command,
        shell=True,
        cwd=cwd,
        capture_output=capture,
        text=True,
        check=check,
    )
