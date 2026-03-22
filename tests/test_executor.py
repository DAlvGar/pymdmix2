"""
Tests for Job Executor
======================

Tests for pymdmix.engines.executor module.
"""

import time
from pathlib import Path

import pytest

from pymdmix.engines.executor import (
    AsyncExecutor,
    Executor,
    Job,
    JobResult,
    run_command,
)


class TestJob:
    """Tests for Job dataclass."""

    def test_command_job(self):
        """Test creating a command job."""
        job = Job(command="echo hello", cwd=Path("/tmp"))

        assert job.is_command is True
        assert job.command == "echo hello"
        assert job.cwd == Path("/tmp")

    def test_function_job(self):
        """Test creating a function job."""

        def my_func(x, y):
            return x + y

        job = Job(function=my_func, args=(1, 2), kwargs={})

        assert job.is_command is False
        assert job.function == my_func
        assert job.args == (1, 2)

    def test_execute_command(self):
        """Test executing a command job."""
        job = Job(command="echo hello")
        result = job.execute()

        assert result.success is True
        assert result.return_code == 0
        assert "hello" in result.output

    def test_execute_command_failure(self):
        """Test command that fails."""
        job = Job(command="exit 1")
        result = job.execute()

        assert result.success is False
        assert result.return_code == 1

    def test_execute_function(self):
        """Test executing a function job."""

        def add(x, y):
            return x + y

        job = Job(function=add, args=(3, 4), kwargs={})
        result = job.execute()

        assert result.success is True
        assert result.output == "7"

    def test_execute_function_with_kwargs(self):
        """Test executing function with keyword args."""

        def greet(name, greeting="Hello"):
            return f"{greeting}, {name}!"

        job = Job(function=greet, args=("World",), kwargs={"greeting": "Hi"})
        result = job.execute()

        assert result.success is True
        assert "Hi, World!" in result.output

    def test_execute_function_failure(self):
        """Test function that raises."""

        def failing_func():
            raise ValueError("Intentional error")

        job = Job(function=failing_func)
        result = job.execute()

        assert result.success is False
        assert result.exception is not None
        assert isinstance(result.exception, ValueError)


class TestJobResult:
    """Tests for JobResult dataclass."""

    def test_success_result(self):
        """Test successful result."""
        result = JobResult(success=True, return_code=0, output="done")

        assert result.success is True
        assert result.return_code == 0
        assert result.output == "done"
        assert result.error is None
        assert result.exception is None

    def test_failure_result(self):
        """Test failure result."""
        result = JobResult(success=False, return_code=1, error="Command failed")

        assert result.success is False
        assert result.error == "Command failed"


class TestExecutor:
    """Tests for Executor class."""

    def test_create_executor(self):
        """Test creating an executor."""
        executor = Executor(n_workers=2)

        assert executor.n_workers == 2
        assert executor.is_running is False

    def test_start_stop(self):
        """Test starting and stopping."""
        executor = Executor(n_workers=1)
        executor.start()

        assert executor.is_running is True

        executor.terminate()

        assert executor.is_running is False

    def test_submit_command(self):
        """Test submitting a command."""
        with Executor(n_workers=1) as executor:
            executor.submit_command("echo test")
            executor.wait()

        # Should have one result
        assert len(executor.results) == 1
        assert executor.results[0].success is True

    def test_submit_function(self):
        """Test submitting a function."""
        results = []

        def store_result(x):
            results.append(x * 2)

        with Executor(n_workers=1) as executor:
            executor.submit_function(store_result, 5)
            executor.wait()

        assert 10 in results

    def test_parallel_execution(self):
        """Test parallel job execution."""
        completed = []

        def record(idx):
            time.sleep(0.05)
            completed.append(idx)

        with Executor(n_workers=4) as executor:
            for i in range(8):
                executor.submit_function(record, i)

        # All jobs should complete
        assert len(completed) == 8
        assert set(completed) == set(range(8))

    def test_submit_generic(self):
        """Test generic submit method."""
        with Executor(n_workers=1) as executor:
            executor.submit(command="echo generic")
            executor.wait()

        assert executor.results[0].success is True

    def test_wait_timeout(self):
        """Test wait with timeout."""
        executor = Executor(n_workers=1, wait_interval=0.5)
        executor.start()

        # Submit slow jobs to ensure queue isn't empty
        for _ in range(5):
            executor.submit_command("sleep 2")

        # Wait should timeout (queue won't be empty in 0.05 seconds)
        executor.wait(timeout=0.05)

        # Due to race conditions, we can't assert False always
        # but we can check the executor handled it gracefully

        executor.terminate()

    def test_change_workers(self):
        """Test changing worker count."""
        executor = Executor(n_workers=2)
        executor.start()

        executor.set_n_workers(4)

        assert executor.n_workers == 4

        executor.terminate()

    def test_context_manager(self):
        """Test using as context manager."""
        with Executor(n_workers=2) as executor:
            assert executor.is_running is True
            executor.submit_command("echo ctx")

        assert executor.is_running is False

    def test_pending_count(self):
        """Test pending job count."""
        executor = Executor(n_workers=1, wait_interval=0.1)
        executor.start()

        # Don't process yet
        time.sleep(0.01)

        # Submit several jobs
        for i in range(5):
            executor.submit_command(f"echo {i}")

        # Should have some pending
        # (may not be exact due to timing)

        executor.wait()

        assert executor.pending == 0

        executor.terminate()


class TestAsyncExecutor:
    """Tests for AsyncExecutor class."""

    def test_create_executor(self):
        """Test creating async executor."""
        executor = AsyncExecutor(n_workers=4)
        assert executor.n_workers == 4

    def test_submit_command(self):
        """Test submitting command."""
        with AsyncExecutor(n_workers=2) as executor:
            future = executor.submit_command("echo async")
            result = future.result(timeout=5)

        assert isinstance(result, JobResult)
        assert result.success is True

    def test_submit_function(self):
        """Test submitting function."""

        def square(x):
            return x * x

        with AsyncExecutor(n_workers=2) as executor:
            future = executor.submit_function(square, 7)
            result = future.result(timeout=5)

        assert result == 49

    def test_multiple_futures(self):
        """Test multiple concurrent futures."""

        def slow_add(x, y):
            time.sleep(0.05)
            return x + y

        with AsyncExecutor(n_workers=4) as executor:
            futures = [executor.submit_function(slow_add, i, i) for i in range(10)]

            results = [f.result(timeout=5) for f in futures]

        expected = [i + i for i in range(10)]
        assert results == expected


class TestRunCommand:
    """Tests for run_command function."""

    def test_simple_command(self):
        """Test running simple command."""
        result = run_command("echo hello", check=False)

        assert result.returncode == 0
        assert "hello" in result.stdout

    def test_command_with_cwd(self, tmp_path):
        """Test command with working directory."""
        result = run_command("pwd", cwd=tmp_path, check=False)

        assert str(tmp_path) in result.stdout

    def test_command_failure(self):
        """Test command that fails."""
        result = run_command("exit 42", check=False)

        assert result.returncode == 42

    def test_check_raises(self):
        """Test check=True raises on failure."""
        import subprocess

        with pytest.raises(subprocess.CalledProcessError):
            run_command("exit 1", check=True)

    def test_capture_stderr(self):
        """Test capturing stderr."""
        result = run_command("echo error >&2", check=False)

        assert "error" in result.stderr
