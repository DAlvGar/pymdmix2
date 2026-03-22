"""
Actions Manager
===============

Parallel execution of analysis actions on multiple replicas.

Replaces legacy ActionsManager with modern async/concurrent approach.

Examples
--------
>>> manager = ActionsManager(ncpus=4)
>>> manager.add_replicas([replica1, replica2])
>>> manager.add_actions([DensityAction, ResidenceAction])
>>> results = manager.run()
"""

from __future__ import annotations

import logging
from abc import ABC, abstractmethod
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Protocol, TypeVar, cast

log = logging.getLogger(__name__)


# Type for replica-like objects
class ReplicaProtocol(Protocol):
    """Protocol for replica objects."""

    name: str
    path: Path

    def get_trajectory(self, **kwargs) -> Any:
        """Get trajectory reader."""
        ...


# Type for action results
ActionResult = dict[str, Any]

# Type variable for actions
T = TypeVar("T", bound="Action")


class Action(ABC):
    """
    Base class for analysis actions.

    Subclasses must implement:
    - action_name: class attribute with action identifier
    - run(): execute the action on a replica
    - postprocess(): process results after run
    """

    action_name: str = "base_action"

    def __init__(self, replica: ReplicaProtocol, **kwargs):
        """
        Initialize action for a replica.

        Parameters
        ----------
        replica : ReplicaProtocol
            Replica to analyze
        **kwargs
            Additional action parameters
        """
        self.replica = replica
        self.params = kwargs
        self.log = logging.getLogger(f"{self.__class__.__name__}")

    @abstractmethod
    def run(self, **kwargs) -> ActionResult:
        """
        Execute the action.

        Returns
        -------
        ActionResult
            Dictionary of results
        """
        ...

    def postprocess(self, results: ActionResult, **kwargs) -> None:
        """
        Post-process results (optional).

        Parameters
        ----------
        results : ActionResult
            Results from run()
        """
        pass


@dataclass
class Job:
    """A job combining a replica and an action."""

    replica: ReplicaProtocol
    action_class: type[Action]
    params: dict = field(default_factory=dict)

    @property
    def name(self) -> str:
        return f"{self.replica.name}:{self.action_class.action_name}"


@dataclass
class JobResult:
    """Result of a job execution."""

    job: Job
    success: bool
    result: ActionResult | None = None
    error: str | None = None


def _run_job(job: Job) -> JobResult:
    """
    Execute a single job (used by process pool).

    Parameters
    ----------
    job : Job
        Job to execute

    Returns
    -------
    JobResult
        Execution result
    """
    try:
        log.info(f"Running {job.name}")
        action = job.action_class(job.replica, **job.params)
        result = action.run(**job.params)
        return JobResult(job=job, success=True, result=result)
    except Exception as e:
        log.error(f"Job {job.name} failed: {e}")
        return JobResult(job=job, success=False, error=str(e))


class ActionsManager:
    """
    Manage parallel execution of analysis actions on replicas.

    Parameters
    ----------
    ncpus : int
        Number of CPUs for parallel execution
    use_threads : bool
        Use threads instead of processes (for I/O-bound work)

    Examples
    --------
    >>> manager = ActionsManager(ncpus=4)
    >>> manager.add_replicas([replica1, replica2])
    >>> manager.add_actions([DensityAction])
    >>> results = manager.run()
    >>> print(results[replica1.name]['DensityAction'])
    """

    def __init__(self, ncpus: int = 1, use_threads: bool = False):
        self.ncpus = ncpus
        self.use_threads = use_threads
        self.replicas: list[ReplicaProtocol] = []
        self.action_classes: list[type[Action]] = []
        self.results: dict[str, dict[str, ActionResult]] = {}
        self.log = logging.getLogger("ActionsManager")

    def add_replicas(self, replicas: ReplicaProtocol | list[ReplicaProtocol]) -> None:
        """
        Add replicas to analyze.

        Parameters
        ----------
        replicas : ReplicaProtocol or list
            Replica(s) to add
        """
        if not isinstance(replicas, list):
            replicas = [replicas]

        self.replicas.extend(replicas)
        self.log.debug(f"Added {len(replicas)} replicas")

    def add_actions(
        self,
        actions: type[Action] | list[type[Action]] | str | list[str],
    ) -> None:
        """
        Add analysis actions to run.

        Parameters
        ----------
        actions : Action class, list, or string name(s)
            Actions to add
        """
        action_list: list[type[Action] | str]
        if isinstance(actions, list):
            action_list = cast(list[type[Action] | str], actions)
        else:
            action_list = [actions]

        for action in action_list:
            if isinstance(action, str):
                # Try to resolve by name (for CLI use)
                raise NotImplementedError("String action names not yet supported")
            elif isinstance(action, type) and issubclass(action, Action):
                self.action_classes.append(action)
            else:
                raise TypeError(f"Expected Action class, got {type(action)}")

        self.log.debug(f"Added {len(action_list)} actions")

    def prepare_jobs(self, **params) -> list[Job]:
        """
        Create jobs from replicas × actions.

        Parameters
        ----------
        **params
            Parameters passed to each action

        Returns
        -------
        list[Job]
            List of jobs to execute
        """
        jobs = []
        for replica in self.replicas:
            for action_class in self.action_classes:
                jobs.append(
                    Job(
                        replica=replica,
                        action_class=action_class,
                        params=params,
                    )
                )

        self.log.info(
            f"Prepared {len(jobs)} jobs ({len(self.replicas)} replicas × {len(self.action_classes)} actions)"
        )
        return jobs

    def run(
        self,
        postprocess: bool = True,
        **params,
    ) -> dict[str, dict[str, ActionResult]]:
        """
        Execute all jobs.

        Parameters
        ----------
        postprocess : bool
            Run postprocessing after each action
        **params
            Parameters passed to actions

        Returns
        -------
        dict[str, dict[str, ActionResult]]
            Nested dict: replica_name -> action_name -> results
        """
        jobs = self.prepare_jobs(**params)

        if not jobs:
            self.log.warning("No jobs to run")
            return {}

        # Choose executor
        if self.ncpus > 1:
            ExecutorClass = ThreadPoolExecutor if self.use_threads else ProcessPoolExecutor
            self.log.info(f"Running {len(jobs)} jobs with {self.ncpus} workers")

            with ExecutorClass(max_workers=self.ncpus) as executor:
                futures = {executor.submit(_run_job, job): job for job in jobs}

                for future in as_completed(futures):
                    job_result = future.result()
                    self._store_result(job_result)

                    if postprocess and job_result.success:
                        self._postprocess_result(job_result)
        else:
            # Serial execution
            self.log.info(f"Running {len(jobs)} jobs serially")
            for job in jobs:
                job_result = _run_job(job)
                self._store_result(job_result)

                if postprocess and job_result.success:
                    self._postprocess_result(job_result)

        return self.results

    def _store_result(self, job_result: JobResult) -> None:
        """Store a job result."""
        replica_name = job_result.job.replica.name
        action_name = job_result.job.action_class.action_name

        if replica_name not in self.results:
            self.results[replica_name] = {}

        if job_result.success:
            assert job_result.result is not None
            self.results[replica_name][action_name] = job_result.result
            self.log.info(f"Completed: {job_result.job.name}")
        else:
            self.results[replica_name][action_name] = {"error": job_result.error}
            self.log.error(f"Failed: {job_result.job.name}: {job_result.error}")

    def _postprocess_result(self, job_result: JobResult) -> None:
        """Run postprocessing for a job result."""
        try:
            if job_result.result is None:
                return
            action = job_result.job.action_class(
                job_result.job.replica,
                **job_result.job.params,
            )
            action.postprocess(job_result.result, **job_result.job.params)
        except Exception as e:
            self.log.error(f"Postprocess failed for {job_result.job.name}: {e}")

    def get_results(
        self,
        replica_name: str | None = None,
        action_name: str | None = None,
    ) -> dict | ActionResult | None:
        """
        Get results with optional filtering.

        Parameters
        ----------
        replica_name : str | None
            Filter by replica name
        action_name : str | None
            Filter by action name

        Returns
        -------
        dict or ActionResult or None
            Filtered results
        """
        if replica_name is None:
            return self.results

        replica_results = self.results.get(replica_name, {})

        if action_name is None:
            return replica_results

        return replica_results.get(action_name)

    def summary(self) -> str:
        """Return a human-readable summary of stored results."""
        if not self.results:
            return "No results available"

        lines = []
        for replica_name in sorted(self.results):
            lines.append(replica_name)
            for action_name, result in sorted(self.results[replica_name].items()):
                status = "failed" if "error" in result else "completed"
                lines.append(f"  - {action_name}: {status}")

        return "\n".join(lines)
