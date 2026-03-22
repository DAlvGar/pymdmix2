"""Tests for ActionsManager."""

import pytest
from dataclasses import dataclass
from pathlib import Path

from pymdmix.analysis.manager import (
    ActionsManager,
    Action,
    ActionResult,
    Job,
    JobResult,
)


# Mock replica for testing
@dataclass
class MockReplica:
    name: str
    path: Path = Path("/tmp/test")
    
    def get_trajectory(self, **kwargs):
        return None


# Mock action for testing
class MockAction(Action):
    action_name = "mock_action"
    
    def run(self, **kwargs) -> ActionResult:
        return {"status": "completed", "replica": self.replica.name}


class FailingAction(Action):
    action_name = "failing_action"
    
    def run(self, **kwargs) -> ActionResult:
        raise RuntimeError("Intentional failure")


class PostprocessAction(Action):
    action_name = "postprocess_action"
    postprocess_called = False
    
    def run(self, **kwargs) -> ActionResult:
        return {"value": 42}
    
    def postprocess(self, results: ActionResult, **kwargs) -> None:
        PostprocessAction.postprocess_called = True


class TestJob:
    """Tests for Job dataclass."""

    def test_job_creation(self):
        replica = MockReplica(name="test_replica")
        job = Job(replica=replica, action_class=MockAction)
        assert job.replica == replica
        assert job.action_class == MockAction

    def test_job_name(self):
        replica = MockReplica(name="replica1")
        job = Job(replica=replica, action_class=MockAction)
        assert job.name == "replica1:mock_action"


class TestActionsManager:
    """Tests for ActionsManager."""

    def test_init_defaults(self):
        manager = ActionsManager()
        assert manager.ncpus == 1
        assert manager.use_threads is False
        assert manager.replicas == []
        assert manager.action_classes == []

    def test_init_with_cpus(self):
        manager = ActionsManager(ncpus=4)
        assert manager.ncpus == 4

    def test_add_single_replica(self):
        manager = ActionsManager()
        replica = MockReplica(name="test")
        manager.add_replicas(replica)
        assert len(manager.replicas) == 1
        assert manager.replicas[0].name == "test"

    def test_add_multiple_replicas(self):
        manager = ActionsManager()
        replicas = [MockReplica(name=f"r{i}") for i in range(3)]
        manager.add_replicas(replicas)
        assert len(manager.replicas) == 3

    def test_add_single_action(self):
        manager = ActionsManager()
        manager.add_actions(MockAction)
        assert len(manager.action_classes) == 1
        assert manager.action_classes[0] == MockAction

    def test_add_multiple_actions(self):
        manager = ActionsManager()
        manager.add_actions([MockAction, PostprocessAction])
        assert len(manager.action_classes) == 2

    def test_add_invalid_action(self):
        manager = ActionsManager()
        # String action names raise NotImplementedError (planned feature)
        with pytest.raises(NotImplementedError):
            manager.add_actions("not_an_action_class")

    def test_prepare_jobs(self):
        manager = ActionsManager()
        manager.add_replicas([MockReplica(name="r1"), MockReplica(name="r2")])
        manager.add_actions([MockAction, PostprocessAction])
        
        jobs = manager.prepare_jobs()
        assert len(jobs) == 4  # 2 replicas × 2 actions

    def test_run_serial(self):
        manager = ActionsManager(ncpus=1)
        manager.add_replicas(MockReplica(name="test_replica"))
        manager.add_actions(MockAction)
        
        results = manager.run()
        
        assert "test_replica" in results
        assert "mock_action" in results["test_replica"]
        assert results["test_replica"]["mock_action"]["status"] == "completed"

    def test_run_empty(self):
        manager = ActionsManager()
        results = manager.run()
        assert results == {}

    def test_run_with_failure(self):
        manager = ActionsManager(ncpus=1)
        manager.add_replicas(MockReplica(name="test"))
        manager.add_actions(FailingAction)
        
        results = manager.run()
        
        assert "test" in results
        assert "error" in results["test"]["failing_action"]

    def test_get_results_all(self):
        manager = ActionsManager(ncpus=1)
        manager.add_replicas(MockReplica(name="r1"))
        manager.add_actions(MockAction)
        manager.run()
        
        results = manager.get_results()
        assert "r1" in results

    def test_get_results_by_replica(self):
        manager = ActionsManager(ncpus=1)
        manager.add_replicas([MockReplica(name="r1"), MockReplica(name="r2")])
        manager.add_actions(MockAction)
        manager.run()
        
        r1_results = manager.get_results(replica_name="r1")
        assert "mock_action" in r1_results
        
    def test_get_results_by_action(self):
        manager = ActionsManager(ncpus=1)
        manager.add_replicas(MockReplica(name="r1"))
        manager.add_actions(MockAction)
        manager.run()
        
        action_result = manager.get_results(replica_name="r1", action_name="mock_action")
        assert action_result["status"] == "completed"

    def test_summary(self):
        manager = ActionsManager(ncpus=1)
        manager.add_replicas(MockReplica(name="test"))
        manager.add_actions(MockAction)
        manager.run()
        
        summary = manager.summary()
        assert "test" in summary
        assert "mock_action" in summary


class TestActionsManagerParallel:
    """Tests for parallel execution."""

    def test_run_parallel_threads(self):
        manager = ActionsManager(ncpus=2, use_threads=True)
        manager.add_replicas([MockReplica(name=f"r{i}") for i in range(4)])
        manager.add_actions(MockAction)
        
        results = manager.run()
        
        assert len(results) == 4
        for name in ["r0", "r1", "r2", "r3"]:
            assert name in results


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
