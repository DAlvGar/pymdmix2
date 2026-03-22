# pyMDMix2 — Development Progress & Resume Guide

**Last Updated:** 2026-03-22
**Goal:** Fix cross-codebase inconsistencies (API signatures, engine classes, CLI, docs)

---

## Current State: ✅ Phases 1–5 Complete

### Latest Validation Snapshot
- `uv run mypy pymdmix` → **Success: no issues found in 42 source files**
- `uv run pytest tests/test_cli.py tests/test_system.py -q` → **50 passed**
- `uv run pytest tests/test_cli.py -q` → **17 passed**

### Historical Test Suite Baseline
```
644 passed, 8 skipped, 2 pre-existing failures
```
Run: `uv run pytest tests/ -q --tb=short`

### Known Pre-existing Failures (do not fix in this track)
| Test | Root Cause |
|------|-----------|
| `test_gromacs.py::test_create_index_groups` | Asserts `"Non-Protein"` but GROMACS outputs `"non-Protein"` |
| `test_project.py::TestConfig::test_config_validation` | Config validation returns no errors (not yet implemented) |

---

## What Was Done

### Tooling
| Item | Status |
|------|--------|
| `uv` as package manager (`uv sync --all-extras`) | ✅ |
| `.pre-commit-config.yaml` (ruff, ruff-format, pre-commit-hooks; mypy at manual stage) | ✅ |
| README updated with `uv` workflow | ✅ |

### Phase 1 — Action.run() Signature Consistency ✅

All `Action` subclasses now match `Action.run(self, trajectory=None, reference=None, output_dir=None, **kwargs)`.
Required parameters moved to keyword-only args with `None` defaults + early-return guards.

| File | Change |
|------|--------|
| `pymdmix/analysis/energy.py` | `EnergyResult.grid: Grid` → `grids: dict[str, Grid]`; `EnergyAction` → proper `Action` subclass decorated `@register_action("energy")` |
| `pymdmix/analysis/align.py` | `AlignmentResult.mean_rmsd` made proper typed field (`float = 0.0`) instead of runtime monkey-patch |
| `pymdmix/analysis/hotspots.py` | `HotspotAction.run()` reordered: trajectory/reference/output_dir first; `grids` optional kwarg |
| `pymdmix/analysis/residence.py` | `ResidenceAction.run()` reordered: trajectory/reference/output_dir first; `hotspot_coords` optional kwarg |
| `pymdmix/analysis/density.py` | `CpptrajDensityAction.run()` full signature rewrite; all required params optional kwargs with validation guards |

### Phase 2 — Engine Writer/Checker Classes ✅

| File | Added |
|------|-------|
| `pymdmix/engines/namd.py` | `NAMDWriterError`, `NAMDWriter(replica)`, `NAMDChecker(replica)` |
| `pymdmix/engines/openmm.py` | `OpenMMWriterError`, `OpenMMWriter(replica)`, `OpenMMChecker(replica)` |

**`NAMDWriter`**: writes `min.namd`, `eq1–2.namd`, per-step `md_{NNN}.namd` files + `COMMANDS.sh`.
**`NAMDChecker`**: checks `.out` files for `"End of program"`; `get_sim_volume()` parses `.xsc` box vectors.
**`OpenMMWriter`**: writes `min.py`, `eq1–5.py` (4× NVT heating + 1× NPT density), per-step `md_{NNN}.py` + `COMMANDS.sh`.
**`OpenMMChecker`**: checks `.log`/`.rst7`/`.chk` existence; `get_sim_volume()` uses parmed to read rst7 box.

`Replica.create_md_input()` already dispatches to `NAMDWriter` / `OpenMMWriter` when `settings.md_program` is `"NAMD"` / `"OPENMM"`.

---

## Remaining Work

### Phase 5 — Module Cleanup + mypy ✅ COMPLETE

Implemented and validated all planned Phase 5 targets:

| Area | Status | Notes |
|------|--------|-------|
| `engines/executor.py` | ✅ | typing/interface cleanup done |
| `analysis/manager.py` | ✅ | dual-Action architecture clarified and typed |
| `utils/tools.py` | ✅ | utility typing and safety fixes |
| `project/browser.py` | ✅ | path/navigation guards fixed |
| `io/off_manager.py` | ✅ | signature/typing cleanup done |
| Global mypy target | ✅ | now clean on full package (`uv run mypy pymdmix`) |

### Legacy API Compatibility Follow-up ✅ COMPLETE

To match `pyMDmix-port` project semantics and avoid CLI/type drift, `Project` compatibility fields and methods were restored:

- Reintroduced `systems` and `groups` project-level attributes.
- Added compatibility methods: `add_system()`, `create_replica()`, `create_group()`, `get_group()`, `remove_group()`.
- Included `systems`/`groups` in project serialization (`to_dict`/`from_dict`).

This addresses the post-refactor mismatch where legacy-style call sites expected attributes/methods that were removed during simplification.

### Phase 3 — Result Dataclasses & Missing Factory Methods ✅

| File | Added |
|------|-------|
| `engines/queue.py` | `QueueConfig.from_file(path)` — loads `.toml` or `.json` |
| `core/solvent.py` | `Solvent.from_file(path)` — dispatches by extension; `_from_cfg()` for INI format |
| `setup/prepare.py` | `StructurePreparationOptions` dataclass; `PrepareResult.modifications: list[str]`; `prepare_structure()` accepts `options=` and `Path` input |
| `setup/solvate.py` | `SolvationOptions` dataclass; `SolvateResult.save_coordinates()`; `SolvateResult.save_topology()`; `solvate_structure()` accepts `options=`; `output_dir` now optional |

---

### Phase 4 — CLI Fixes ✅ COMPLETE (644 tests green)

**File:** `pymdmix/cli.py`

| Location | Bug | Fix Applied |
|----------|-----|-------------|
| `setup_solvate` | `library.get()` in `except KeyError` — `.get()` returns `None`, never raises | Changed to `if solv is None:` |
| `setup_solvate` | `result.topology.write(str(top_path))` — `topology` is a `Path`, not parmed | Changed to `result.save_topology(top_path)` |
| `analyze_density` | `DensityAction(spacing=, nprocs=)` — no custom `__init__` → `TypeError` | `DensityAction()`, pass `spacing`/`n_workers` to `run()`, get trajectory/reference/probe_selections from replica |
| `analyze_density` | `action.run(replica)` — passes `Replica` as `TrajectoryReader` | `action.run(trajectory, reference=..., probe_selections=..., spacing=..., n_workers=...)` |
| `analyze_density` | `result.grids.items()` — `ActionResult` has no `.grids` | `result.metadata.get("probe_names", [])` |
| `analyze_energy` | `action.run(replica)` — `replica` must be keyword-only arg (after `*`) | `action.run(replica=replica)` |
| `analyze_energy` | `result.grids.items()` — `ActionResult` has no `.grids` | `result.metadata.get("energy_result").grids.items()` |
| `analyze_hotspots` | `HotspotAction(energy_cutoff=, min_size=)` — no custom `__init__` → `TypeError` | `HotspotAction()`, pass `grids`/`energy_cutoff`/`min_points` to `run()` |
| `analyze_hotspots` | `hotspots = action.run(replica)` treated as iterable list | Load grids from replica, `result = action.run(grids=..., ...)`, use `result.metadata["n_hotspots"]` |
| `analyze_residence` | `ResidenceAction(cutoff=, min_time=)` — no custom `__init__` → `TypeError` | `ResidenceAction()`, pass `trajectory` and `tolerance=cutoff` to `run()` |
| `analyze_residence` | `result.mean_residence` — `ActionResult` has no such attr | `result.metadata.get("mean_residence", 0.0)` |
| `queue_generate` | `generate_queue_script(replica, queue_config, template=template)` — wrong signature entirely | `generate_queue_script(config=queue_config, job_name=replica.name, commands=commands, work_dir=...)`, read commands from `COMMANDS.sh` |

---

### Optional Next Work (Post-Phase-5)

- Expand regression matrix beyond focused suites (`cli`, `system`) to full test suite refresh.
- Decide whether to keep lightweight compatibility (`systems` metadata dict) or reintroduce richer typed `System` objects in Project.

---

## How to Resume

```bash
cd ~/clawd/pymdmix2
uv run mypy pymdmix
uv run pytest tests/test_cli.py tests/test_system.py -q
```

Then implement phases in order. Each phase can be committed independently.

### Quick File Reference
| Phase | Primary File | Notes |
|-------|-------------|-------|
| 3.1 | `engines/queue.py` | Check `examples/replica.toml` for expected shape |
| 3.2 | `core/solvent.py` | JSON files live in `pymdmix/data/solvents/` |
| 3.3 | `setup/prepare.py` | Check `pyMDmix-port/pyMDMix/AutoPrepare.py` for context |
| 3.4 | `setup/solvate.py` | Check `pyMDmix-port/pyMDMix/Systems.py` for context |
| 4 | `cli.py` | `grep -n "TODO\|FIXME\|raise\|NotImplemented" pymdmix/cli.py` to find broken sites |
| 5.1 | `engines/executor.py` | — |
| 5.2 | `analysis/manager.py` | — |

---

## Architecture Quick Reference

### Action Plugin Pattern
```python
from pymdmix.analysis.base import Action, ActionResult, register_action

@register_action("my_action")
class MyAction(Action):
    def run(self, trajectory=None, reference=None, output_dir=None,
            *, my_param=None, **kwargs) -> ActionResult:
        if trajectory is None:
            return ActionResult(success=False, error="trajectory required")
        # ... do work ...
        return ActionResult(success=True, output_files=[...])
```

### ReplicaState Machine
```
CREATED → SETUP → READY → RUNNING → COMPLETE → ALIGNED → ANALYZED
                                                               ↓
                                                            ERROR (any stage)
```

### Engine Writer/Checker Dispatch (in `Replica`)
```python
# replica.create_md_input() dispatches on settings.md_program:
#   "AMBER"  → AmberWriter  / AmberChecker
#   "NAMD"   → NAMDWriter   / NAMDChecker
#   "OPENMM" → OpenMMWriter / OpenMMChecker
```

### Data Flow
```
PDB → AutoPrepare → OFF/prmtop
                          ↓
                    System.solvate()
                          ↓
                    Replica (state machine) ← MDSettings
                          ↓
               AmberEngine / NAMDEngine / OpenMMEngine
                          ↓
                    trajectory (.nc / .dcd / .xtc)
                          ↓
              CpptrajDensityAction  ← primary (shells out to cpptraj)
              DensityAction         ← fallback (pure Python)
                          ↓
                       Grid (DX)
                          ↓
                    HotspotAction
```

---

## Environment
```
Python:     3.10.12 (CPython)
uv:         0.10.0
parmed:     4.3.1
MDAnalysis: 2.9.0
numpy:      2.2.6
scipy:      1.15.3
click:      8.3.1
pytest:     9.0.2
mypy:       1.19.1
ruff:       0.15.7
```

```bash
uv sync --all-extras                              # install / sync
uv run pytest tests/ -q --tb=short               # run tests
uv run ruff check pymdmix && uv run ruff format pymdmix  # lint + format
uv run pre-commit run mypy --hook-stage manual   # type check
```
