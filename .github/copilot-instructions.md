# pyMDMix2 – Copilot Instructions

## Repositories
- `~/clawd/pymdmix2/` — **canonical repo**, all new development goes here
- `~/clawd/pyMDmix-port/` — **legacy reference only** (Python 2.7, Biskit). Read it for algorithmic intent; never commit to it.

## Architecture

Layered package under `pymdmix/`:

| Layer | Path | Purpose |
|-------|------|---------|
| Core types | `core/` | `Grid`, `Trajectory`/`Frame`, `parmed.Structure` wrappers, `Solvent`/`Probe` |
| Analysis | `analysis/` | Plugin actions: `DensityAction`, `ResidenceAction`, `HotspotAction`, `AlignAction` |
| Engines | `engines/` | MD engine interfaces: **`LeapSession`/`AmberEngine`** (primary), `GromacsEngine`, `OpenMMEngine`, `NAMDEngine` (secondary) |
| I/O | `io/` | File formats: DX/MRC grids, DCD parser, OFF manager, Amber parsers |
| Project | `project/` | `Project`, `Replica`, `ReplicaState`, `MDSettings`, `Config` |
| Setup | `setup/` | `AutoPrepare` (structure prep), `solvate` (LEaP solvation) |
| CLI | `cli.py` | Click-based `pymdmix` command with `create`/`add`/`analyze`/`info` subgroups |

**Data flow:** PDB → `AutoPrepare` → OFF/prmtop → `System.solvate()` → `Replica` (state machine) → `AmberEngine` → trajectory → `open_trajectory()` → **`CpptrajDensityAction`** (primary) / `DensityAction` (pure-Python fallback) → `Grid` → `HotspotAction`.

**Primary density workflow is `cpptraj_density`** (registered name), which shells out to `cpptraj` from AmberTools. It is faster for large systems and handles trajectory batching natively. Locate it with `AMBER_PTRAJ` env var or `cpptraj` on `PATH`. The pure-Python `DensityAction` (registered as `"density"`) is the fallback when cpptraj is unavailable, and supports the same parameters plus `n_workers` for multiprocessing.

## Key Conventions

### Dataclasses everywhere
All domain objects use `@dataclass` (not attrs, not pydantic). Use `field(default_factory=...)` for mutable defaults.

### parmed as the structure type
`parmed.Structure` is the core structure type throughout. Helper masks live in `core/structure.py` (e.g., `PROTEIN_RESIDUES`, `get_protein_mask()`, `find_disulfides()`). Never import Biskit.

### Analysis plugin registration
New analysis actions must use the `@register_action` decorator from `analysis/base.py`:
```python
from pymdmix.analysis.base import Action, ActionResult, register_action

@register_action("my_action")
class MyAction(Action):
    def run(self, trajectory, **kwargs) -> ActionResult:
        ...
```
There is a **second** `Action` base class in `analysis/manager.py` (for parallel dispatch). These are distinct—the registry pattern lives in `base.py`.

### Optional dependency guards
MDAnalysis, netCDF4, matplotlib are optional. Always guard:
```python
try:
    import MDAnalysis as mda
    HAS_MDANALYSIS = True
except ImportError:
    HAS_MDANALYSIS = False
```
Check `HAS_*` before use; raise `ImportError` with an install hint.

### `tomllib` compat
`tomllib` is stdlib in Python 3.11+. For 3.10 support use the existing pattern in `project/settings.py`:
```python
try:
    import tomllib
except ImportError:
    import tomli as tomllib
```

### TYPE_CHECKING imports
Use `if TYPE_CHECKING:` to avoid circular imports. All cross-layer type refs go there.

### Solvent library
Solvent definitions are JSON files in `pymdmix/data/solvents/` (e.g., `eta.json`), paired with AmberTools `.off` boxes (e.g., `ETAWAT20.off`). Load via `SolventLibrary` from `core/solvent.py`—never hardcode solvent data.

### ReplicaState machine
Replicas progress: `CREATED → SETUP → READY → RUNNING → COMPLETE → ALIGNED → ANALYZED → ERROR`. Only advance state after verifying files/output exist.

### HPC Queue Submission
Queue scripts are generated via `engines/queue.py`. `QueueConfig` supports `slurm`, `sge`, `pbs`, `local`. Generation is triggered from `Replica.create_queue_input(queue_system, ...)` or directly:
```python
from pymdmix.engines.queue import QueueConfig, generate_queue_script
config = QueueConfig(system="slurm", partition="gpu", n_gpus=1, time_hours=24)
script = generate_queue_script(config=config, job_name="md_production", commands=[...])
```
Queue templates are currently defined inline in `pymdmix/engines/queue.py` (`SLURM_TEMPLATE`, `SGE_TEMPLATE`, `PBS_TEMPLATE`, `LOCAL_TEMPLATE`). When adding a new queue system, update both the `QueueConfig` logic **and** the template mapping in `engines/queue.py`.

## Developer Workflows

**Package manager:** [`uv`](https://docs.astral.sh/uv/) — all dependency and environment management goes through `uv`.

```bash
# Install / sync all deps (from pymdmix2/)
uv sync --all-extras

# Run all tests
uv run pytest tests/

# Type-check
uv run mypy pymdmix

# Lint / auto-fix
uv run ruff check pymdmix
uv run ruff format pymdmix

# Run a specific test file
uv run pytest tests/test_grid.py -v

# Add a new runtime dependency
uv add <package>

# Add a dev-only dependency
uv add --dev <package>

# Update the lock file after editing pyproject.toml
uv lock
```

**AmberTools** must be on `PATH` (and `AMBERHOME` set) for engine/solvation tests; those are skipped automatically when absent.

## Critical Files
- [pymdmix/core/grid.py](../pymdmix/core/grid.py) — `Grid` dataclass with DX/XPLOR I/O, free-energy conversion, all grid math
- [pymdmix/core/structure.py](../pymdmix/core/structure.py) — residue sets (`PROTEIN_RESIDUES`, etc.) and mask helpers replacing Biskit `PDBModel`
- [pymdmix/analysis/base.py](../pymdmix/analysis/base.py) — `Action` ABC, `ActionResult`, `@register_action` registry
- [pymdmix/project/replica.py](../pymdmix/project/replica.py) — `Replica`, `ReplicaState`, folder layout constants
- [pymdmix/engines/amber.py](../pymdmix/engines/amber.py) — `LeapSession` subprocess control, `AmberEngine` input generation
- [tests/conftest.py](../tests/conftest.py) — shared fixtures (`tmp_output_dir`, `sample_coordinates`, `sample_dx_content`, etc.)
