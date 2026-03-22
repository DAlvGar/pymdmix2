# pyMDMix 2.0 - Migration Plan

**Goal**: Port pyMDMix from Python 2.7 to Python 3.10+, removing Biskit dependency.

**Strategy**: Extract needed Biskit functionality, rely on parmed/MDAnalysis for PDB handling.

---

## Repository Structure

```
~/clawd/pymdmix2/          ← CANONICAL REPO (all new development)
~/clawd/pyMDmix-port/      ← REFERENCE ONLY (legacy Python 2.7 code)
    └── pyMDMix/           ← Original source to read from
```

---

## Biskit Dependency Removal

### What Biskit provided:
| Biskit Module | Usage | Replacement |
|--------------|-------|-------------|
| `Biskit.PDBModel` | PDB handling, atom selection | `parmed.Structure` |
| `Biskit.tools` | Utility functions | `utils/tools.py` (extracted) |
| `Biskit.SettingsParser` | Config parsing | `utils/settings_parser.py` (extracted) |
| `Biskit.PDBParseFile` | PDB parsing | `parmed.load_file()` |
| `Biskit.mathUtils` | Math utilities | `numpy` |
| `Biskit.LogFile` | Logging | `logging` (stdlib) |
| `Biskit.test` | Testing | `pytest` |

### Extraction strategy:
1. Copy only needed functions from Biskit modules
2. Modernize (type hints, pathlib, dataclasses)
3. Replace PDBModel with parmed throughout

---

## Module Migration Status

### ✅ COMPLETE (32 modules)

| Legacy | New Location | Notes |
|--------|--------------|-------|
| `Align.py` | `analysis/align.py` | - |
| `Amber.py` | `engines/amber.py` | LeapSession, AmberChecker, AmberWriter added |
| `AutoPrepare.py` | `setup/prepare.py` | PDB2PQR interface |
| `Browser.py` | `project/browser.py` | Filesystem navigation |
| `containers.py` | `core/containers.py` | Atom, Residue, Probe |
| `Energy.py` | `analysis/energy.py` | Density→ΔG conversion |
| `Executor.py` | `engines/executor.py` | Parallel job execution |
| `GridData.py` | `core/grid.py` | DX/XPLOR I/O, radial ops |
| `GridsManager.py` | `core/grid.py` | Merged into Grid class |
| `NamdDCDParser.py` | `io/dcd_parser.py` | DCD trajectory reader |
| `NAMD.py` | `engines/namd.py` | NAMD input generation |
| `OFFManager.py` | `io/off_manager.py` | Amber OFF file handling |
| `OpenMM.py` | `engines/openmm.py` | OpenMM interface |
| `OpenMMUtils.py` | `engines/openmm.py` | Merged |
| `Parsers.py` | `io/parsers.py` | Config file parsing |
| `Plotter.py` | `io/plotting.py` | Visualization |
| `Projects.py` | `project/project.py` | Project container |
| `QueueWriting.py` | `engines/queue.py` | SLURM/PBS/SGE scripts |
| `Replicas.py` | `project/replica.py` | Replica management |
| `SettingsParser.py` | `utils/settings_parser.py` | Config management |
| `settings.py` | (simplified) | No Biskit, env detection |
| `Solvents.py` | `core/solvent.py` | JSON-based solvent library |
| `Structures.py` | `core/structure.py` | parmed-based structure ops |
| `Systems.py` | `core/system.py` | System class |
| `tools.py` | `utils/tools.py` | Utility functions |
| `Trajectory.py` | `core/trajectory.py` | MDAnalysis-based |
| - | `cli.py` | New Click-based CLI |
| - | `engines/gromacs.py` | New GROMACS support |
| - | `setup/solvate.py` | Solvation workflow |
| - | `analysis/density.py` | Density calculation |
| - | `analysis/residence.py` | Residence time analysis |
| - | `analysis/hotspots.py` | Basic hotspot detection |

### ❌ NOT YET PORTED (4 modules)

| Legacy | Status | Priority |
|--------|--------|----------|
| `MDSettings.py` | Not started | P1 |
| `PDB.py` | Not started | P2 |
| `Analysis.py` | Partial (ActionsManager missing) | P1 |
| `HotSpotsManager.py` | Partial (clustering missing) | P1 |

---

## Implementation Tasks

### Phase 1: Complete Core (P0) ✅ DONE
- [x] Grid with DX/XPLOR I/O
- [x] Replica with folder creation, status checks
- [x] Amber with LeapSession, AmberChecker
- [x] CLI with all commands wired

### Phase 2: Missing Modules (P1)

#### 2.1 MDSettings Class
**Source**: `pyMDmix-port/pyMDMix/MDSettings.py`
**Target**: `pymdmix2/pymdmix/project/settings.py`

```python
@dataclass
class MDSettings:
    """MD simulation parameters."""
    solvent: str
    nanos: int = 20
    temperature: float = 300.0
    timestep: float = 2.0  # fs
    restraint_mode: str = "FREE"  # FREE, BB, HA, CUSTOM
    restraint_force: float = 0.0
    restraint_mask: str = ""
    align_mask: str = ""
    # ... parsing from config files
```

#### 2.2 ActionsManager
**Source**: `pyMDmix-port/pyMDMix/Analysis.py`
**Target**: `pymdmix2/pymdmix/analysis/manager.py`

```python
class ActionsManager:
    """Run analysis actions on multiple replicas in parallel."""
    def __init__(self, ncpus: int = 1):
        ...
    def add_replicas(self, replicas: list[Replica]):
        ...
    def add_actions(self, actions: list[Action]):
        ...
    def run(self):
        # ThreadPoolExecutor-based parallel execution
```

#### 2.3 HotSpotSet & Clustering
**Source**: `pyMDmix-port/pyMDMix/HotSpotsManager.py`
**Target**: `pymdmix2/pymdmix/analysis/hotspots.py` (expand)

```python
@dataclass
class HotSpotSet:
    """Collection of hotspots with clustering."""
    hotspots: list[Hotspot]
    
    def cluster(self, cutoff: float = 2.5) -> list[HotSpotCluster]:
        # scipy.cluster.hierarchy
        ...
    
    def prune_by_energy(self, max_energy: float):
        ...
    
    def prune_by_volume(self, min_volume: float):
        ...
```

### Phase 3: Secondary Modules (P2)

#### 3.1 SolvatedPDB
**Source**: `pyMDmix-port/pyMDMix/PDB.py`  
**Decision**: May not be needed - parmed handles this

The old SolvatedPDB was a Biskit.PDBModel subclass. Review if `parmed.Structure` + our `Replica.get_pdb()` covers the use cases.

---

## Testing

Current: **549 tests passing**

Each new module needs:
- Unit tests in `tests/test_<module>.py`
- Integration test if it interacts with files

Run: `cd ~/clawd/pymdmix2 && ~/.local/bin/pytest tests/ -q`

---

## File Locations

| Purpose | Location |
|---------|----------|
| New code | `~/clawd/pymdmix2/pymdmix/` |
| Tests | `~/clawd/pymdmix2/tests/` |
| This plan | `~/clawd/pymdmix2/docs/PLAN.md` |
| Legacy reference | `~/clawd/pyMDmix-port/pyMDMix/` |

---

## Quick Reference

```bash
# Run tests
cd ~/clawd/pymdmix2 && ~/.local/bin/pytest tests/ -q

# Check legacy code
grep -n "def something" ~/clawd/pyMDmix-port/pyMDMix/SomeFile.py

# Line counts
wc -l ~/clawd/pymdmix2/pymdmix/**/*.py | tail -5
```

---

## Progress Tracking

Update this section as work progresses:

| Date | Work Done |
|------|-----------|
| 2026-03-20 | Initial scaffolding, core modules |
| 2026-03-21 | Solvent library, CLI, tests |
| 2026-03-22 | Expanded replica, amber, grid; repo consolidation |
| Next | MDSettings, ActionsManager, HotSpotSet clustering |
