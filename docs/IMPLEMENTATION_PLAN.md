# pyMDMix Implementation Plan

**Project:** Complete rewrite of pyMDMix for Python 3.10+  
**Start Date:** 2026-03-20  
**Target:** Production-ready v0.3.0

---

## Overview

### Goals
1. Replace Biskit dependency with parmed + MDAnalysis
2. Modern Python 3.10+ codebase with type hints
3. Support external trajectories (Gromacs, NAMD, etc.)
4. Plugin architecture for analysis actions
5. Clean, maintainable, documented code

### Non-Goals (for v0.3.0)
- GUI interface
- Cloud/HPC integration beyond queue scripts
- Real-time visualization

---

## Phase 1: Core Infrastructure ✅ COMPLETE

**Status:** Done  
**Completed:** 2026-03-20

### Deliverables
- [x] Package structure (`pymdmix/`)
- [x] `core/grid.py` - 3D density grid with DX I/O
- [x] `core/trajectory.py` - MDAnalysis + Amber NetCDF readers
- [x] `core/structure.py` - parmed utilities and atom masks
- [x] `core/solvent.py` - JSON-based solvent definitions
- [x] `pyproject.toml` - modern packaging

### Files Created
```
pymdmix/
├── __init__.py
├── core/
│   ├── __init__.py
│   ├── grid.py          (490 lines)
│   ├── trajectory.py    (380 lines)
│   ├── structure.py     (370 lines)
│   └── solvent.py       (340 lines)
└── pyproject.toml
```

---

## Phase 2: Analysis Actions 🔲 IN PROGRESS

**Status:** Not started  
**Estimated:** 1-2 days

### 2.1 Base Action Framework
**File:** `pymdmix/analysis/base.py`

- [ ] `Action` abstract base class
- [ ] `ActionResult` dataclass
- [ ] `@register_action` decorator for plugin discovery
- [ ] `get_action()`, `list_actions()` registry functions

```python
# Target API
@register_action("density")
class DensityAction(Action):
    def run(self, replica, trajectory, **kwargs) -> ActionResult:
        ...
```

### 2.2 Density Grid Calculation
**File:** `pymdmix/analysis/density.py`

Port from: `pyMDMix/Actions/Density.py` (855 lines)

- [ ] `DensityAction` class
- [ ] Per-probe grid calculation
- [ ] Multiprocessing support (optional)
- [ ] Progress reporting
- [ ] Output to DX format

**Key logic to preserve:**
```python
# From old code - core algorithm
for frame in trajectory:
    for probe in probes:
        coords = frame[probe_mask]
        for coord in coords:
            grid.add_count(coord)
density = grid / n_frames
```

### 2.3 Residence Time Analysis
**File:** `pymdmix/analysis/residence.py`

Port from: `pyMDMix/Actions/Residence.py` (310 lines)

- [ ] `ResidenceAction` class
- [ ] Hotspot proximity tracking
- [ ] Per-residue residence times
- [ ] Output format (CSV/JSON)

### 2.4 Hotspot Detection
**File:** `pymdmix/analysis/hotspots.py`

Port from: `pyMDMix/HotSpotsManager.py` (1050 lines)

- [ ] `HotspotAction` class
- [ ] Grid clustering
- [ ] Free energy calculation
- [ ] Hotspot ranking
- [ ] PDB output with B-factors

### 2.5 Tests
**File:** `tests/test_analysis.py`

- [ ] Test density calculation with mock trajectory
- [ ] Test residence time calculation
- [ ] Test hotspot detection
- [ ] Compare outputs with known-good results

---

## Phase 3: Project & Replica Management 🔲 NOT STARTED

**Status:** Not started  
**Estimated:** 1 day

### 3.1 Configuration
**File:** `pymdmix/project/config.py`

Port from: `settings.py`, `SettingsParser.py`, `Parsers.py`, `MDSettings.py` (~1538 lines total)

- [ ] `Config` dataclass for global settings
- [ ] `MDSettings` dataclass for simulation parameters
- [ ] YAML/JSON config loading (replace old .cfg format)
- [ ] Validation with sensible defaults
- [ ] Environment variable overrides (AMBERHOME, etc.)
- [ ] Path resolution utilities (from `tools.py`)

### 3.2 Replica
**File:** `pymdmix/project/replica.py`

Port from: `pyMDMix/Replicas.py` (960 lines)

- [ ] `Replica` dataclass
- [ ] State tracking (setup, running, complete, analyzed)
- [ ] Path management
- [ ] Topology/coordinate file references
- [ ] Serialization (JSON)

### 3.3 Project
**File:** `pymdmix/project/project.py`

Port from: `pyMDMix/Projects.py` (569 lines)

- [ ] `Project` class
- [ ] Multi-replica management
- [ ] Batch operations
- [ ] Project save/load

### 3.4 Tests
**File:** `tests/test_project.py`

- [ ] Test replica creation
- [ ] Test project save/load
- [ ] Test state transitions

---

## Phase 4: Structure Preparation 🔲 NOT STARTED

**Status:** Not started  
**Estimated:** 1-2 days

### 4.1 PDB Preparation
**File:** `pymdmix/setup/prepare.py`

Port from: `pyMDMix/AutoPrepare.py` (830 lines)

- [ ] `prepare_structure()` function
- [ ] Protonation (PDB2PQR integration or warn)
- [ ] ACE/NME capping
- [ ] Disulfide bond handling (CYS → CYX)
- [ ] Remove clashing waters
- [ ] Chain handling

**Key logic to preserve:**
- Cap fitting via superposition
- Disulfide detection by SG-SG distance
- Water clash removal

### 4.2 Solvation
**File:** `pymdmix/setup/solvate.py`

Port from: `pyMDMix/Systems.py` (518 lines)

- [ ] `solvate_structure()` function
- [ ] Solvent box creation
- [ ] Ion addition
- [ ] Box size calculation

### 4.3 Topology Generation
**File:** `pymdmix/setup/topology.py`

Port from: `pyMDMix/Amber.py` (1239 lines - partial), `pyMDMix/OFFManager.py` (357 lines)

- [ ] LEaP script generation
- [ ] OFF file handling (from OFFManager)
- [ ] Topology/coordinate output
- [ ] Force field handling
- [ ] Solvent box object file management

### 4.4 Tests
**File:** `tests/test_setup.py`

- [ ] Test capping
- [ ] Test disulfide detection
- [ ] Test solvation

---

## Phase 5: MD Engine Interfaces 🔲 NOT STARTED

**Status:** Not started  
**Estimated:** 1 day

### 5.1 Base Engine
**File:** `pymdmix/engines/base.py`

- [ ] `MDEngine` abstract base class
- [ ] Common interface for all engines
- [ ] Input file generation
- [ ] Job submission (queue scripts)

### 5.2 Amber Engine
**File:** `pymdmix/engines/amber.py`

Port from: `pyMDMix/Amber.py` (1239 lines)

- [ ] `AmberEngine` class
- [ ] Minimization input
- [ ] Equilibration input
- [ ] Production input
- [ ] Restraint handling

### 5.3 OpenMM Engine
**File:** `pymdmix/engines/openmm.py`

Port from: `pyMDMix/OpenMM.py` (615 lines)

- [ ] `OpenMMEngine` class
- [ ] System setup
- [ ] Simulation parameters

### 5.4 NAMD Engine (optional)
**File:** `pymdmix/engines/namd.py`

Port from: `pyMDMix/NAMD.py` (587 lines)

- [ ] `NAMDEngine` class
- [ ] Config file generation

### 5.5 Queue Scripts
**File:** `pymdmix/engines/queue.py`

Port from: `pyMDMix/QueueWriting.py` (178 lines)

- [ ] Template-based queue script generation
- [ ] Support: Slurm, SGE, PBS, local

---

## Phase 6: CLI & I/O 🔲 NOT STARTED

**Status:** Not started  
**Estimated:** 0.5 days

### 6.1 CLI
**File:** `pymdmix/cli.py`

- [ ] Click-based CLI
- [ ] `mdmix create` - create project/replica
- [ ] `mdmix setup` - prepare structures
- [ ] `mdmix analyze` - run analysis
- [ ] `mdmix info` - show project status

### 6.2 I/O Utilities
**File:** `pymdmix/io/`

- [ ] `readers.py` - file format readers
- [ ] `writers.py` - file format writers
- [ ] Grid format conversions (DX, MRC, CCP4)

### 6.3 Plotting
**File:** `pymdmix/io/plotting.py`

Port from: `pyMDMix/Plotter.py` (345 lines)

- [ ] Energy plot generation
- [ ] Density visualization helpers
- [ ] Matplotlib-based plotting utilities
- [ ] Optional (matplotlib not required for core)

---

## Phase 7: Data Migration & Testing 🔲 NOT STARTED

**Status:** Not started  
**Estimated:** 0.5 days

### 7.1 Solvent Library Migration
**Script:** `scripts/migrate_solvents.py`

- [ ] Convert old `.config` files to JSON
- [ ] Copy `.off` files to new location
- [ ] Validate all solvents load correctly

### 7.2 Integration Tests
**File:** `tests/test_integration.py`

- [ ] Full workflow: PDB → setup → (mock) → analysis
- [ ] Compare density outputs with old version
- [ ] Performance benchmarks

### 7.3 Documentation
**Directory:** `docs/`

- [ ] API documentation (docstrings → Sphinx)
- [ ] User guide
- [ ] Migration guide from v0.2.x
- [ ] Examples

---

## Phase 8: Polish & Release 🔲 NOT STARTED

**Status:** Not started  
**Estimated:** 0.5 days

### 8.1 Code Quality
- [ ] Run ruff/mypy, fix issues
- [ ] Ensure test coverage > 80%
- [ ] Review all TODOs

### 8.2 Packaging
- [ ] Verify pip install works
- [ ] Test conda install
- [ ] Create release notes

### 8.3 Release
- [ ] Tag v0.3.0
- [ ] Push to PyPI
- [ ] Update documentation site

---

## Timeline Summary

| Phase | Description | Estimate | Status |
|-------|-------------|----------|--------|
| 1 | Core Infrastructure | 0.5 days | ✅ Done |
| 2 | Analysis Actions | 1-2 days | 🔲 Next |
| 3 | Project Management | 1 day | 🔲 |
| 4 | Structure Preparation | 1-2 days | 🔲 |
| 5 | MD Engines | 1 day | 🔲 |
| 6 | CLI & I/O | 0.5 days | 🔲 |
| 7 | Migration & Testing | 0.5 days | 🔲 |
| 8 | Polish & Release | 0.5 days | 🔲 |
| **Total** | | **6-8 days** | |

---

## File Mapping: Old → New

### Core & Analysis
| Old File | Lines | New Location | Status |
|----------|-------|--------------|--------|
| `GridData.py` | 1341 | `core/grid.py` | ✅ |
| `GridsManager.py` | 968 | `core/grid.py` | ✅ (merged) |
| `PDB.py` | 362 | `core/structure.py` | ✅ |
| `Solvents.py` | 740 | `core/solvent.py` | ✅ |
| `Trajectory.py` | 176 | `core/trajectory.py` | ✅ |
| `NamdDCDParser.py` | 227 | `core/trajectory.py` | ✅ (MDAnalysis) |
| `Actions/Density.py` | 855 | `analysis/density.py` | 🔲 |
| `Actions/Residence.py` | 310 | `analysis/residence.py` | 🔲 |
| `HotSpotsManager.py` | 1050 | `analysis/hotspots.py` | 🔲 |

### Project Management
| Old File | Lines | New Location | Status |
|----------|-------|--------------|--------|
| `Projects.py` | 569 | `project/project.py` | 🔲 |
| `Replicas.py` | 960 | `project/replica.py` | 🔲 |
| `MDSettings.py` | 246 | `project/config.py` | 🔲 |
| `Parsers.py` | 476 | `project/config.py` | 🔲 |
| `SettingsParser.py` | 590 | `project/config.py` | 🔲 |
| `settings.py` | 226 | `project/config.py` | 🔲 |

### Setup & Engines
| Old File | Lines | New Location | Status |
|----------|-------|--------------|--------|
| `AutoPrepare.py` | 830 | `setup/prepare.py` | 🔲 |
| `Systems.py` | 518 | `setup/solvate.py` | 🔲 |
| `OFFManager.py` | 357 | `setup/topology.py` | 🔲 |
| `Amber.py` | 1239 | `engines/amber.py` | 🔲 |
| `OpenMM.py` | 615 | `engines/openmm.py` | 🔲 |
| `NAMD.py` | 587 | `engines/namd.py` | 🔲 |
| `QueueWriting.py` | 178 | `engines/queue.py` | 🔲 |

### I/O & Utilities
| Old File | Lines | New Location | Status |
|----------|-------|--------------|--------|
| `Plotter.py` | 345 | `io/plotting.py` | 🔲 |
| `tools.py` | 328 | Distributed / `utils.py` | 🔲 |
| `test.py` | 160 | `tests/` | 🔲 |

---

## Dependencies

### Required
```
numpy>=1.21
scipy>=1.7
parmed>=4.0
```

### Optional (full install)
```
MDAnalysis>=2.0
matplotlib>=3.5
netCDF4>=1.5
click>=8.0
```

### Development
```
pytest>=7.0
pytest-cov>=4.0
mypy>=1.0
ruff>=0.1
```
