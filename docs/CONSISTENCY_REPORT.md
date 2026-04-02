# pyMDMix 2.0 Consistency Report

**Date:** 2026-03-22
**Status:** ✅ VERIFIED

## Module Coverage

| Legacy Module | New Location | Status |
|--------------|--------------|--------|
| `Align.py` | `analysis/align.py` | ✅ |
| `Amber.py` | `engines/amber.py` | ✅ |
| `Analysis.py` | `analysis/manager.py` | ✅ |
| `AutoPrepare.py` | `setup/prepare.py` | ✅ |
| `Browser.py` | `project/browser.py` | ✅ |
| `containers.py` | `core/containers.py` | ✅ |
| `Energy.py` | `analysis/energy.py` | ✅ |
| `Executor.py` | `engines/executor.py` | ✅ |
| `GridData.py` | `core/grid.py` | ✅ |
| `GridsManager.py` | `core/grid.py` | ✅ (merged) |
| `HotSpotsManager.py` | `analysis/hotspots.py` | ✅ |
| `MDSettings.py` | `project/settings.py` | ✅ |
| `NAMD.py` | `engines/namd.py` | ✅ (`NAMDEngine` + `NAMDWriter` + `NAMDChecker`) |
| `NamdDCDParser.py` | `io/dcd_parser.py` | ✅ |
| `OFFManager.py` | `io/off_manager.py` | ✅ |
| `OpenMM.py` | `engines/openmm.py` | ✅ (`OpenMMEngine` + `OpenMMWriter` + `OpenMMChecker`) |
| `Parsers.py` | `io/parsers.py` | ✅ |
| `PDB.py` | `core/structure.py` | ✅ (replaced by parmed) |
| `Plotter.py` | `io/plotting.py` | ✅ |
| `Projects.py` | `project/project.py` | ✅ |
| `QueueWriting.py` | `engines/queue.py` | ✅ |
| `Replicas.py` | `project/replica.py` | ✅ |
| `SettingsParser.py` | `utils/settings_parser.py` | ✅ |
| `Solvents.py` | `core/solvent.py` | ✅ |
| `Structures.py` | `core/structure.py` | ✅ |
| `Systems.py` | `core/system.py` | ✅ |
| `tools.py` | `utils/tools.py` | ✅ |
| `Trajectory.py` | `core/trajectory.py` | ✅ |

## CLI Coverage

| Command | Subcommands | Status |
|---------|-------------|--------|
| `create` | project, solvent | ✅ |
| `add` | system, replica, group | ✅ |
| `setup` | prepare, solvate | ✅ |
| `analyze` | align, density, energy, hotspots, filter-hotspots, residence | ✅ |
| `info` | project, systems, replicas, solvents, settings, analysis | ✅ |
| `plot` | rmsd, energy, density | ✅ |
| `queue` | generate, submit | ✅ |
| `tools` | diffgrids, sumgrids, avggrids, energy, grid-info, grid-math, convert | ✅ |

## Core Functionality Verified

### Structure Handling (Biskit replacement)
- [x] Load PDB files (parmed)
- [x] Protein/water/nucleic masks
- [x] Backbone/heavy atom masks
- [x] Hydrogen mask
- [x] Atom selection
- [x] Chain operations
- [x] Disulfide detection
- [x] Structure alignment (Kabsch)
- [x] Profile conversions

### MD Preparation
- [x] tleap integration
- [x] Solvation with organic solvents
- [x] Topology/coordinate generation
- [x] Minimization input generation

### Analysis
- [x] Trajectory alignment
- [x] Density grid calculation
- [x] Free energy conversion
- [x] Hotspot detection
- [x] Hotspot clustering (HotSpotSet)
- [x] Hotspot filtering

### Project Management
- [x] Project creation
- [x] System management
- [x] Replica management
- [x] MDSettings configuration

## Real Workflow Test

Tested on 2026-03-22 with `pep.pdb` (8-residue peptide):

1. **Structure loading:** ✅ 122 atoms, 8 residues
2. **tleap check:** ✅ No errors
3. **Ethanol solvation:** ✅ 89 ETA + 1425 WAT
4. **Topology generation:** ✅ 5198 atoms
5. **Minimization:** ✅ 100 steps in 2 seconds

## Test Count

- **698 tests passing** (updated after PR #1 and PR #2 merges)
- 8 skipped (optional dependencies)
- 2 pre-existing failures (unrelated to port work)

### Known Pre-existing Failures
| Test | Reason |
|------|--------|
| `test_gromacs.py::test_create_index_groups` | Asserts `"Non-Protein"` but GROMACS outputs `"non-Protein"` (capitalisation) |
| `test_project.py::TestConfig::test_config_validation` | Config validation not yet implemented; returns empty error list |

## Recent Changes (PR #1 and PR #2)

### PR #1 — Fix Public API Inconsistencies

- **`project/__init__.py`**: Fixed `MDSettings` export — now correctly re-exports `MDSettings` from `replica.py` (workflow-level: `nanos`, `temperature`, `restraint_mode`, …) instead of the Amber-specific `MDSettings` from `config.py`.
- **`analysis/__init__.py`**: Added missing public exports: `ActionsManager`, `Job`, `JobResult`, and `HotSpotSet` are now importable directly from `pymdmix.analysis`.
- **Deleted** stale editor artifact `pymdmix/analysis/manager.py.backup`.
- **`tests/test_utils.py`**: Suppressed spurious pytest collection warning on `test_root` helper.

### PR #2 — End-to-End and Integration Tests

Added 77 new tests covering full cross-module pipelines and real bundled data:

| File | Tests | Description |
|------|-------|-------------|
| `tests/test_e2e_workflow.py` | 40 | Synthetic-data e2e pipelines (density→hotspot, grid chain, project lifecycle, etc.) |
| `tests/test_integration_real_data.py` | 37 | Real bundled data (`pep.pdb`, `pep.prmtop`, solvent JSONs, OFF files, legacy configs) |

## Dependencies

- Python 3.10+
- parmed 4.3.1
- numpy, scipy
- click (CLI)
- AmberTools (tleap, cpptraj)
- Optional: MDAnalysis, matplotlib
