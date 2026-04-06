# Running the Test Suite

pyMDMix2 has three test tiers. Each tier requires progressively more external
resources — from nothing, to Git LFS binary files, to a full AmberTools
installation.

---

## Overview

| Tier | Marker | Requires | Count |
|------|--------|----------|-------|
| [Standard](#1-standard-tests-no-marker) | *(none)* | Python + optional extras (MDAnalysis, netCDF4…) | ~750 |
| [Real data](#2-real-data-tests-real_data) | `real_data` | Git LFS binary files in `tests/data/` | ~68 |
| [AmberTools](#3-ambertools-tests) | `ambertools` / `cpptraj` | AmberTools installation (tleap **and** cpptraj) | ~40 |

> **Note on markers:** `cpptraj` is part of AmberTools — it ships in the same
> conda package and the same Docker image.  The `cpptraj` marker is a
> sub-classification of `ambertools` tests (it tags tests that specifically
> call the cpptraj binary).  Both markers skip under the **same** condition:
> AmberTools not installed.  There is no scenario where `ambertools` tests pass
> but `cpptraj` tests skip.

All three tiers auto-skip gracefully when their prerequisites are missing —
no configuration needed.

---

## Prerequisites

### Python environment

```bash
git clone https://github.com/DAlvGar/pymdmix2.git
cd pymdmix2
uv sync --all-extras    # install all optional dependencies
```

### Git LFS (for `real_data` tier)

```bash
git lfs install         # one-time setup
git lfs pull            # download binary test data (~20 MB)
```

Files fetched:

| File | Size | Description |
|------|------|-------------|
| `tests/data/amber/traj.nc` | 1.7 MB | 20-frame Amber NetCDF trajectory (7079 atoms) |
| `tests/data/amber/pep_WAT_WAT_1.pdb` | 576 KB | Pre-solvated peptide reference structure |
| `tests/data/grids/ETA_CT.dx` | 17 MB | Pre-computed ETA CT-probe density grid (132³, 0.5 Å) |

### AmberTools (for `ambertools` and `cpptraj` tiers)

Both `ambertools` and `cpptraj` markers require a single AmberTools
installation — cpptraj ships as part of AmberTools.

Install from [ambermd.org](https://ambermd.org) or via conda-forge:

```bash
conda install -c conda-forge ambertools
```

Then set `AMBERHOME`:

```bash
export AMBERHOME=/opt/conda    # adjust to your install path
```

Or, use the pre-built Docker image (see [Docker workflow](#docker-workflow)).

---

## 1. Standard tests (no marker)

These tests use only synthetic data generated in-memory. They run without any
binary files or external tools.

```bash
# Run all standard tests
uv run pytest tests/

# Run specific test files
uv run pytest tests/test_grid.py tests/test_trajectory.py -v

# Exclude the optional tiers explicitly
uv run pytest tests/ -m "not ambertools and not cpptraj and not real_data"
```

---

## 2. Real-data tests (`real_data`)

These tests load the actual binary files from `tests/data/` and exercise the
full I/O and analysis pipeline.  They are automatically skipped when the LFS
objects are not checked out.

```bash
# Pull the LFS data first
git lfs pull

# Run all real-data tests
uv run pytest -m real_data -v tests/

# Run individual test modules
uv run pytest tests/test_trajectory_real.py -v    # trajectory reading
uv run pytest tests/test_align_real.py     -v    # MDAnalysis alignment
uv run pytest tests/test_density_real.py   -v    # density calculation (Python)
uv run pytest tests/test_grid_real.py      -v    # Grid I/O and free energy
```

### What is tested

| File | What it tests |
|------|--------------|
| `test_trajectory_real.py` | `open_trajectory` with both MDAnalysis and native Amber NetCDF backends; `FrameSliceReader` |
| `test_align_real.py` | `align_trajectory` — MDAnalysis alignment; RMSD regression guard |
| `test_density_real.py` | `DensityAction` (pure Python) and `CpptrajDensityAction` (cpptraj-gated) with real trajectory |
| `test_grid_real.py` | Load `ETA_CT.dx`; DX round-trip; free energy conversion; `HotspotAction` |

> **Tip:** You can visually inspect the generated DX grids in VMD, PyMOL, or UCSF
> ChimeraX to verify that probe density appears in chemically sensible locations on
> the peptide surface.

---

## 3. AmberTools tests

Both markers below require the **same** thing: a working AmberTools
installation.  `cpptraj` is not a separate tool — it is bundled in AmberTools.

### `ambertools` marker

Tags tests that use tleap (LEaP) for solvation and project-creation workflows.

```bash
# Directly, if AmberTools is installed
uv run pytest -m ambertools -v tests/test_ambertools_e2e.py

# Or with the Docker helper
./scripts/run_ambertools_tests.sh
```

### `cpptraj` marker

Tags tests that specifically call the cpptraj binary for density calculations.
Because cpptraj ships with AmberTools, these tests skip under the **same**
condition as `ambertools` tests.

```bash
uv run pytest -m cpptraj -v tests/
```

### Combined: cpptraj + real trajectory data

These tests run cpptraj against the real 20-frame `traj.nc` trajectory. They
require AmberTools **and** the Git LFS data.

```bash
# Pull data first
git lfs pull

# Run combined tests
uv run pytest -m "cpptraj and real_data" -v tests/
```

---

## Docker workflow

The Dockerfile bundles AmberTools (including cpptraj) + the full pymdmix
package. Use it to run AmberTools and real-data tests without modifying your
local environment.

### Build the image

```bash
./scripts/build_docker.sh
# or force a clean rebuild:
REBUILD=1 ./scripts/build_docker.sh
```

### Run all AmberTools tests

```bash
./scripts/run_ambertools_tests.sh
```

### Run AmberTools + real-data tests

```bash
# Pull LFS data so it is available when the container mounts tests/
git lfs pull

MARKER="ambertools or real_data" ./scripts/run_ambertools_tests.sh
```

### Run cpptraj tests on the real trajectory

```bash
git lfs pull
MARKER="cpptraj and real_data" ./scripts/run_ambertools_tests.sh
```

### Script environment variables

| Variable | Default | Description |
|----------|---------|-------------|
| `IMAGE_NAME` | `pymdmix` | Docker image name |
| `IMAGE_TAG` | `latest` | Docker image tag |
| `REBUILD` | `0` | Set to `1` to force a fresh `docker build` |
| `PYTEST_ARGS` | `-v --tb=short` | Extra arguments passed to `pytest` |
| `MARKER` | `ambertools` | pytest marker expression |

---

## Running a combined check

To run the entire test suite (all tiers) in one command, ensure both LFS data
and AmberTools are available, then:

```bash
git lfs pull
uv run pytest tests/ -v
```

Tests that lack their prerequisites are reported as `SKIPPED` (not failures).

---

## Test data layout

```
tests/
├── conftest.py              # shared fixtures and auto-skip hooks
├── data/
│   ├── pep/                 # dry peptide (pep.pdb, pep.off, pep.prmtop, pep.prmcrd)
│   ├── amber/               # real solvated-system files (Git LFS)
│   │   ├── traj.nc          # 20-frame Amber NetCDF trajectory (1.7 MB)
│   │   └── pep_WAT_WAT_1.pdb# solvated reference structure (576 KB)
│   └── grids/               # pre-computed density grids (Git LFS)
│       └── ETA_CT.dx        # ETA CT-probe density (17 MB, 132³ grid)
└── test_*.py                # test modules
```

Binary files in `tests/data/` are tracked via Git LFS (see `.gitattributes`).
Run `git lfs pull` to download them before running `real_data` or combined
`cpptraj and real_data` tests.

---

## Adding new tests

### Standard tests
Write a new `tests/test_<name>.py`.  Use synthetic data from `conftest.py`
fixtures (`sample_coordinates`, `sample_trajectory_coords`, etc.).

### Real-data tests
1. Add `@pytest.mark.real_data` to the test class or function.
2. Use the session-scoped path fixtures from `conftest.py`:
   - `pep_pdb_path`, `pep_prmtop_path`, `pep_prmcrd_path`, `pep_off_path`
   - `solvated_pep_pdb_path`, `amber_traj_nc_path`, `eta_ct_dx_path`

### AmberTools tests
1. Add `@pytest.mark.ambertools` to tests that need tleap/LEaP.
2. Add `@pytest.mark.cpptraj` (in addition to `@pytest.mark.ambertools`) to
   tests that specifically call the cpptraj binary.  Both markers skip under
   the same condition — they just document which tool is exercised.
3. Use the `solvated_pep` session fixture from `test_ambertools_e2e.py` if you
   need a solvated topology generated by tleap.
4. If your test needs both cpptraj and real data, add both markers:
   `@pytest.mark.cpptraj` + `@pytest.mark.real_data`.
