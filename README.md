# pyMDMix2

A modern Python toolkit for molecular dynamics simulations with organic solvent mixtures, used in drug discovery for identifying binding hotspots on protein surfaces.

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![Tests](https://img.shields.io/badge/tests-739%20passed-brightgreen.svg)]()

## Overview

pyMDMix2 enables drug discovery researchers to:

- **Identify binding hotspots** on protein surfaces using organic solvent mixtures
- **Run MD simulations** with various solvent mixtures (ethanol, acetamide, acetonitrile, etc.)
- **Analyze probe density** to map favorable binding regions
- **Calculate free energies** at potential binding sites
- **Detect and cluster hotspots** for fragment-based drug design

### Key Features

- Modern Python 3.10+ codebase with full type hints
- Clean rewrite of the original pyMDMix (Python 2.7)
- Built on parmed and MDAnalysis for robust MD handling
- Support for multiple MD engines (Amber, OpenMM, NAMD, GROMACS)
- **Single-file project initialization** from a combined `project.cfg`
- CLI and Python API interfaces
- Comprehensive test suite (700+ tests)

## Installation

```bash
# Clone the repository
git clone https://github.com/DAlvGar/pymdmix2.git
cd pymdmix2

# Install with uv (recommended)
uv sync --all-extras

# Or with pip
pip install -e ".[dev]"
```

**Requirements:** Python 3.10+, AmberTools (for solvation and analysis)

## Quick Start: Project from a Single Config File

The primary workflow uses a single annotated `project.cfg` file that defines
the system and all MD settings in one place — the same approach as the original
pyMDMix.

### Step 1 — Get an annotated template

```bash
pymdmix create template               # saves project.cfg in current folder
# or
pymdmix create template -o my_sim.cfg
```

### Step 2 — Edit the template

```ini
[SYSTEM]
NAME = MyProtein
OFF  = myprotein.off   # or PDB = myprotein.pdb

[MDSETTINGS]
SOLVENTS = ETA, MAM, WAT
NREPL    = 3           # 3 replicas per solvent
NANOS    = 20          # 20 ns production
TEMP     = 300         # 300 K
RESTR    = FREE        # no positional restraints
```

### Step 3 — Create the project

```bash
pymdmix create project -n myproject -f project.cfg
```

This single command:
1. Creates the project folder structure
2. Registers the system
3. Creates all replicas (here: 3 × 3 solvents = 9 replicas)

### Step 4 — Run simulations

```bash
# Check what was created
pymdmix info project

# Run on your HPC cluster
cd MyProtein_ETA_1
cat COMMANDS.sh        # review generated commands
./COMMANDS.sh          # or submit queue.sh to your scheduler
```

### Step 5 — Analyze results

```bash
# Align trajectories
pymdmix analyze align all

# Calculate density grids
pymdmix analyze density all

# Convert density to free energy
pymdmix analyze energy all

# Find binding hotspots
pymdmix analyze hotspots all --threshold -0.5
```

---

## Alternative: Step-by-step Setup

If you prefer to build the project incrementally:

```bash
# 1. Create empty project
pymdmix create project -n myproject

# 2. Add the system
cat > system.cfg << EOF
[SYSTEM]
NAME = MyProtein
OFF  = myprotein.off
EOF
pymdmix add system -f system.cfg

# 3. Add replicas using an MDSETTINGS file
cat > settings.cfg << EOF
[MDSETTINGS]
SOLVENTS = ETA, MAM
NREPL    = 3
NANOS    = 20
RESTR    = FREE
EOF
pymdmix add replica -f settings.cfg
```

---

## Configuration Reference

### Project config file (`project.cfg`)

| Section | Key | Default | Description |
|---------|-----|---------|-------------|
| `[SYSTEM]` | `NAME` | — | System identifier (required) |
| `[SYSTEM]` | `OFF` | — | Amber Object File path (or `PDB`) |
| `[SYSTEM]` | `EXTRARES` | — | Non-standard residue names |
| `[MDSETTINGS]` | `SOLVENTS` | — | Comma-separated solvent list (required) |
| `[MDSETTINGS]` | `NREPL` | 1 | Replicas per solvent |
| `[MDSETTINGS]` | `NANOS` | 20 | Production nanoseconds |
| `[MDSETTINGS]` | `TEMP` | 300 | Temperature (K) |
| `[MDSETTINGS]` | `RESTR` | `FREE` | Restraint mode (`FREE`, `HA`, `BB`) |
| `[MDSETTINGS]` | `FORCE` | 0.0 | Restraint force (kcal/mol·Å²) |

### Restraint Modes

| Mode | Atoms | Amber Mask | Typical Force |
|------|-------|------------|---------------|
| `FREE` | None | - | 0 |
| `BB` | Backbone | `@CA,C,N,O` | 0.1 - 5.0 |
| `HA` | Heavy atoms | `!@H=` | 0.01 - 1.0 |

### Solvent Mixtures

| Solvent | Description | Probes |
|---------|-------------|--------|
| ETA | 20% Ethanol / 80% Water | OH, CT, WAT |
| MAM | 20% Acetamide / 80% Water | N, O, CT, WAT |
| ANT | 20% Acetonitrile / 80% Water | N, C, WAT |
| WAT | Pure water | O |

---

## Python API

```python
from pathlib import Path
from pymdmix.io.parsers import parse_project_config
from pymdmix.project import Project

# Load a full project config file
cfg = parse_project_config("project.cfg")
print(cfg.system.name)        # → 'MyProtein'
print(len(cfg.settings))      # → 9  (3 solvents × 3 replicas)

# Or work with a project directly
project = Project.load("myproject")
for replica in project.replicas:
    print(f"{replica.name}: {replica.settings.nanos} ns")
```

---

## Development Setup

```bash
# Sync all dev + optional deps
uv sync --all-extras

# Install pre-commit hooks
uv run pre-commit install

# Run tests
uv run pytest tests/

# Lint / format
uv run ruff check pymdmix
uv run ruff format pymdmix

# Type-check
uv run mypy pymdmix
```

---

## Project Structure

```
pymdmix/
├── core/           # Grid, trajectory, structure, solvent
├── analysis/       # Density, energy, hotspots, alignment
├── project/        # Replica, settings, config
├── engines/        # Amber, OpenMM, NAMD, GROMACS
├── io/             # Parsers, grid I/O, DCD reader
├── cli.py          # Command-line interface
└── data/
    ├── solvents/   # Solvent definitions (JSON)
    ├── defaults/   # Default settings files
    └── templates/  # Config file templates
```

## Documentation

See `docs/` and `examples/` for detailed guides.

- [Quick Start](docs/quickstart.md)
- [User Guide: Projects](docs/user-guide/projects.md)
- [User Guide: MD Settings](docs/user-guide/md-settings.md)
- [Tutorials](docs/tutorials/toy-project-setup.md)

## License

MIT License — see [LICENSE](LICENSE)

## Citation

> Alvarez-Garcia, D., & Barril, X. (2014). Molecular simulations with solvent competition quantify water displaceability and provide accurate interaction maps of protein binding sites. *J. Med. Chem.*, 57(20), 8530-8539.

## Author

**Daniel Alvarez-Garcia** — [algarcia.daniel@gmail.com](mailto:algarcia.daniel@gmail.com)
