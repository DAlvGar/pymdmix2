# pyMDMix2

A modern Python toolkit for molecular dynamics simulations with organic solvent mixtures, used in drug discovery for identifying binding hotspots on protein surfaces.

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![Tests](https://img.shields.io/badge/tests-621%20passed-brightgreen.svg)]()

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
- JSON/TOML-based configuration for easy customization
- CLI and Python API interfaces
- Comprehensive test suite (621+ tests)

## Installation

```bash
# Clone the repository
git clone https://github.com/DAlvGar/pymdmix2.git
cd pymdmix2

# Install with uv (recommended)
uv sync --all-extras

# Activate the managed venv (optional — uv run works without it)
source .venv/bin/activate
```

**Requirements:** Python 3.10+, AmberTools (for solvation)

## Development Setup

```bash
# Sync all dev + optional deps
uv sync --all-extras

# Install pre-commit hooks (runs on every git commit)
uv run pre-commit install

# Run hooks manually across all files
uv run pre-commit run --all-files

# Run tests
uv run pytest tests/

# Lint / format
uv run ruff check pymdmix
uv run ruff format pymdmix

# Type-check
uv run mypy pymdmix

# Add a dependency
uv add <package>           # runtime
uv add --dev <package>     # dev-only
```

## Quick Start: Creating a Project

### Python API

```python
from pathlib import Path
from pymdmix.core.system import System
from pymdmix.project.replica import Replica, MDSettings

# =============================================================================
# 1. Load your system from an Amber OFF file
# =============================================================================
system = System.from_off("protein.off", unit_name="mol", name="myprotein")

# =============================================================================
# 2. Solvate with organic solvent mixture
# =============================================================================
solvated = system.solvate("ETA")  # 20% Ethanol / 80% Water

# =============================================================================
# 3. Configure MD settings
# =============================================================================
settings = MDSettings(
    nanos=50,                   # 50 ns production
    temperature=300.0,          # 300 K
    restraint_mode="HA",        # Restrain heavy atoms
    restraint_force=0.01,       # Soft: 0.01 kcal/mol·Å²
    align_mask="@CA,C,N",       # Alignment mask for analysis
)

# =============================================================================
# 4. Create replicas
# =============================================================================
replicas = []
for i in range(1, 3):  # 2 replicates
    replica = Replica(
        name=f"ETA_{i}",
        solvent="ETA",
        path=Path(f"replicas/ETA_{i}"),
        settings=settings,
    )

    # Create folder structure with topology/coordinates
    replica.create_folder(system=solvated)

    # Generate Amber input files (min, eq, prod)
    replica.create_md_input()

    replicas.append(replica)
    print(f"Created: {replica.path}")

# =============================================================================
# 5. Run MD, then analyze
# =============================================================================
# After running ./COMMANDS.sh in each replica folder:

from pymdmix.core.grid import Grid

# Load density and convert to free energy
density = Grid.read_dx("replicas/ETA_1/density/OH_density.dx")
energy = density.to_free_energy(temperature=300.0)
energy.write_dx("OH_energy.dx")
print(f"Best binding: {energy.data.min():.2f} kcal/mol")
```

### CLI Workflow

```bash
# 1. Prepare structure
pymdmix setup prepare protein.pdb -o protein_clean.pdb

# 2. Create project
pymdmix create project -n my_project

# 3. Add system (system.cfg)
cat > system.cfg << EOF
[SYSTEM]
NAME = myprotein
PDB = protein_clean.pdb
EOF
pymdmix add system -f system.cfg

# 4. Add replicas with restraints (replica.cfg)
cat > replica.cfg << EOF
[REPLICA]
SYSTEM = myprotein
SOLVENT = ETA
NANOS = 50
RESTRMODE = HA
RESTRFORCE = 0.01
EOF
pymdmix add replica -f replica.cfg --count 2

# 5. After MD - run analysis
pymdmix analyze align all
pymdmix analyze density all
pymdmix analyze energy all
pymdmix analyze hotspots all --threshold -0.5
```

## Configuration Reference

### Restraint Modes

| Mode | Atoms | Amber Mask | Typical Force |
|------|-------|------------|---------------|
| `FREE` | None | - | 0 |
| `BB` | Backbone | `@CA,C,N,O` | 0.1 - 5.0 |
| `HA` | Heavy atoms | `!@H=` | 0.01 - 1.0 |
| `CUSTOM` | User mask | `RESTRAINMASK` | varies |

### Solvent Mixtures

| Solvent | Description | Probes |
|---------|-------------|--------|
| ETA | 20% Ethanol / 80% Water | OH, CT, WAT |
| MAM | 20% Acetamide / 80% Water | N, O, CT, WAT |
| ANT | 20% Acetonitrile / 80% Water | N, C, WAT |
| WAT | Pure water | O |

## Project Structure

```
pymdmix/
├── core/           # Grid, trajectory, structure, solvent
├── analysis/       # Density, energy, hotspots, alignment
├── project/        # Replica, settings, config
├── engines/        # Amber, OpenMM, NAMD, GROMACS
├── cli.py          # Command-line interface
└── data/solvents/  # Solvent definitions
```

## Documentation

See `docs/` and `examples/` for detailed guides and configuration templates.

## License

MIT License - see [LICENSE](LICENSE)

## Citation

> Alvarez-Garcia, D., & Barril, X. (2014). Molecular simulations with solvent competition quantify water displaceability and provide accurate interaction maps of protein binding sites. *J. Med. Chem.*, 57(20), 8530-8539.

## Author

**Daniel Alvarez-Garcia** - [algarcia.daniel@gmail.com](mailto:algarcia.daniel@gmail.com)
