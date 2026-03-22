# pyMDMix2

A modern Python toolkit for molecular dynamics simulations with organic solvent mixtures, used in drug discovery for identifying binding hotspots on protein surfaces.

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![Tests](https://img.shields.io/badge/tests-627%20passed-brightgreen.svg)]()

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
- Comprehensive test suite (627 tests)

## Installation

### Requirements

- Python 3.10+
- AmberTools (for solvation and MD prep)

### Using pip

```bash
# Clone the repository
git clone https://github.com/DAlvGar/pymdmix2.git
cd pymdmix2

# Install
pip install -e .

# With optional dependencies
pip install -e ".[full]"

# With development tools
pip install -e ".[dev]"
```

### Using conda

```bash
# Create environment with AmberTools
conda create -n mdmix python=3.10 ambertools -c conda-forge
conda activate mdmix

# Install pyMDMix2
pip install -e .
```

## Quick Start

### Command Line

```bash
# List available solvents
mdmix info solvents

# Prepare a structure
mdmix setup prepare protein.pdb -o prepared.pdb

# Solvate with ethanol mixture
mdmix setup solvate prepared.pdb -s ETA -o solvated

# Create a project
mdmix create project -n my_project

# Add replicas
mdmix add replica -f replica.cfg --count 3

# Run analysis (after MD)
mdmix analyze align all
mdmix analyze density all
mdmix analyze energy all
mdmix analyze hotspots all --threshold -0.5
```

### Replica Configuration

```ini
[REPLICA]
SYSTEM = myprotein
SOLVENT = ETA
NANOS = 50
RESTRMODE = HA          # FREE | BB | HA | CUSTOM
RESTRFORCE = 0.01       # kcal/mol·Å² (soft restraints)
```

### Python API

```python
from pymdmix.core.solvent import SolventLibrary
from pymdmix.core.grid import Grid
from pymdmix.project.replica import Replica, MDSettings

# Load available solvents
library = SolventLibrary()
ethanol = library.get("ETA")
print(ethanol.probes)  # [OH, CT, WAT]

# Configure MD settings with soft restraints
settings = MDSettings(
    nanos=50,
    temperature=300.0,
    restraint_mode="HA",
    restraint_force=0.01,
)

# Create a replica
replica = Replica(
    name="ETA_1",
    solvent="ETA",
    settings=settings,
)

# Work with grids
grid = Grid.read_dx("density.dx")
energy = grid.to_free_energy(temperature=300.0)
energy.write_dx("energy.dx")
```

## Solvent Mixtures

Pre-defined solvent mixtures with different probe types:

| Solvent | Description | Probes |
|---------|-------------|--------|
| ETA | 20% Ethanol / 80% Water | OH (H-bond), CT (hydrophobic), WAT |
| MAM | 20% Acetamide / 80% Water | N (donor), O (acceptor), CT, WAT |
| ANT | 20% Acetonitrile / 80% Water | N (acceptor), C (hydrophobic), WAT |
| WAT | Pure water | O |

## Restraint Modes

| Mode | Atoms Restrained | Typical Force |
|------|------------------|---------------|
| `FREE` | None | 0 |
| `BB` | Backbone (CA,C,N,O) | 0.1 - 5.0 |
| `HA` | Heavy atoms (!@H=) | 0.01 - 1.0 |
| `CUSTOM` | User mask | varies |

## Documentation

Full documentation available in the `docs/` directory:

- [User Guide](docs/user-guide/) - Getting started and workflows
- [API Reference](docs/api/) - Python API documentation
- [Tutorials](docs/tutorials/) - Step-by-step examples

## Development

```bash
# Install dev dependencies
pip install -e ".[dev]"

# Run tests
pytest tests/ -q

# Run linting
ruff check pymdmix/

# Type checking
mypy pymdmix/
```

## Project Structure

```
pymdmix/
├── core/           # Core modules (grid, trajectory, structure, solvent)
├── analysis/       # Analysis actions (density, energy, hotspots, align)
├── project/        # Project management (replica, config, settings)
├── setup/          # Structure preparation and solvation
├── engines/        # MD engine interfaces (amber, openmm, namd, gromacs)
├── io/             # I/O utilities (grid formats, parsers)
├── cli.py          # Command-line interface
└── data/solvents/  # Solvent definitions (JSON + OFF files)
```

## License

MIT License - see [LICENSE](LICENSE) for details.

## Citation

If you use pyMDMix2 in your research, please cite:

> Alvarez-Garcia, D., & Barril, X. (2014). Molecular simulations with solvent competition quantify water displaceability and provide accurate interaction maps of protein binding sites. *Journal of Medicinal Chemistry*, 57(20), 8530-8539.

## Author

**Daniel Alvarez-Garcia** - [algarcia.daniel@gmail.com](mailto:algarcia.daniel@gmail.com)

## See Also

- [Original pyMDMix](https://github.com/DAlvGar/pyMDmix) (deprecated, Python 2.7)
