# pyMDMix

A Python toolkit for molecular dynamics simulations with organic solvent mixtures, used in drug discovery for identifying binding hotspots on protein surfaces.

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)

## Overview

pyMDMix enables drug discovery researchers to:

- **Identify binding hotspots** on protein surfaces using organic solvent mixtures
- **Run MD simulations** with various solvent mixtures (ethanol, acetonitrile, etc.)
- **Analyze probe density** to map favorable binding regions
- **Calculate free energies** at potential binding sites

### Key Features

- Modern Python 3.10+ codebase with type hints
- Built on parmed and MDAnalysis for robust MD handling
- Support for multiple MD engines (Amber, OpenMM, NAMD)
- JSON-based solvent definitions for easy customization
- CLI and Python API interfaces
- Integration tests with 190+ passing tests

## Installation

### Requirements

- Python 3.10+
- AmberTools (AMBER MD suite)

### Using pip

```bash
# Clone the repository
git clone https://github.com/DAlvGar/pyMDmix.git
cd pyMDmix

# Install in development mode
pip install -e ".[dev]"
```

### Using conda

```bash
# Create environment with AmberTools
conda create -n mdmix python=3.10 ambertools -c conda-forge
conda activate mdmix

# Install pyMDMix
pip install -e .
```

## Quick Start

### Command Line

```bash
# List available solvents
mdmix info --solvents

# Prepare a structure
mdmix setup prepare structure.pdb -o prepared.pdb

# Create a project
mdmix create project -n my_project -s ETA

# Run density analysis
mdmix analyze density -p my_project
```

### Python API

```python
from pymdmix.core.solvent import SolventLibrary
from pymdmix.core.grid import Grid
from pymdmix.project import Project, Replica

# Load available solvents
library = SolventLibrary()
ethanol = library.get("ETA")

# Create a project
project = Project(name="binding_study", path="./my_project")

# Add a replica
replica = Replica(name="rep1", solvent="ETA")
project.add_replica(replica)

# Analyze density
from pymdmix.analysis import DensityAction
action = DensityAction()
result = action.run(trajectory, probes=ethanol.probes)
```

## Solvent Mixtures

Pre-defined solvent mixtures for different probe types:

| Solvent | Description | Probe Types |
|---------|-------------|-------------|
| ETA | 20% Ethanol | Hydrophobic, H-bond donor/acceptor |
| MAM | 20% Methylamine | Positive charge |
| ANT | 20% Acetonitrile | H-bond acceptor |
| ION | Ammonium acetate | Positive/negative charge |
| ISO | 20% Isopropanol | Hydrophobic |

## Documentation

- [Implementation Plan](docs/IMPLEMENTATION_PLAN.md)
- [Progress](docs/PROGRESS.md)
- [Rewrite Plan](docs/REWRITE_PLAN.md)

## Development

```bash
# Install dev dependencies
pip install -e ".[dev]"

# Run tests
pytest tests/

# Run linting
ruff check pymdmix/

# Type checking
mypy pymdmix/
```

## Project Structure

```
pymdmix/
├── core/           # Core modules (grid, trajectory, structure, solvent)
├── analysis/       # Analysis actions (density, residence, hotspots)
├── project/        # Project management (config, replica, project)
├── setup/          # Structure preparation and solvation
├── engines/        # MD engine interfaces (amber, queue)
├── io/             # I/O utilities (grid formats, plotting)
├── cli.py          # Command-line interface
└── data/solvents/  # Solvent definitions (JSON + OFF files)
```

## License

GNU General Public License v3.0 - see [LICENSE](LICENSE) for details.

## Citation

If you use pyMDMix in your research, please cite:

> Alvarez-Garcia, D. et al. pyMDMix: A Python toolkit for mixed solvent molecular dynamics. 

## Author

- **Daniel Alvarez-Garcia** - [algarcia.daniel@gmail.com](mailto:algarcia.daniel@gmail.com)

## Acknowledgments

- Original pyMDMix development supported by various research grants
- Built on the excellent work of the parmed and MDAnalysis teams
