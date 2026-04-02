# pyMDMix Documentation

**pyMDMix** is a Python toolkit for molecular dynamics simulations with organic solvent mixtures, used in drug discovery for identifying binding hotspots on protein surfaces.

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## What is MDMix?

The MDMix method uses molecular dynamics simulations with organic solvent/water mixtures to probe protein surfaces and identify favorable binding regions (hotspots). By tracking the density of organic probe molecules around the protein, we can:

- **Identify binding hotspots** on protein surfaces
- **Map druggable pockets** with atomic resolution
- **Calculate binding free energies** at potential interaction sites
- **Guide drug design** by revealing preferred chemical interactions

## Key Features

- 🧪 **Multiple Solvent Mixtures**: Pre-equilibrated boxes for ethanol, acetonitrile, isopropanol, acetamide, and more
- 🔬 **Multiple MD Engines**: Support for Amber, OpenMM, NAMD, and GROMACS
- 📊 **Comprehensive Analysis**: Density grids, energy maps, hotspot detection, residence times
- 🖥️ **CLI & Python API**: Both command-line tools and programmatic access
- 📁 **Project Management**: Organize multiple systems and replicas efficiently

## Quick Start

```bash
# Install
git clone https://github.com/DAlvGar/pymdmix2.git
cd pymdmix2
uv sync --all-extras

# Create a new project
pymdmix create project myproject

# Check available solvents
pymdmix info solvents

# Get help
pymdmix --help
```

## Documentation

### Getting Started
- [Installation](installation.md) - Install pyMDMix and dependencies
- [Quick Start](quickstart.md) - Your first MDMix simulation

### User Guide
- [System Preparation](user-guide/system-preparation.md) - Prepare your protein structure
- [Projects & Systems](user-guide/projects.md) - Manage MDMix projects
- [Creating Replicas](user-guide/replicas.md) - Set up simulation replicas
- [Solvent Mixtures](user-guide/solvents.md) - Available solvents and custom mixtures
- [MD Settings](user-guide/md-settings.md) - Simulation parameters
- [Queue Scripts](user-guide/queue-scripts.md) - HPC cluster submission

### Analysis
- [Analysis Overview](analysis/overview.md) - Analysis workflow
- [Trajectory Alignment](analysis/alignment.md) - Align and image trajectories
- [Density Calculation](analysis/density.md) - Compute probe density grids
- [Energy Conversion](analysis/energy.md) - Convert density to free energy
- [Hotspot Detection](analysis/hotspots.md) - Identify binding hotspots
- [Residence Analysis](analysis/residence.md) - Solvent exchange kinetics

### Tutorials
- [Tutorial: Toy Project Setup](tutorials/toy-project-setup.md) - Step-by-step project creation
- [Tutorial: Running Simulations](tutorials/toy-project-simulation.md) - Execute and analyze

### Reference
- [CLI Reference](api/cli-reference.md) - Command-line interface
- [Tools](tools.md) - Additional utilities

## Workflow Overview

```
┌─────────────────────────────────────────────────────────────────┐
│                     pyMDMix Workflow                            │
├─────────────────────────────────────────────────────────────────┤
│                                                                 │
│  1. PREPARATION                                                 │
│     ├── Prepare structure (PDB/OFF)                            │
│     ├── Create project                                          │
│     └── Add system                                              │
│                                                                 │
│  2. REPLICA SETUP                                               │
│     ├── Choose solvent mixture                                  │
│     ├── Configure MD settings                                   │
│     └── Create replicas                                         │
│                                                                 │
│  3. SIMULATION                                                  │
│     ├── Run MD (your HPC cluster)                              │
│     └── Return trajectories to project                          │
│                                                                 │
│  4. ANALYSIS                                                    │
│     ├── Align trajectories                                      │
│     ├── Calculate density grids                                 │
│     ├── Convert to energy maps                                  │
│     └── Identify hotspots                                       │
│                                                                 │
└─────────────────────────────────────────────────────────────────┘
```

## Architecture

pyMDMix is organized into several core modules:

| Module | Purpose |
|--------|---------|
| `pymdmix.core` | Grid handling, trajectory I/O, solvent library |
| `pymdmix.analysis` | Density, residence, hotspot analysis actions |
| `pymdmix.project` | Project, replica, and configuration management |
| `pymdmix.setup` | Structure preparation and solvation |
| `pymdmix.engines` | MD engine interfaces (Amber, OpenMM, NAMD, GROMACS) |
| `pymdmix.io` | File format handling and plotting |
| `pymdmix.cli` | Command-line interface |

## Citation

If you use pyMDMix in your research, please cite:

> Alvarez-Garcia, D., & Barril, X. (2014). Molecular simulations with solvent competition quantify water displaceability and provide accurate interaction maps of protein binding sites. *Journal of Medicinal Chemistry*, 57(20), 8530-8539.

## License

pyMDMix is released under the [MIT License](https://opensource.org/licenses/MIT).

## Support

- **Issues**: [GitHub Issues](https://github.com/DAlvGar/pymdmix2/issues)
- **Author**: Daniel Alvarez-Garcia
