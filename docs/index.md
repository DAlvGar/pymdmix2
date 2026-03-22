# pyMDMix Documentation

**pyMDMix** is a Python toolkit for identifying binding hotspots on protein surfaces using molecular dynamics simulations with organic solvent/water mixtures.

## Quick Start

```bash
# Install
pip install pymdmix

# Create a project
mdmix create project -n myproject

# Add a system
mdmix add system -f system.cfg

# Run analysis
mdmix analyze density bysolvent -s ETA
```

## Contents

- [Installation](installation.md)
- [Quickstart Tutorial](quickstart.md)
- [Configuration Files](configuration.md)
- [Solvents](solvents.md)
- [Analysis](analysis.md)
- [CLI Reference](cli.md)
- [API Reference](api.md)

## Overview

pyMDMix uses mixed-solvent molecular dynamics (MSMD) to probe protein binding sites. The approach involves:

1. **Solvation**: Place the protein in a box of organic solvent/water mixture
2. **Simulation**: Run MD simulations (typically 20ns × 3 replicas)
3. **Analysis**: Calculate solvent density maps and convert to binding free energies
4. **Hotspots**: Identify regions with favorable binding energy

### Key Concepts

- **Probes**: Functional groups (OH, NH, etc.) that sample the protein surface
- **Density Grids**: 3D maps of probe occupancy during simulation
- **Energy Grids**: Free energy maps (ΔG = -RT ln ρ/ρ₀)
- **Hotspots**: Clustered regions with favorable binding energy

## Supported Solvents

| Solvent | Probes | Description |
|---------|--------|-------------|
| ETA | OH, CT | Ethanol (hydroxyl + hydrophobic) |
| MAM | NH, CO | Methylacetamide (H-bond donor/acceptor) |
| WAT | O | Water reference |

## Requirements

- Python ≥ 3.10
- AmberTools (tleap, cpptraj)
- MD engine: Amber, OpenMM, or GROMACS

## License

pyMDMix is released under the GPL-3.0 License.

## Citation

If you use pyMDMix in your research, please cite:

```
Alvarez-Garcia D, Barril X. Molecular Simulations with Solvent Mixtures.
J. Chem. Theory Comput. 2014, 10, 6, 2608-2614.
```
