# Installation

## Requirements

- **Python 3.10+**
- **AmberTools** (for LEaP and cpptraj)
- Optional: OpenMM, NAMD, GROMACS (for alternative MD engines)

## Installation Methods

### Using pip (recommended)

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

# Clone and install pyMDMix
git clone https://github.com/DAlvGar/pyMDmix.git
cd pyMDmix
pip install -e ".[dev]"
```

## Verify Installation

```bash
# Check CLI is available
python -m pymdmix.cli --help

# Or if installed with entry point
pymdmix --help

# Check version
pymdmix --version

# List available solvents
pymdmix info solvents

# Run tests
pytest tests/ -q
```

## Environment Variables

pyMDMix uses the following environment variables:

| Variable | Purpose | Default |
|----------|---------|---------|
| `AMBERHOME` | Path to Amber/AmberTools installation | Auto-detected |
| `MDMIX_DATA` | Custom solvent library path | Package default |

## Dependencies

Core dependencies (installed automatically):

- `numpy` - Numerical operations
- `parmed` - Molecular structure handling
- `mdanalysis` - Trajectory analysis
- `click` - CLI framework
- `pyyaml` - Configuration files

Optional dependencies:

- `openmm` - OpenMM MD engine support
- `matplotlib` - Plotting capabilities
- `scipy` - Advanced analysis functions

## Troubleshooting

### AmberTools not found

Ensure `$AMBERHOME` is set correctly:

```bash
export AMBERHOME=/path/to/amber
source $AMBERHOME/amber.sh
```

### Import errors

If you encounter import errors, ensure all dependencies are installed:

```bash
pip install -e ".[dev]" --force-reinstall
```

### Permission issues with solvent database

If you can't modify the system solvent database, pyMDMix will automatically create a user-local copy in `~/.mdmix/SOLVENTS.db`.
