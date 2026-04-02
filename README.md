# pyMDMix2

A modern Python toolkit for molecular dynamics simulations with organic solvent mixtures, used in drug discovery for identifying binding hotspots on protein surfaces.

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![Tests](https://img.shields.io/badge/tests-698%20passed-brightgreen.svg)]()

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

## Cloud Execution on AWS EC2 GPU Instances

pyMDMix2 can launch MD replicas directly on AWS EC2 GPU instances (`g4dn.xlarge` and
similar), upload input data to S3, and automatically download results when the simulation
finishes — all from the CLI.

> **Note:** All existing local and HPC workflows are completely unaffected.
> The cloud feature is opt-in and requires the `cloud` extras.

### Prerequisites

1. **AWS account** with permissions for EC2, S3, and STS.
2. **An S3 bucket** for staging data (e.g. `my-pymdmix-staging`).
3. **An EC2 key pair** in the target region (download the `.pem` private key).
4. **A security group** that allows inbound SSH (TCP 22) from your IP.
5. Optional but recommended: an IAM instance profile with S3 read/write permissions,
   so the instance never needs static credentials.

### Installation

```bash
# Install the cloud extras (boto3 + paramiko)
pip install "pymdmix[cloud]"
# or with uv:
uv sync --extra cloud
```

### AWS Credentials

pyMDMix2 **never stores AWS credentials** in the project config. It relies on the
standard boto3 credential chain:

| Method | How |
|--------|-----|
| Environment variables | `export AWS_ACCESS_KEY_ID=… AWS_SECRET_ACCESS_KEY=…` |
| Shared credentials file | `~/.aws/credentials` (created by `aws configure`) |
| IAM instance profile | Automatic when running on EC2 |

Store the path to your SSH private key in an environment variable:

```bash
export PYMDMIX_AWS_KEY_PATH=/path/to/my-keypair.pem
```

### Verify credentials

```bash
pymdmix cloud test-connection
# ✓ AWS credentials are valid
#   Account:  123456789012
#   User ARN: arn:aws:iam::123456789012:user/me
```

### Configure a project for cloud execution

Run the interactive wizard once per project:

```bash
cd my_project
pymdmix cloud configure
```

Or pass all options as flags (suitable for CI/automation):

```bash
pymdmix cloud configure \
  --region      us-east-1 \
  --instance-type g4dn.xlarge \
  --key-pair    my-keypair \
  --s3-bucket   my-pymdmix-staging \
  --s3-prefix   pymdmix/ \
  --ami-id      ami-0abc123def456     # optional: skip bootstrap install
```

Settings are saved to the project config file. To inspect available AMIs and
instance type pricing:

```bash
pymdmix cloud list-amis
pymdmix cloud estimated-cost          # rough cost estimate for all pending replicas
pymdmix cloud estimated-cost --hours 6
```

#### AMI Strategy

| Option | How | Startup time |
|--------|-----|-------------|
| **Pre-baked AMI** (recommended) | Set `--ami-id` to an Ubuntu 22.04 + CUDA + AmberTools 23 AMI | ~3 min |
| **Bootstrap AMI** | Leave `ami-id` blank — AmberTools is installed on first boot | ~13 min |

AMI IDs are maintained in `pymdmix/data/cloud/amis.json` (publication in progress).

### Launching replicas

```bash
# Launch a single replica and wait (optionally stream the log)
pymdmix run replica ETA_1
pymdmix run replica ETA_1 --watch          # stream bootstrap.log via SSH

# Launch all pending replicas in parallel
pymdmix run all
pymdmix run all --max-parallel 4           # limit to 4 concurrent instances

# Dry run — show what would happen without spending money
pymdmix run all --dry-run
```

What happens under the hood:
1. Input files (`*.prmtop`, `*.inpcrd`, `*.in`, `COMMANDS.sh`) are uploaded to S3.
2. An EC2 Spot instance is launched with a cloud-init bootstrap script.
3. The instance downloads the inputs, runs `bash COMMANDS.sh`, uploads results, then writes a `DONE` sentinel to S3.
4. Optionally, the instance self-terminates.

### Monitoring progress

```bash
# Show cached status (fast, no AWS calls)
pymdmix run status

# Poll AWS + S3 for live updates and auto-download completed results
pymdmix run status --poll
```

Example output:

```
REPLICA                        STATE        INSTANCE             IP               ELAPSED
------------------------------------------------------------------------
ETA_1                          running      i-0a1b2c3d4e5f678   1.2.3.4          00:47:12
ETA_2                          done         i-0b2c3d4e5f6789a   —                01:03:45
MAM_1                          launching    i-0c3d4e5f678901b   —                00:01:02
```

### Streaming logs

```bash
# Tail the bootstrap log via SSH (requires PYMDMIX_AWS_KEY_PATH)
pymdmix run logs ETA_1

# Fetch the last-synced log from S3 (no SSH key needed)
pymdmix run logs ETA_1 --s3
```

### Retrieving results

Results are automatically downloaded when `pymdmix run status --poll` detects a
`DONE` sentinel. To trigger a manual download:

```bash
pymdmix run fetch ETA_1
```

### Cancelling a job

```bash
pymdmix run cancel ETA_1            # terminates the EC2 instance (default)
pymdmix run cancel ETA_1 --no-terminate  # stops without terminating (cheaper resume)
```

### Full cloud workflow example

```bash
# 0. One-time setup
export AWS_ACCESS_KEY_ID=...
export AWS_SECRET_ACCESS_KEY=...
export PYMDMIX_AWS_KEY_PATH=~/keys/my-keypair.pem

pip install "pymdmix[cloud]"
pymdmix cloud test-connection

# 1. Create project and replicas (same as local workflow)
pymdmix create project -n my_project
cd my_project
pymdmix add system -f system.cfg
pymdmix add replica -f replica.cfg --count 3

# 2. Configure AWS settings
pymdmix cloud configure \
  --region us-east-1 \
  --key-pair my-keypair \
  --s3-bucket my-staging-bucket \
  --instance-type g4dn.xlarge

# 3. Launch all replicas on EC2
pymdmix run all

# 4. Monitor and auto-fetch results
watch -n 60 "pymdmix run status --poll"

# 5. Analyse locally once all replicas are done
pymdmix analyze align all
pymdmix analyze density all
pymdmix analyze hotspots all --threshold -0.5
```

### Security notes

- **AWS credentials** are never stored in project files — use env vars or `~/.aws/credentials`.
- **SSH private key** path is stored only in `PYMDMIX_AWS_KEY_PATH` (env var); never committed.
- SSH host key validation is enforced (`RejectPolicy`). Add the instance's key with
  `ssh-keyscan <ip> >> ~/.ssh/known_hosts` before using `--watch` or `pymdmix run logs`.
- Use `iam_instance_profile` (set during `cloud configure`) to give the EC2 instance
  S3 access without any static AWS credentials on the instance.
- S3 bucket server-side encryption (SSE-S3 or SSE-KMS) is strongly recommended.
- Set your security group to allow SSH only from your current public IP.

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
├── cloud/          # AWS EC2/S3 cloud execution (boto3, paramiko)
├── project/        # Replica, settings, config
├── engines/        # Amber, OpenMM, NAMD, GROMACS
├── cli.py          # Command-line interface
└── data/
    ├── solvents/   # Solvent definitions
    └── cloud/      # AMI catalog and instance type pricing
```

## Documentation

See `docs/` and `examples/` for detailed guides and configuration templates.

## License

MIT License - see [LICENSE](LICENSE)

## Citation

> Alvarez-Garcia, D., & Barril, X. (2014). Molecular simulations with solvent competition quantify water displaceability and provide accurate interaction maps of protein binding sites. *J. Med. Chem.*, 57(20), 8530-8539.

## Author

**Daniel Alvarez-Garcia** - [algarcia.daniel@gmail.com](mailto:algarcia.daniel@gmail.com)
