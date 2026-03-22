# Quick Start

This guide walks you through setting up your first pyMDMix project.

## Overview

The typical pyMDMix workflow involves:

1. **Prepare** your protein structure (PDB or Amber Object File)
2. **Create** a project and add your system
3. **Generate** replicas with chosen solvents and settings
4. **Run** MD simulations (on your cluster)
5. **Analyze** trajectories to identify binding hotspots

## Step 1: Prepare Your Structure

You need either:
- An **Amber Object File (OFF)** with your parameterized system, or
- A **clean PDB file** with correct protonation states

For a simple protein:

```bash
# In tLeap
source leaprc.protein.ff14SB
system = loadPdb myprotein.pdb
check system
saveOff system myprotein.off
quit
```

## Step 2: Create a Project

```bash
# Create a new project directory
pymdmix create project myproject
cd myproject
```

## Step 3: Add Your System

Create a system configuration file (`system.cfg`):

```ini
[SYSTEM]
NAME = MyProtein
OFF = /path/to/myprotein.off
```

Add it to the project:

```bash
pymdmix add system -f system.cfg
```

## Step 4: Create Replicas

Create a replica configuration file (`replica.cfg`):

```ini
[REPLICA]
SYSTEM = MyProtein
SOLVENT = ETA
NANOS = 20
RESTRMODE = FREE
```

Generate the replica:

```bash
pymdmix add replica -f replica.cfg
```

This creates:
- Input files for minimization, equilibration, and production
- Queue submission scripts
- Reference structures

## Step 5: Run Simulations

The replica folder contains a `COMMANDS.sh` script with all commands:

```bash
cd MyProtein_ETA_1/
cat COMMANDS.sh

# Submit to your cluster or run locally
# Modify queue scripts as needed for your HPC system
```

## Step 6: Analyze Results

After simulation completes, bring trajectories back and analyze:

```bash
# Align all trajectories
pymdmix analyze align all

# Calculate density grids
pymdmix analyze density all

# Convert to energy maps
pymdmix analyze energy all

# Find hotspots
pymdmix analyze hotspots all
```

## View Results

```bash
# Check project status
pymdmix info project

# List replicas
pymdmix info replicas

# Plot RMSD
pymdmix plot rmsd all
```

## Output Files

After analysis, key outputs include:

| File | Description |
|------|-------------|
| `*_DG.dx` | Energy grid (free energy in kcal/mol) |
| `*_density.dx` | Raw density grid |
| `hotspots.pdb` | Detected binding hotspots |
| `*_rmsd.dat` | RMSD trajectories |

## Next Steps

- [System Preparation](user-guide/system-preparation.md) - Detailed structure prep
- [Solvent Mixtures](user-guide/solvents.md) - Available solvents
- [Analysis Guide](analysis/overview.md) - Full analysis workflow
- [Tutorials](tutorials/toy-project-setup.md) - Complete walkthrough
