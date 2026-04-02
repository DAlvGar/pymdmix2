# Quick Start

This guide walks you through setting up your first pyMDMix project using the
primary single-file workflow.

## Overview

The typical pyMDMix workflow:

1. **Prepare** your protein structure and create an Amber Object File (OFF)
2. **Get a template** config file and edit it
3. **Create** the project — system and all replicas in one command
4. **Run** MD simulations (on your cluster)
5. **Analyze** trajectories to identify binding hotspots

---

## Step 1: Prepare Your Structure

You need either:
- An **Amber Object File (OFF)** with your parameterised system (recommended), or
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

---

## Step 2: Get a Project Template

```bash
# Save an annotated project template to the current directory
pymdmix create template
# → writes project.cfg
```

---

## Step 3: Edit the Template

Open `project.cfg` and fill in your system and simulation parameters:

```ini
[SYSTEM]
NAME = MyProtein
OFF  = /path/to/myprotein.off

[MDSETTINGS]
SOLVENTS = ETA, MAM, WAT
NREPL    = 3        # 3 replicas per solvent
NANOS    = 20       # 20 ns production
TEMP     = 300      # 300 K
RESTR    = FREE     # no positional restraints
```

All options have sensible defaults — only `NAME`, the structure path, and
`SOLVENTS` are strictly required.

---

## Step 4: Create the Project

```bash
pymdmix create project -n myproject -f project.cfg
```

This single command:
- Creates the project directory structure
- Registers the system
- Creates all replicas (here: 3 solvents × 3 replicas = **9 replicas**)

Check what was created:

```bash
pymdmix info project
```

---

## Step 5: Run Simulations

Each replica folder contains a `COMMANDS.sh` script with all simulation steps:

```bash
cd myproject/MyProtein_ETA_1/
cat COMMANDS.sh            # review generated commands
./COMMANDS.sh              # run locally, or:
# sbatch queue.sh          # submit to SLURM
```

---

## Step 6: Analyze Results

After simulations complete, run analysis from the project root:

```bash
# Align all trajectories
pymdmix analyze align all

# Calculate density grids
pymdmix analyze density all

# Convert density to free energy maps
pymdmix analyze energy all

# Find binding hotspots
pymdmix analyze hotspots all
```

---

## Output Files

After analysis, key outputs include:

| File | Description |
|------|-------------|
| `*_DG.dx` | Energy grid (free energy in kcal/mol) |
| `*_density.dx` | Raw probe density grid |
| `hotspots.pdb` | Detected binding hotspots |
| `*_rmsd.dat` | RMSD over time |

---

## Alternative: Step-by-step Setup

If you prefer to add systems and replicas separately:

```bash
pymdmix create project -n myproject
pymdmix add system  -f system.cfg
pymdmix add replica -f settings.cfg
```

See [User Guide: Projects](user-guide/projects.md) for the full reference.

---

## Next Steps

- [User Guide: Projects](user-guide/projects.md) — project management
- [User Guide: MD Settings](user-guide/md-settings.md) — simulation parameters
- [Tutorials: Toy Project](tutorials/toy-project-setup.md) — complete walkthrough
