# Tutorial: Setting Up a Toy Project

This tutorial walks through setting up a complete pyMDMix project using a small protein.

## Overview

We'll:
1. Download a test protein (BPTI - 58 residues)
2. Prepare it for simulation
3. Create a pyMDMix project
4. Set up replicas with different solvents

---

## Prerequisites

- pyMDMix installed
- AmberTools available (`tleap`, `cpptraj`)
- ~30 minutes

---

## Step 1: Get the Protein

Download BPTI (Bovine Pancreatic Trypsin Inhibitor):

```bash
# Create working directory
mkdir mdmix_tutorial && cd mdmix_tutorial

# Download from PDB
curl -o 5pti.pdb "https://files.rcsb.org/download/5PTI.pdb"
```

---

## Step 2: Prepare the Structure

### Clean the PDB

```bash
# Keep only protein atoms, remove waters and ligands
grep "^ATOM" 5pti.pdb > bpti_clean.pdb
```

### Create Amber Object File

```bash
# Start tleap
tleap
```

In tleap:

```
# Load force field
source leaprc.protein.ff14SB

# Load structure
protein = loadPdb bpti_clean.pdb

# Check for problems
check protein

# Add disulfide bonds (BPTI has 3)
bond protein.5.SG protein.55.SG
bond protein.14.SG protein.38.SG
bond protein.30.SG protein.51.SG

# Check again
check protein

# Save object file
saveOff protein bpti.off

# Verify it works
saveAmberParm protein bpti_test.prmtop bpti_test.inpcrd

# Exit
quit
```

---

## Step 3: Create the Project

```bash
# Create project
pymdmix create project bpti_mdmix
cd bpti_mdmix

# Copy the OFF file
cp ../bpti.off .
```

---

## Step 4: Add the System

Create `system.cfg`:

```ini
[SYSTEM]
NAME = BPTI
OFF = bpti.off
```

Add to project:

```bash
pymdmix add system -f system.cfg
```

Verify:

```bash
pymdmix info systems
```

---

## Step 5: Create Replicas

### Ethanol Replica

Create `replica_eta.cfg`:

```ini
[REPLICA]
SYSTEM = BPTI
SOLVENT = ETA
NANOS = 5
RESTRMODE = FREE
```

```bash
pymdmix add replica -f replica_eta.cfg
```

### Acetamide Replica

Create `replica_mam.cfg`:

```ini
[REPLICA]
SYSTEM = BPTI
SOLVENT = MAM
NANOS = 5
RESTRMODE = FREE
```

```bash
pymdmix add replica -f replica_mam.cfg
```

### Multiple Ethanol Replicas

```bash
# Create 2 more ETA replicas
pymdmix add replica -f replica_eta.cfg --count 2
```

---

## Step 6: Check Project Status

```bash
pymdmix info project
```

Expected output:

```
Project: bpti_mdmix
Systems: 1
  - BPTI

Replicas: 4
  - BPTI_ETA_1 (ETA, 5 ns, FREE)
  - BPTI_ETA_2 (ETA, 5 ns, FREE)
  - BPTI_ETA_3 (ETA, 5 ns, FREE)
  - BPTI_MAM_1 (MAM, 5 ns, FREE)
```

---

## Step 7: Examine Replica Contents

```bash
ls BPTI_ETA_1/
```

You should see:

```
BPTI_ETA_1/
├── BPTI_ETA_1.prmtop      # Solvated topology
├── BPTI_ETA_1.inpcrd      # Starting coordinates
├── BPTI_ETA_1_ref.pdb     # Reference structure
├── COMMANDS.sh            # Execution commands
├── queue.sh               # Queue submission
├── min/                   # Minimization inputs
├── eq/                    # Equilibration inputs
└── md/                    # Production inputs
```

---

## Step 8: Review Commands

```bash
cat BPTI_ETA_1/COMMANDS.sh
```

Shows the sequence of simulation steps:
1. Minimization (min1, min2)
2. Equilibration (eq1, eq2, ...)
3. Production (md1, md2, md3, md4, md5)

---

## Project Structure

```
bpti_mdmix/
├── .mdmix/
│   └── project.json
├── systems/
│   └── BPTI.json
├── bpti.off
├── system.cfg
├── replica_eta.cfg
├── replica_mam.cfg
├── BPTI_ETA_1/
├── BPTI_ETA_2/
├── BPTI_ETA_3/
└── BPTI_MAM_1/
```

---

## Next Steps

You're ready to run simulations! See:
- [Tutorial: Running Simulations](toy-project-simulation.md)
- [Queue Scripts](../user-guide/queue-scripts.md) - For HPC submission

---

## Quick Test (Optional)

Run a very short test locally:

```bash
cd BPTI_ETA_1

# Run just minimization
cd min
pmemd -O -i min1.in -o min1.out -p ../BPTI_ETA_1.prmtop \
      -c ../BPTI_ETA_1.inpcrd -r min1.rst -ref ../BPTI_ETA_1.inpcrd

# Check output
tail min1.out
```

If minimization completes without errors, your setup is correct!
