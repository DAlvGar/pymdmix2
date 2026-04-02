# Tutorial: Setting Up a Toy Project

This tutorial walks through setting up a complete pyMDMix project using a
small protein, demonstrating the primary single-file workflow.

## Overview

We'll:
1. Download a test protein (BPTI - 58 residues)
2. Prepare it for simulation
3. Create a `project.cfg` file
4. Bootstrap the entire project from that single file

---

## Prerequisites

- pyMDMix2 installed
- AmberTools available (`tleap`, `cpptraj`)
- ~30 minutes

---

## Step 1: Get the Protein

Download BPTI (Bovine Pancreatic Trypsin Inhibitor):

```bash
mkdir mdmix_tutorial && cd mdmix_tutorial
curl -o 5pti.pdb "https://files.rcsb.org/download/5PTI.pdb"
```

---

## Step 2: Prepare the Structure

### Clean the PDB

```bash
grep "^ATOM" 5pti.pdb > bpti_clean.pdb
```

### Create Amber Object File

```bash
tleap
```

In tLeap:

```
source leaprc.protein.ff14SB

protein = loadPdb bpti_clean.pdb
check protein

# BPTI has 3 disulfide bonds
bond protein.5.SG  protein.55.SG
bond protein.14.SG protein.38.SG
bond protein.30.SG protein.51.SG

check protein
saveOff protein bpti.off
saveAmberParm protein bpti_test.prmtop bpti_test.inpcrd

quit
```

---

## Step 3: Get a Project Template

```bash
pymdmix create template
```

This writes `project.cfg` to the current directory.

---

## Step 4: Edit the Template

Open `project.cfg` and set:

```ini
[SYSTEM]
NAME = BPTI
OFF  = bpti.off

[MDSETTINGS]
SOLVENTS = ETA, MAM
NREPL    = 3
NANOS    = 5
TEMP     = 300
RESTR    = FREE
```

---

## Step 5: Create the Project

```bash
pymdmix create project -n bpti_mdmix -f project.cfg
```

Expected output:

```
Creating project 'bpti_mdmix' in .../bpti_mdmix
✓ Project 'bpti_mdmix' created
  System:   BPTI
  Solvents: ETA, MAM
  Replicas: 6
    - BPTI_ETA_1
    - BPTI_ETA_2
    - BPTI_ETA_3
    - BPTI_MAM_1
    - BPTI_MAM_2
    - BPTI_MAM_3
```

---

## Step 6: Check Project Status

```bash
cd bpti_mdmix
pymdmix info project
```

---

## Step 7: Examine Replica Contents

```bash
ls BPTI_ETA_1/
```

```
BPTI_ETA_1/
├── BPTI_ETA_1.prmtop
├── BPTI_ETA_1.inpcrd
├── BPTI_ETA_1_ref.pdb
├── COMMANDS.sh
├── queue.sh
├── min/
├── eq/
└── md/
```

---

## Step 8: Review Commands

```bash
cat BPTI_ETA_1/COMMANDS.sh
```

The script shows:
1. Minimization (min1, min2)
2. Equilibration (eq1, eq2, ...)
3. Production (md1, md2, md3, md4, md5)

---

## Project Structure

```
bpti_mdmix/
├── .mdmix/
│   └── project.json
├── inputs/
│   └── project.cfg
├── BPTI_ETA_1/
├── BPTI_ETA_2/
├── BPTI_ETA_3/
├── BPTI_MAM_1/
├── BPTI_MAM_2/
└── BPTI_MAM_3/
```

---

## Quick Test (Optional)

Run a very short test locally:

```bash
cd BPTI_ETA_1/min
pmemd -O -i min1.in -o min1.out -p ../BPTI_ETA_1.prmtop \
      -c ../BPTI_ETA_1.inpcrd -r min1.rst -ref ../BPTI_ETA_1.inpcrd
tail min1.out
```

---

## Next Steps

You're ready to run simulations! See:

- [Tutorial: Running Simulations](toy-project-simulation.md)
- [Queue Scripts](../user-guide/queue-scripts.md) — HPC submission
