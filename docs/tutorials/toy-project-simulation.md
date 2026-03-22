# Tutorial: Running and Analyzing Simulations

Continuing from the [setup tutorial](toy-project-setup.md), we'll run the simulations and analyze results.

## Overview

We'll:
1. Run MD simulations (or use pre-computed data)
2. Align trajectories
3. Calculate density grids
4. Convert to energy maps
5. Detect hotspots
6. Visualize results

---

## Step 1: Run Simulations

### Option A: Local Testing (Short)

For testing, modify inputs for a very short run:

```bash
cd BPTI_ETA_1

# Edit production to run only 1000 steps
sed -i 's/nstlim=500000/nstlim=1000/' md/md1.in

# Run (will take ~5 minutes with GPU)
bash COMMANDS.sh
```

### Option B: Full Simulation (HPC)

Edit queue script for your cluster and submit:

```bash
cd BPTI_ETA_1
vim queue.sh  # Adjust for your system
sbatch queue.sh
```

### Option C: Use Pre-computed Data

For this tutorial, you can download example trajectories:

```bash
# Download example data (if available)
curl -O https://example.com/bpti_tutorial_data.tar.gz
tar xzf bpti_tutorial_data.tar.gz
```

---

## Step 2: Verify Simulation Output

After simulation completes:

```bash
ls BPTI_ETA_1/md/
```

Expected files:
```
md1.nc  md1.out  md1.rst
md2.nc  md2.out  md2.rst
...
```

Check for errors:

```bash
grep -i error BPTI_ETA_1/md/*.out
tail BPTI_ETA_1/md/md5.out  # Should show "Total wall time"
```

---

## Step 3: Align Trajectories

```bash
cd bpti_mdmix

# Align all replicas
pymdmix analyze align all

# Or with parallel execution
pymdmix analyze align all -C 4
```

Check alignment quality:

```bash
pymdmix plot rmsd all
```

This creates RMSD plots in each replica's `align/` folder.

---

## Step 4: Calculate Density

```bash
# Calculate probe density for all replicas
pymdmix analyze density all -C 4
```

This creates density grids in `grids/`:

```bash
ls BPTI_ETA_1/grids/
# ETA_CT_density.dx  ETA_OH_density.dx  ETA_WAT_density.dx
```

---

## Step 5: Convert to Energy

```bash
pymdmix analyze energy all
```

Creates energy grids:

```bash
ls BPTI_ETA_1/grids/
# ETA_CT_DG.dx  ETA_OH_DG.dx  ETA_WAT_DG.dx
```

---

## Step 6: Detect Hotspots

```bash
pymdmix analyze hotspots all
```

Results in `hotspots/`:

```bash
ls BPTI_ETA_1/hotspots/
# hotspots.pdb  hotspots.csv  summary.txt

cat BPTI_ETA_1/hotspots/summary.txt
```

---

## Step 7: Visualize Results

### PyMOL Visualization

```bash
pymol BPTI_ETA_1/BPTI_ETA_1_ref.pdb
```

In PyMOL:

```python
# Load energy grids
load BPTI_ETA_1/grids/ETA_CT_DG.dx, hydrophobic
load BPTI_ETA_1/grids/ETA_OH_DG.dx, hbond

# Show protein surface
show surface, BPTI_ETA_1_ref
set transparency, 0.7

# Show favorable hydrophobic regions
isosurface ct_surf, hydrophobic, -1.0
color green, ct_surf

# Show favorable H-bond regions
isosurface oh_surf, hbond, -1.0
color cyan, oh_surf

# Load hotspots
load BPTI_ETA_1/hotspots/hotspots.pdb, hotspots
show spheres, hotspots
spectrum b, red_white_green, hotspots
```

### Save PyMOL Session

```python
save bpti_analysis.pse
```

---

## Step 8: Compare Solvents

Compare hotspots from different solvents:

```bash
# Combine hotspots
pymdmix tools combine-hotspots \
    BPTI_ETA_1/hotspots/hotspots.pdb \
    BPTI_MAM_1/hotspots/hotspots.pdb \
    -o combined_hotspots.pdb
```

---

## Step 9: Generate Report

```bash
pymdmix info analysis BPTI_ETA_1
```

Output:

```
Replica: BPTI_ETA_1
Solvent: ETA (Ethanol 20%)
Simulation: 5 ns
Status: Complete

Alignment:
  Mean RMSD: 1.23 Å
  Max RMSD: 2.01 Å

Density Grids: 3
  - ETA_CT: max density 4.52
  - ETA_OH: max density 3.21
  - ETA_WAT: max density 2.87

Energy Grids: 3
  - ETA_CT: min ΔG = -1.82 kcal/mol
  - ETA_OH: min ΔG = -1.45 kcal/mol

Hotspots: 5 detected
  1. (-2.31 kcal/mol) at binding pocket
  2. (-1.87 kcal/mol) at surface groove
  ...
```

---

## Results Interpretation

### BPTI Binding Sites

BPTI is a trypsin inhibitor with a well-characterized binding loop. You should see:

1. **Primary hotspot** at the reactive site loop (residues 14-18)
2. **Secondary hotspots** at surface grooves
3. **Hydrophobic hotspots** (green) at core-exposed regions
4. **Polar hotspots** (cyan) near charged residues

### Comparing Probes

| Probe | Chemical Type | Expected Location |
|-------|---------------|-------------------|
| ETA_CT | Hydrophobic | Core, loops |
| ETA_OH | H-bond D/A | Near Asp, Glu, Lys |
| MAM_N | H-bond donor | Backbone interactions |
| MAM_O | H-bond acceptor | Near NH groups |

---

## Summary

You've completed a full MDMix analysis:

✅ Created project and replicas
✅ Ran MD simulations
✅ Aligned trajectories
✅ Calculated density grids
✅ Converted to energy maps
✅ Detected binding hotspots
✅ Visualized results

---

## Next Steps

- Try different solvents (ANT, ISO)
- Run longer simulations (20+ ns)
- Apply to your protein of interest
- See [Analysis Guide](../analysis/overview.md) for advanced options
