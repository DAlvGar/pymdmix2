# Trajectory Alignment

Before analysis, trajectories must be aligned to a common reference structure. This removes translational and rotational motion, ensuring probe density is calculated in the protein's frame of reference.

## What Alignment Does

1. **Center and image** - Fix periodic boundary artifacts
2. **RMSD alignment** - Superimpose on reference structure
3. **Save aligned trajectory** - For density calculation
4. **Calculate RMSD** - Quality metrics
5. **Generate averages** - Per-step average structures

---

## Running Alignment

### Basic Usage

```bash
# Align all replicas
pymdmix analyze align all

# Specific solvents
pymdmix analyze align bysolvent -s ETA MAM

# Parallel execution
pymdmix analyze align all -C 4
```

### Options

```bash
# Select nanoseconds
pymdmix analyze align all -N 1:10

# Custom alignment mask
pymdmix analyze align all --mask 1-100

# Custom reference structure
pymdmix analyze align all --ref custom_ref.pdb
```

---

## Alignment Mask

By default, alignment uses backbone atoms of the protein. Customize with:

```bash
# Residues 10-150
pymdmix analyze align all --mask 10-150

# Specific atoms
pymdmix analyze align all --mask "1-100@CA,C,N,O"
```

---

## Output Files

```
align/
├── md1.nc              # Aligned trajectory
├── md2.nc
├── md1_bb_rmsd.dat     # Backbone RMSD vs reference
├── md1_ha_rmsd.dat     # Heavy atom RMSD
├── prot_avg_1.pdb      # Average structure (step 1)
└── prot_avg_2.pdb
```

---

## Checking Alignment Quality

### Plot RMSD

```bash
pymdmix plot rmsd all
pymdmix plot rmsd byname -s MyProtein_ETA_1
```

Good alignment shows:
- RMSD plateaus (not drifting)
- No sudden jumps
- Reasonable values (typically 1-3 Å for stable proteins)

### Visual Inspection

Load aligned trajectory in VMD/PyMOL and verify:
- Protein stays centered
- No imaging artifacts (atoms jumping across box)
- Smooth motion

---

## Multi-chain Proteins

For multi-chain proteins, imaging can be complex. pyMDMix generates per-chain imaging commands, but manual adjustment may be needed.

### Edit Alignment Scripts

```bash
# View generated script
cat MyProtein_ETA_1/align/md1.ptraj

# Edit if needed
vim MyProtein_ETA_1/align/md1.ptraj

# Re-run manually
cpptraj -i md1.ptraj
```

### Common Imaging Fixes

```
# Center on first chain
center :1-200 mass origin

# Image each chain separately
image :1-200 origin center familiar
image :201-400 origin center familiar
```

---

## Python API

```python
from pymdmix.analysis import AlignAction
from pymdmix.project import Project

project = Project.load("myproject")
replica = project.get_replica("MyProtein_ETA_1")

action = AlignAction(
    mask="1-150",
    reference="custom_ref.pdb"
)
results = action.run(replica, nprocs=4)

# Check RMSD
print(f"Mean RMSD: {results.mean_rmsd:.2f} Å")
print(f"Max RMSD: {results.max_rmsd:.2f} Å")
```

---

## Troubleshooting

### High RMSD

- Check simulation stability
- Verify correct reference structure
- Consider different alignment mask

### Imaging Artifacts

- Edit ptraj scripts manually
- Check chain definitions
- Try different imaging keywords

### Memory Issues

```bash
# Process fewer steps at once
pymdmix analyze align all -N 1:5
pymdmix analyze align all -N 6:10
```

---

## Next Steps

- [Density Calculation](density.md) - Calculate probe grids
