# Analysis Overview

After running MD simulations, pyMDMix provides tools to extract binding information from the trajectories.

## Analysis Workflow

```
┌────────────────────────────────────────────────────────────────┐
│                    Analysis Pipeline                           │
├────────────────────────────────────────────────────────────────┤
│                                                                │
│  1. ALIGNMENT                                                  │
│     └── Align trajectories to reference structure              │
│                                                                │
│  2. DENSITY CALCULATION                                        │
│     └── Calculate probe density grids                          │
│                                                                │
│  3. ENERGY CONVERSION                                          │
│     └── Convert density to free energy (ΔG)                   │
│                                                                │
│  4. HOTSPOT DETECTION                                          │
│     └── Identify and rank binding hotspots                     │
│                                                                │
│  5. OPTIONAL: RESIDENCE ANALYSIS                               │
│     └── Study solvent exchange kinetics                        │
│                                                                │
└────────────────────────────────────────────────────────────────┘
```

---

## Prerequisites

Before analysis:

1. **Simulations must be complete**
2. **Trajectories copied** back to replica folders
3. **File names match** expected patterns (md1.nc, md2.nc, etc.)

Verify with:

```bash
ls MyProtein_ETA_1/md/*.nc
```

---

## Replica Selection Syntax

Most analysis commands use a selection syntax to specify which replicas to process:

### Selection Methods

| Method | Syntax | Description |
|--------|--------|-------------|
| `all` | `all` | All replicas in project |
| `bysolvent` | `bysolvent -s ETA MAM` | Replicas with specified solvents |
| `byname` | `byname -s MyProtein_ETA_1 MyProtein_ETA_2` | Specific replica names |
| `group` | `group -s mygroup` | Pre-defined replica group |

### Examples

```bash
# All replicas
pymdmix analyze align all

# Only ethanol replicas
pymdmix analyze density bysolvent -s ETA

# Specific replicas
pymdmix analyze energy byname -s MyProtein_ETA_1 MyProtein_ETA_2

# Named group
pymdmix analyze hotspots group -s ethanol_runs
```

---

## Common Options

### Step Selection

Select which nanoseconds (or trajectory files) to analyze:

```bash
# Analyze only ns 1-10
pymdmix analyze align all -N 1:10

# Analyze ns 5, 10, 15, 20
pymdmix analyze density all -N 5,10,15,20
```

### Parallel Execution

Use multiple CPUs for faster analysis:

```bash
# Use 8 processors
pymdmix analyze density all -C 8
```

---

## Quick Analysis Pipeline

Run the complete analysis:

```bash
cd myproject

# 1. Align all trajectories
pymdmix analyze align all -C 4

# 2. Check alignment (plot RMSD)
pymdmix plot rmsd all

# 3. Calculate density grids
pymdmix analyze density all -C 4

# 4. Convert to energy
pymdmix analyze energy all

# 5. Detect hotspots
pymdmix analyze hotspots all
```

---

## Output Files

### Per-Replica Outputs

Each replica folder will contain:

```
MyProtein_ETA_1/
├── align/
│   ├── md1.nc              # Aligned trajectory
│   ├── md1_bb_rmsd.dat     # Backbone RMSD
│   ├── md1_ha_rmsd.dat     # Heavy atom RMSD
│   └── prot_avg_1.pdb      # Average structure
├── grids/
│   ├── ETA_CT_density.dx   # Raw density
│   ├── ETA_CT_DG.dx        # Free energy (kcal/mol)
│   ├── ETA_OH_density.dx
│   ├── ETA_OH_DG.dx
│   └── ...
└── hotspots/
    └── hotspots.pdb        # Detected hotspots
```

### Grid File Format

Density and energy grids are saved in OpenDX format (`.dx`), compatible with:
- VMD
- PyMOL
- Chimera/ChimeraX

---

## Visualization

### In PyMOL

```python
# Load protein and grid
load MyProtein_ETA_1_ref.pdb
load grids/ETA_CT_DG.dx, ct_energy

# Show isosurface at -1 kcal/mol
isosurface ct_surf, ct_energy, -1.0
color green, ct_surf
```

### In VMD

```tcl
# Load structure
mol new MyProtein_ETA_1_ref.pdb

# Load energy grid
mol addfile grids/ETA_CT_DG.dx
mol modstyle 0 1 Isosurface -1.0 0 0 0 1 1
```

---

## Python API

```python
from pymdmix.analysis import DensityAction, HotspotAction
from pymdmix.project import Project

# Load project
project = Project.load("myproject")

# Run density analysis on one replica
replica = project.get_replica("MyProtein_ETA_1")
density_action = DensityAction()
results = density_action.run(replica, nprocs=4)

# Run hotspot detection
hotspot_action = HotspotAction()
hotspots = hotspot_action.run(replica)
for hs in hotspots:
    print(f"Hotspot at {hs.center}: {hs.energy:.2f} kcal/mol")
```

---

## Detailed Guides

- [Trajectory Alignment](alignment.md) - Detailed alignment procedures
- [Density Calculation](density.md) - Probe density computation
- [Energy Conversion](energy.md) - Density to free energy
- [Hotspot Detection](hotspots.md) - Identifying binding sites
- [Residence Analysis](residence.md) - Solvent exchange kinetics

---

## Troubleshooting

### Missing Trajectory Files

```
Error: Trajectory file md1.nc not found
```

Ensure simulation output is in the correct location with expected names.

### Memory Issues

For large systems, reduce memory usage:

```bash
# Process fewer frames
pymdmix analyze density all -N 1:10

# Use fewer parallel processes
pymdmix analyze density all -C 2
```

### Poor Alignment

Check RMSD plots for jumps or drift:

```bash
pymdmix plot rmsd all
```

If imaging is incorrect (multi-chain proteins), edit the ptraj scripts manually.
