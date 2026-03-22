# Density Grid Calculation

Density calculation is the core of MDMix analysis. For each probe atom defined in the solvent mixture, pyMDMix calculates a 3D grid representing the spatial distribution of that probe around the protein.

## Theory

The **density grid** represents the normalized probability of finding a probe atom at each grid point:

```
ρ(x,y,z) = N_observed(x,y,z) / N_expected(x,y,z)
```

Where:
- `N_observed` = Number of probe visits to grid cell during simulation
- `N_expected` = Expected visits based on bulk concentration

A density > 1 indicates favorable binding; < 1 indicates unfavorable.

---

## Running Density Calculation

### Basic Usage

```bash
# Calculate density for all replicas
pymdmix analyze density all

# For specific solvents
pymdmix analyze density bysolvent -s ETA MAM

# For specific replicas
pymdmix analyze density byname -s MyProtein_ETA_1
```

### With Options

```bash
# Use 8 processors
pymdmix analyze density all -C 8

# Only analyze first 10 ns
pymdmix analyze density all -N 1:10

# Custom grid spacing (default: 0.5 Å)
pymdmix analyze density all --spacing 0.25
```

### Command Help

```bash
pymdmix analyze density --help
```

---

## Grid Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `spacing` | 0.5 Å | Grid resolution |
| `buffer` | 3.0 Å | Padding around protein |
| `mask` | auto | Atoms to exclude (protein) |

Higher resolution (smaller spacing) gives more detail but requires more memory and time.

---

## Output Files

For each probe in the solvent, a density grid is created:

```
grids/
├── ETA_CT_density.dx    # Carbon probe density
├── ETA_OH_density.dx    # Hydroxyl probe density
├── ETA_WAT_density.dx   # Water probe density
└── ...
```

---

## Understanding Density Values

| Density | Interpretation |
|---------|----------------|
| > 2.0 | Strong favorable binding |
| 1.5 - 2.0 | Moderate favorable binding |
| 1.0 | Bulk-like (neutral) |
| 0.5 - 1.0 | Slightly unfavorable |
| < 0.5 | Strongly unfavorable (excluded) |

---

## Combining Replicas

Average densities across replicas for better statistics:

```bash
# Average all ETA replicas
pymdmix analyze density bysolvent -s ETA --average

# This creates averaged grids in:
# grids/ETA_CT_density_avg.dx
```

---

## Python API

```python
from pymdmix.analysis import DensityAction
from pymdmix.project import Project
from pymdmix.core import Grid

# Run density calculation
project = Project.load("myproject")
replica = project.get_replica("MyProtein_ETA_1")

action = DensityAction(
    spacing=0.5,
    buffer=3.0,
    nprocs=4
)
results = action.run(replica)

# Access results
for probe_name, grid in results.grids.items():
    print(f"{probe_name}: max density = {grid.data.max():.2f}")
    grid.write_dx(f"{probe_name}_density.dx")

# Load and manipulate grids
grid = Grid.from_dx("grids/ETA_CT_density.dx")
print(f"Grid dimensions: {grid.shape}")
print(f"Grid spacing: {grid.spacing}")
print(f"Grid origin: {grid.origin}")
```

---

## Advanced Options

### Custom Probe Tracking

Track specific atoms not defined in solvent:

```python
from pymdmix.analysis import DensityAction

action = DensityAction()
action.add_probe("my_probe", residue="LIG", atoms=["C1", "C2"])
results = action.run(replica)
```

### Excluding Regions

Exclude specific protein regions from density calculation:

```bash
pymdmix analyze density all --exclude-mask ":100-120"
```

### Frame Selection

Analyze specific trajectory portions:

```python
action = DensityAction()
results = action.run(
    replica,
    start_frame=100,
    end_frame=1000,
    stride=2  # Every 2nd frame
)
```

---

## Performance Tips

1. **Use multiple CPUs** - Density calculation is parallelizable
2. **Reduce grid size** - Increase spacing for initial tests
3. **Select frames wisely** - Skip equilibration, use stride
4. **Check memory** - Large grids require significant RAM

Estimated memory usage:
```
Memory ≈ (box_size / spacing)³ × 8 bytes × n_probes
```

---

## Troubleshooting

### Out of Memory

```bash
# Reduce resolution
pymdmix analyze density all --spacing 1.0

# Process fewer frames
pymdmix analyze density all -N 10:20
```

### Slow Performance

```bash
# Use more CPUs
pymdmix analyze density all -C 16

# Skip frames
pymdmix analyze density all --stride 10
```

### Unexpected Density Values

- Check trajectory alignment was successful
- Verify probe definitions in solvent
- Ensure sufficient sampling (enough nanoseconds)

---

## Next Steps

- [Energy Conversion](energy.md) - Convert density to free energy
- [Hotspot Detection](hotspots.md) - Find binding sites
