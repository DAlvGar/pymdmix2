# Energy Conversion

After calculating density grids, they are converted to **free energy grids** (ΔG) using the Boltzmann relationship. These energy values represent the binding affinity of each probe type at each grid point.

## Theory

The free energy is calculated from density using:

```
ΔG = -RT × ln(ρ)
```

Where:
- `R` = Gas constant (1.987 cal/mol·K)
- `T` = Temperature (300 K by default)
- `ρ` = Density (relative to bulk)

| Density | ΔG (kcal/mol) | Interpretation |
|---------|---------------|----------------|
| 2.0 | -0.41 | Favorable |
| 1.0 | 0.0 | Neutral (bulk) |
| 0.5 | +0.41 | Unfavorable |
| 0.1 | +1.37 | Strongly unfavorable |

More negative ΔG = more favorable binding.

---

## Running Energy Conversion

### Basic Usage

```bash
# Convert all density grids to energy
pymdmix analyze energy all

# For specific solvents
pymdmix analyze energy bysolvent -s ETA

# For specific replicas
pymdmix analyze energy byname -s MyProtein_ETA_1 MyProtein_ETA_2
```

### With Options

```bash
# Custom temperature
pymdmix analyze energy all --temp 310

# Set minimum density threshold (avoid log(0))
pymdmix analyze energy all --min-density 0.001
```

---

## Output Files

Energy grids are created alongside density grids:

```
grids/
├── ETA_CT_density.dx    # Input density
├── ETA_CT_DG.dx         # Output free energy
├── ETA_OH_density.dx
├── ETA_OH_DG.dx
└── ...
```

---

## Averaging Across Replicas

For robust results, average energy grids from multiple replicas:

```bash
# Average energy grids from all ETA replicas
pymdmix analyze energy bysolvent -s ETA --average
```

This creates:
- `ETA_CT_DG_avg.dx` - Averaged across replicas
- `ETA_CT_DG_std.dx` - Standard deviation (optional)

---

## Energy Interpretation

### Favorable Binding

| ΔG (kcal/mol) | Strength |
|---------------|----------|
| < -1.5 | Very strong hotspot |
| -1.5 to -1.0 | Strong hotspot |
| -1.0 to -0.5 | Moderate hotspot |
| -0.5 to 0.0 | Weak preference |

### Typical Thresholds

- **Visualization**: -0.5 to -1.0 kcal/mol isosurface
- **Hotspot detection**: -1.0 to -1.5 kcal/mol
- **Druggability**: Sites with ΔG < -1.5 kcal/mol

---

## Python API

```python
from pymdmix.analysis import EnergyAction
from pymdmix.core import Grid

# Run energy conversion
project = Project.load("myproject")
replica = project.get_replica("MyProtein_ETA_1")

action = EnergyAction(temperature=300)
results = action.run(replica)

# Access energy grids
for probe_name, grid in results.grids.items():
    print(f"{probe_name}: min ΔG = {grid.data.min():.2f} kcal/mol")

# Manual conversion from density
density_grid = Grid.from_dx("grids/ETA_CT_density.dx")
energy_grid = density_grid.to_energy(temperature=300)
energy_grid.write_dx("grids/ETA_CT_DG.dx")

# Average multiple grids
from pymdmix.core import Grid
grids = [
    Grid.from_dx(f"MyProtein_ETA_{i}/grids/ETA_CT_DG.dx")
    for i in range(1, 4)
]
avg_grid = Grid.average(grids)
avg_grid.write_dx("ETA_CT_DG_avg.dx")
```

---

## Combining Probe Types

Combine grids from the same chemical type across different solvents:

```python
from pymdmix.core import Grid

# Load hydrophobic probes from different solvents
eta_ct = Grid.from_dx("MyProtein_ETA_1/grids/ETA_CT_DG.dx")
iso_ct = Grid.from_dx("MyProtein_ISO_1/grids/ISO_CT_DG.dx")

# Take minimum (most favorable) at each point
combined = Grid.minimum([eta_ct, iso_ct])
combined.write_dx("combined_hydrophobic_DG.dx")
```

---

## Advanced Options

### Smoothing

Apply Gaussian smoothing to reduce noise:

```bash
pymdmix analyze energy all --smooth 1.0  # Sigma in Angstroms
```

### Normalization Methods

Different normalization approaches:

```python
action = EnergyAction(
    normalization="bulk"  # "bulk", "local", "raw"
)
```

### Custom Reference

Use a custom reference state:

```python
action = EnergyAction(
    reference_density=0.033  # Custom bulk density (molecules/Å³)
)
```

---

## Visualization Tips

### PyMOL

```python
# Load energy grid
load grids/ETA_CT_DG.dx, ct_energy

# Favorable regions (negative ΔG)
isosurface favorable, ct_energy, -1.0
color green, favorable

# Unfavorable regions (positive ΔG)
isosurface unfavorable, ct_energy, 1.0
color red, unfavorable
```

### VMD

```tcl
mol addfile grids/ETA_CT_DG.dx type dx
mol modstyle 0 1 Isosurface -1.0 0 0 0 1 1
mol modcolor 0 1 ColorID 7  # Green
```

---

## Troubleshooting

### All Energies Near Zero

- Check density calculation succeeded
- Verify there's sufficient probe density
- Ensure trajectory was properly aligned

### Very Negative Energies

- May indicate clipping artifacts
- Check for extremely high density spikes
- Consider smoothing

### Inconsistent Results Across Replicas

- Normal variation expected
- Use replica averaging
- Check for convergence (compare 10ns vs 20ns)

---

## Next Steps

- [Hotspot Detection](hotspots.md) - Identify binding sites
- [Residence Analysis](residence.md) - Solvent kinetics
