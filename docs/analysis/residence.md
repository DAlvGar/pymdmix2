# Residence Analysis

Residence analysis studies **how long** probe molecules stay near protein sites, providing kinetic information complementary to the thermodynamic data from density grids.

## Theory

Residence time measures solvent exchange kinetics:
- **Long residence** → Tight binding, slow exchange
- **Short residence** → Weak binding, fast exchange

This helps distinguish:
- Stable binding sites (long residence + favorable energy)
- Transient contacts (short residence despite favorable energy)

---

## Running Residence Analysis

```bash
# Calculate residence times
pymdmix analyze residence all

# For specific replicas
pymdmix analyze residence byname -s MyProtein_ETA_1

# With options
pymdmix analyze residence all --cutoff 4.0 --min-time 10
```

### Options

| Option | Default | Description |
|--------|---------|-------------|
| `--cutoff` | 3.5 Å | Distance cutoff for "bound" |
| `--min-time` | 5 ps | Minimum residence to count |
| `--residues` | all | Specific residues to analyze |

---

## Output Files

```
residence/
├── residence_summary.csv    # Per-residue statistics
├── residence_ETA_CT.dat     # Time series for each probe
├── residence_ETA_OH.dat
└── exchange_rates.csv       # Kinetic parameters
```

---

## Understanding Results

### Residence Summary

```csv
residue,probe,mean_residence_ps,max_residence_ps,n_events
ALA_25,ETA_CT,45.2,320.5,127
GLU_42,ETA_OH,128.7,890.2,43
...
```

### Key Metrics

| Metric | Description |
|--------|-------------|
| `mean_residence` | Average time probe stays bound |
| `max_residence` | Longest single binding event |
| `n_events` | Number of binding/unbinding events |
| `occupancy` | Fraction of time with probe bound |

---

## Combining with Energy Data

The most druggable sites show **both**:
- Favorable energy (negative ΔG)
- Long residence times

```python
from pymdmix.analysis import ResidenceAction, HotspotAction

# Get hotspots
hotspots = HotspotAction().run(replica)

# Get residence data
residence = ResidenceAction().run(replica)

# Correlate
for hs in hotspots:
    nearby_residues = hs.protein_contacts
    avg_residence = residence.get_mean(nearby_residues)
    print(f"Hotspot at {hs.center}: ΔG={hs.energy:.1f}, τ={avg_residence:.0f} ps")
```

---

## Python API

```python
from pymdmix.analysis import ResidenceAction
from pymdmix.project import Project

project = Project.load("myproject")
replica = project.get_replica("MyProtein_ETA_1")

action = ResidenceAction(
    cutoff=3.5,
    min_time=5.0
)
results = action.run(replica)

# Per-residue analysis
for res, data in results.by_residue.items():
    if data.mean_residence > 100:  # ps
        print(f"{res}: τ = {data.mean_residence:.0f} ps")

# Plot residence distribution
results.plot_histogram("residence_dist.png")
```

---

## Visualization

### Residence-Colored Structure

```python
from pymdmix.io import write_residence_pdb

write_residence_pdb(
    structure="MyProtein_ETA_1_ref.pdb",
    residence_data=results,
    output="residence_colored.pdb",
    probe="ETA_CT"
)
```

In PyMOL:
```python
load residence_colored.pdb
spectrum b, blue_white_red, minimum=0, maximum=200
```

---

## Applications

1. **Binding site validation** - True sites show long residence
2. **Allosteric sites** - May have different kinetics
3. **Water displacement** - Compare water vs organic probe residence
4. **Selectivity** - Different probes may have different kinetics

---

## Next Steps

- [Hotspot Detection](hotspots.md) - Identify binding sites
- [Tools](../tools.md) - Additional utilities
