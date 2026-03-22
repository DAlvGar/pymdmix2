# Hotspot Detection

Hotspots are localized regions with favorable binding energy. pyMDMix identifies and ranks these sites automatically from energy grids.

## Theory

Hotspot detection involves:

1. **Thresholding** - Find regions below energy cutoff
2. **Clustering** - Group nearby favorable points
3. **Ranking** - Order by total binding energy
4. **Characterization** - Identify chemical preferences

---

## Running Hotspot Detection

### Basic Usage

```bash
# Detect hotspots in all replicas
pymdmix analyze hotspots all

# For specific solvent
pymdmix analyze hotspots bysolvent -s ETA
```

### With Options

```bash
# Custom energy threshold (default: -1.0 kcal/mol)
pymdmix analyze hotspots all --threshold -1.5

# Minimum cluster size
pymdmix analyze hotspots all --min-size 3

# Maximum number of hotspots to report
pymdmix analyze hotspots all --max-hotspots 20
```

---

## Output Files

```
hotspots/
├── hotspots.pdb         # Hotspot centers as pseudo-atoms
├── hotspots.csv         # Ranked list with energies
├── hotspots_CT.dx       # Per-probe contribution maps
├── hotspots_OH.dx
└── summary.txt          # Human-readable summary
```

### Hotspots PDB

The `hotspots.pdb` file contains pseudo-atoms at hotspot centers:

- **B-factor** = Energy (more negative = stronger)
- **Occupancy** = Cluster size
- **Residue name** = Dominant probe type (HYD, DON, ACC, etc.)

---

## Hotspot Properties

Each hotspot includes:

| Property | Description |
|----------|-------------|
| `center` | XYZ coordinates of centroid |
| `energy` | Total binding energy (kcal/mol) |
| `size` | Number of grid points |
| `probes` | Contributing probe types |
| `contacts` | Nearby protein residues |

---

## Interpreting Results

### Energy Values

| Energy (kcal/mol) | Druggability |
|-------------------|--------------|
| < -3.0 | Highly druggable |
| -3.0 to -2.0 | Druggable |
| -2.0 to -1.0 | Potentially druggable |
| > -1.0 | Challenging target |

### Probe Composition

Ideal drug-binding sites typically show:
- Mixed probe types (hydrophobic + polar)
- Multiple converging hotspots
- Enclosure by protein atoms

---

## Python API

```python
from pymdmix.analysis import HotspotAction
from pymdmix.project import Project

# Detect hotspots
project = Project.load("myproject")
replica = project.get_replica("MyProtein_ETA_1")

action = HotspotAction(
    threshold=-1.0,
    min_size=3,
    max_hotspots=20
)
hotspots = action.run(replica)

# Analyze results
for i, hs in enumerate(hotspots):
    print(f"\nHotspot {i+1}:")
    print(f"  Center: {hs.center}")
    print(f"  Energy: {hs.energy:.2f} kcal/mol")
    print(f"  Size: {hs.size} points")
    print(f"  Probes: {hs.probe_contributions}")
    print(f"  Contacts: {hs.protein_contacts}")

# Save results
action.write_pdb(hotspots, "hotspots.pdb")
action.write_csv(hotspots, "hotspots.csv")
```

---

## Combining Across Solvents

Merge hotspots from different solvent simulations:

```python
from pymdmix.analysis import HotspotAction

# Get hotspots from each solvent
eta_hotspots = action.run(project.get_replica("MyProtein_ETA_1"))
mam_hotspots = action.run(project.get_replica("MyProtein_MAM_1"))

# Merge nearby hotspots
merged = HotspotAction.merge_hotspots(
    [eta_hotspots, mam_hotspots],
    distance_cutoff=3.0  # Å
)

# Rank by combined energy
merged.sort(key=lambda h: h.energy)
```

---

## Clustering Methods

### Default: Grid-based Clustering

Fast, deterministic clustering based on connected grid points.

### Alternative: DBSCAN

Density-based clustering for better shape detection:

```python
action = HotspotAction(
    clustering="dbscan",
    eps=2.0,      # Maximum point distance
    min_samples=5  # Minimum points per cluster
)
```

### Surface Distance Clustering

For identifying pocket-specific hotspots:

```python
# Requires surface calculation
action = HotspotAction(
    clustering="surface",
    surface_distance_cutoff=6.0
)
```

---

## Visualization

### PyMOL Script

```python
# Load protein and hotspots
load MyProtein_ETA_1_ref.pdb, protein
load hotspots/hotspots.pdb, hotspots

# Color by energy (B-factor)
spectrum b, red_white_green, hotspots, minimum=-3, maximum=0

# Show as spheres
show spheres, hotspots
set sphere_scale, 0.5, hotspots

# Label top hotspots
label hotspots and name CA, "%.1f" % b
```

### Generate Publication Figure

```python
from pymdmix.io import plot_hotspots

plot_hotspots(
    protein="MyProtein_ETA_1_ref.pdb",
    hotspots="hotspots/hotspots.pdb",
    energy_grid="grids/ETA_CT_DG.dx",
    output="figure_hotspots.png",
    view="pocket"  # or "global"
)
```

---

## Druggability Assessment

Score binding site druggability:

```python
from pymdmix.analysis import DruggabilityScore

scorer = DruggabilityScore()
score = scorer.calculate(hotspots)

print(f"Druggability Score: {score.value:.2f}")
print(f"Classification: {score.classification}")
print(f"Confidence: {score.confidence:.0%}")
```

### Score Components

- **Hotspot energy** - Total binding energy available
- **Hotspot clustering** - Proximity of multiple hotspots
- **Probe diversity** - Mix of chemical types
- **Enclosure** - Degree of protein burial

---

## Troubleshooting

### No Hotspots Detected

- Lower the threshold (e.g., -0.5 kcal/mol)
- Check energy grids exist
- Verify probe density was sufficient

### Too Many Small Hotspots

- Increase `min_size` parameter
- Raise energy threshold
- Apply smoothing to energy grids

### Hotspots at Surface Only

- Normal for flat proteins
- Look for convergence of probe types
- Consider surface clustering method

---

## Next Steps

- [Residence Analysis](residence.md) - Kinetics of binding
- [Tools](../tools.md) - Additional analysis utilities
