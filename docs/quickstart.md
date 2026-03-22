# Quickstart Tutorial

This tutorial walks through a complete pyMDMix workflow: from protein preparation to hotspot identification.

## Prerequisites

- Python 3.10+
- AmberTools installed (`tleap`, `cpptraj`)
- A protein structure (PDB file)

## Installation

```bash
pip install pymdmix
```

## Step 1: Create a Project

```bash
# Create a new project
mdmix create project -n my_protein

cd my_protein
```

This creates the project structure:
```
my_protein/
├── .mdmix/
│   └── project.toml
├── systems/
├── replicas/
└── analysis/
```

## Step 2: Prepare the Protein

```bash
# Clean and prepare the PDB
mdmix setup prepare protein.pdb -o protein_clean.pdb

# This will:
# - Add ACE/NME caps to termini
# - Detect disulfide bonds
# - Remove waters and ions
# - Fix residue naming for Amber
```

## Step 3: Add a System

Create a system configuration:

```toml
# system.toml
[system]
name = "1abc"
pdb = "protein_clean.pdb"
solvent = "ETA"

[solvation]
buffer = 14.0
```

Add it to the project:

```bash
mdmix add system -f system.toml
```

## Step 4: Create Replicas

```bash
# Create 3 replicas with ethanol
mdmix add replica --system 1abc --solvent ETA --count 3 --nanos 20
```

This creates:
```
replicas/
├── 1abc_ETA_1/
│   ├── min/
│   ├── eq/
│   └── md/
├── 1abc_ETA_2/
└── 1abc_ETA_3/
```

## Step 5: Run Simulations

Generate input files:

```bash
# For Amber
mdmix setup amber --replica 1abc_ETA_1

# For OpenMM
mdmix setup openmm --replica 1abc_ETA_1
```

Run the simulation (example with Amber):

```bash
cd replicas/1abc_ETA_1
./run_minimization.sh
./run_equilibration.sh
./run_production.sh
```

## Step 6: Align Trajectories

After simulations complete:

```bash
mdmix analyze align --replica 1abc_ETA_1
# Or align all replicas:
mdmix analyze align all
```

## Step 7: Calculate Density Grids

```bash
# Calculate density for all probes
mdmix analyze density --replica 1abc_ETA_1

# Or by solvent across all replicas
mdmix analyze density bysolvent -s ETA
```

This produces:
- `OH_density.dx` - Hydroxyl probe density
- `CT_density.dx` - Methyl probe density

## Step 8: Detect Hotspots

```bash
mdmix analyze hotspots --energy-cutoff -0.5 --cluster-distance 2.0

# Output:
# - hotspots.json
# - hotspots.pdb
# - hotspots_summary.txt
```

## Step 9: Visualize Results

Load the density grids and hotspots in PyMOL or VMD:

```python
# In PyMOL
load protein.pdb
load OH_density.dx, OH_density
isosurface OH_surface, OH_density, -1.0
load hotspots.pdb
```

## Python API

You can also use pyMDMix as a Python library:

```python
from pymdmix.core import Grid, load_structure
from pymdmix.analysis import HotspotAction, HotSpotSet

# Load density grid
density = Grid.read_dx("OH_density.dx")

# Convert to free energy
dg = density.to_free_energy(temperature=300.0)

# Detect hotspots
action = HotspotAction()
result = action.run(
    grids={"OH": density},
    energy_cutoff=-0.5,
    cluster_distance=2.0,
)

# Work with hotspot set
from pymdmix.analysis.hotspots import HotSpotSet
hs_set = HotSpotSet(probe="OH", name="hydroxyl")
hs_set.add_hotspots(hotspots)

# Filter
filtered = hs_set.prune_by_energy(-1.0)

# Cluster
representatives = hs_set.get_cluster_representatives(cutoff=2.5)
```

## Next Steps

- See [Configuration](configuration.md) for detailed settings
- See [Solvents](solvents.md) for available solvent mixtures
- See [Analysis](analysis.md) for advanced analysis options
- See [CLI Reference](cli.md) for all commands
