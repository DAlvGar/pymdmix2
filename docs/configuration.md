# Configuration Files

pyMDMix uses TOML configuration files for projects, systems, and replicas.

## Project Configuration

```toml
# project.toml
[project]
name = "my_protein_study"
description = "Binding site analysis of MyProtein"

[defaults]
temperature = 300.0
nanos = 20
md_program = "AMBER"
```

## System Configuration

```toml
# system.toml
[system]
name = "1abc"
pdb = "protein.pdb"
solvent = "ETA"

[preparation]
cap_termini = true
detect_disulfides = true
remove_waters = true
extra_residues = []  # Non-standard residues to keep

[solvation]
buffer = 14.0  # Box buffer in Angstroms
```

## MD Settings

```toml
# mdsettings.toml
[mdsettings]
solvent = "ETA"
nanos = 20
temperature = 300.0
timestep = 2.0  # femtoseconds

# Restraints
restraint_mode = "FREE"  # FREE, BB (backbone), HA (heavy atoms), CUSTOM
restraint_force = 0.0
restraint_mask = ""

# MD engine
md_program = "AMBER"  # AMBER, OPENMM, GROMACS, NAMD
production_ensemble = "NVT"  # NVT or NPT

# Step counts
trajectory_frequency = 500
minimization_steps = 5000
heating_steps = 100000
npt_equilibration_steps = 500000
production_steps = 500000

# Other
heating_initial_temp = 100.0
npt_pressure = 1.0
use_netcdf = true
wrap_coordinates = true

# Force fields
force_fields = ["leaprc.protein.ff14SB", "leaprc.gaff"]
```

## Replica Configuration

```toml
# replica.toml
[replica]
name = "replica_ETA_1"
solvent = "ETA"

[md]
nanos = 20
temperature = 300.0
restraint_mode = "FREE"

[analysis]
align_mask = ":1-150@CA"  # Amber mask for alignment
density_probes = ["OH", "CT"]  # Probes to analyze
```

## Solvent Definition

Custom solvents can be defined:

```toml
# solvent.toml
[solvent]
name = "ETA"
description = "Ethanol"
box_file = "etawat.off"
unit_name = "ETAWAT"

[[residues]]
name = "ETA"
atoms = ["C1", "C2", "O", "H1", "H2", "H3", "H4", "H5", "HO"]
charge = 0.0

[[probes]]
name = "OH"
residue = "ETA"
atoms = ["O", "HO"]
description = "Hydroxyl probe"

[[probes]]
name = "CT"
residue = "ETA"
atoms = ["C2", "H4", "H5"]
description = "Methyl probe"
```

## Analysis Configuration

```toml
# analysis.toml
[density]
spacing = 0.5  # Grid spacing in Angstroms
cutoff = 10.0  # Distance from protein

[hotspots]
energy_cutoff = -0.5  # kcal/mol
cluster_distance = 2.0  # Angstroms
min_points = 3

[energy]
temperature = 300.0
reference_density = 1.0  # Reference for free energy
```

## Environment Variables

| Variable | Description | Default |
|----------|-------------|---------|
| `AMBERHOME` | Path to AmberTools | Auto-detect |
| `MDMIX_DATA` | Custom solvent library | `~/.mdmix/solvents` |
| `MDMIX_NCPUS` | Default CPU count | 1 |

## Legacy Configuration Support

pyMDMix 2.0 can read legacy `.cfg` files from pyMDMix 1.x:

```python
from pymdmix.project import MDSettings

# Load from legacy config
settings = MDSettings.from_legacy_cfg("old_settings.cfg")
```
