# Solvents

pyMDMix uses organic solvent/water mixtures to probe protein binding sites. Each solvent provides different chemical probes.

## Built-in Solvents

### ETA - Ethanol

Ethanol/water mixture for probing hydroxyl and hydrophobic interactions.

| Property | Value |
|----------|-------|
| Box file | `etawat.off` |
| Residue | ETA |
| Water ratio | ~5:1 water:ethanol |

**Probes:**
- **OH** - Hydroxyl oxygen (hydrogen bond donor/acceptor)
- **CT** - Methyl carbon (hydrophobic)

### MAM - Methylacetamide

Methylacetamide/water mixture for probing amide interactions.

| Property | Value |
|----------|-------|
| Box file | `mamwat.off` |
| Residue | MAM |
| Water ratio | ~5:1 water:MAM |

**Probes:**
- **N** - Amide nitrogen (hydrogen bond donor)
- **O** - Carbonyl oxygen (hydrogen bond acceptor)

### WAT - Water

Pure water reference for normalization and hydration analysis.

| Property | Value |
|----------|-------|
| Box file | `watbox.off` |
| Residue | WAT |

**Probes:**
- **O** - Water oxygen

## Using Solvents

### List Available Solvents

```bash
mdmix info solvents
```

```python
from pymdmix.core import SolventLibrary

library = SolventLibrary()
for name in library.list():
    print(name)
```

### Get Solvent Details

```bash
mdmix info solvents --name ETA
```

```python
from pymdmix.core import SolventLibrary

library = SolventLibrary()
eta = library.get("ETA")

print(f"Name: {eta.name}")
print(f"Residues: {[r.name for r in eta.residues]}")
print(f"Probes: {[p.name for p in eta.probes]}")
```

## Creating Custom Solvents

### 1. Prepare the Solvent Box

Create a solvent box using tleap:

```
# In tleap
source leaprc.gaff
loadamberparams mysolvent.frcmod
loadoff mysolvent.lib
solvateBox MYSOL TIP3PBOX 14.0
saveamberparm MYSOL mysolvent.prmtop mysolvent.inpcrd
saveoff MYSOL mysolvwat.off
```

### 2. Create Solvent Definition

```toml
# mysolvent.toml
[solvent]
name = "MYSOL"
description = "My custom solvent mixture"
box_file = "mysolvwat.off"
unit_name = "MYSOLVWAT"

[[residues]]
name = "MSOL"
atoms = ["C1", "C2", "O1", "H1", "H2", "H3"]
charge = 0.0

[[probes]]
name = "O1"
residue = "MSOL"
atoms = ["O1"]
description = "Custom probe 1"

[[probes]]
name = "CT"
residue = "MSOL"
atoms = ["C1", "H1", "H2", "H3"]
description = "Methyl probe"
```

### 3. Register the Solvent

```bash
mdmix create solvent -f mysolvent.toml
```

Or programmatically:

```python
from pymdmix.core import Solvent, Probe, SolventLibrary

solvent = Solvent(
    name="MYSOL",
    description="My custom solvent",
    box_file="mysolvwat.off",
    residues=[...],
    probes=[
        Probe(name="O1", residue="MSOL", atoms=["O1"]),
        Probe(name="CT", residue="MSOL", atoms=["C1", "H1", "H2", "H3"]),
    ],
)

library = SolventLibrary()
library.add(solvent)
library.save()
```

## Probe Selection

When running analysis, you can select specific probes:

```bash
# Analyze only hydroxyl probes
mdmix analyze density --probes OH

# Analyze multiple probes
mdmix analyze density --probes OH --probes CT
```

```python
from pymdmix.analysis import DensityAction

action = DensityAction()
result = action.run(
    replica=replica,
    probes=["OH", "CT"],
)
```

## Center-of-Mass Probes

For residue-based analysis, you can use center-of-mass (COM) probes:

```toml
[[probes]]
name = "ETA_COM"
residue = "ETA"
atoms = []  # Empty = use COM
description = "Ethanol center of mass"
is_com = true
```

## Best Practices

1. **Choose appropriate solvents** for your target:
   - Kinases: ETA (ATP binding) + MAM (peptide binding)
   - Hydrophobic pockets: ETA
   - Polar sites: MAM

2. **Run multiple replicas** (at least 3) for statistical significance

3. **Use consistent conditions** across solvents for comparison

4. **Normalize densities** when comparing different probes
