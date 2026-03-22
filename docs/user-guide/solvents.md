# Solvent Mixtures

Solvent mixtures are central to the MDMix methodology. pyMDMix includes a library of pre-equilibrated organic solvent/water mixtures ready for use.

## Available Solvents

Query available solvents at any time:

```bash
pymdmix info solvents
```

### Pre-packaged Mixtures

| Solvent | Description | Probes | Use Case |
|---------|-------------|--------|----------|
| **ETA** | Ethanol 20% | CT (hydrophobic), OH (H-bond donor/acceptor), WAT | General screening |
| **MOH** | Methanol 20% | CT, OH, WAT | Smaller H-bond probe |
| **ANT** | Acetonitrile 20% | C, N, WAT | Acceptor-rich regions |
| **ISO** | Isopropanol 20% | CT, OH, WAT | Larger hydrophobic |
| **ISO5** | Isopropanol 5% | CT, OH, WAT | Lower concentration |
| **MAM** | Acetamide 20% | N, O, CT, WAT | Amide interactions |
| **ION** | Ammonium Acetate | Charged probes | Charged hotspots |
| **WAT** | Pure Water | WAT | Control/comparison |

### Probe Types

Each solvent defines **probes** - specific atoms tracked during analysis:

![Probe diagram](../assets/probes.png)

| Chemical Type | Description | Example Probes |
|---------------|-------------|----------------|
| **Hyd** | Hydrophobic | ETA_CT, ISO_CT |
| **Don** | H-bond donor | ETA_OH, MAM_N |
| **Acc** | H-bond acceptor | ANT_N, MAM_O |
| **Don,Acc** | Both donor/acceptor | ETA_OH |
| **Wat** | Water reference | ETA_WAT, ANT_WAT |

---

## Using Solvents

### In Replica Configuration

Specify the solvent when creating replicas:

```ini
[REPLICA]
SYSTEM = MyProtein
SOLVENT = ETA    # Use ethanol mixture
NANOS = 20
```

### Multiple Solvents

Run the same system with different solvents for comprehensive screening:

```bash
# Create replicas for each solvent
for solv in ETA MAM ANT; do
    cat > replica_${solv}.cfg << EOF
[REPLICA]
SYSTEM = MyProtein
SOLVENT = $solv
NANOS = 20
EOF
    pymdmix add replica -f replica_${solv}.cfg
done
```

---

## Creating Custom Solvents

### Workflow

1. **Parameterize** the organic molecule (using antechamber or similar)
2. **Build** the mixture box (using packmol)
3. **Equilibrate** the box (300K, 1 atm)
4. **Create** the OFF file with LEaP
5. **Register** the solvent with pyMDMix

### Step 1: Parameterize the Molecule

```bash
# Using antechamber for a small molecule
antechamber -i molecule.pdb -fi pdb -o molecule.mol2 -fo mol2 \
            -c bcc -s 2 -nc 0

parmchk2 -i molecule.mol2 -f mol2 -o molecule.frcmod
```

### Step 2: Build the Mixture Box

Using packmol to create a 20% mixture:

```
# packmol input
tolerance 2.0
filetype pdb
output mixture_box.pdb

structure water.pdb
  number 2000
  inside box 0. 0. 0. 50. 50. 50.
end structure

structure molecule.pdb
  number 200
  inside box 0. 0. 0. 50. 50. 50.
end structure
```

### Step 3: Equilibrate

Run a short MD equilibration (NPT, 300K) to get proper density.

### Step 4: Create the OFF File

```bash
# In tLeap
source leaprc.water.tip3p
source leaprc.gaff2

# Load molecule parameters
MOL = loadMol2 molecule.mol2
loadAmberParams molecule.frcmod

# Load equilibrated box
box = loadPdb equilibrated_box.pdb
setBox box centers

# Save everything
saveOff box MOLWAT20.off
saveOff MOL MOLWAT20.off
```

### Step 5: Register the Solvent

Create a solvent configuration file (`my_solvent.cfg`):

```ini
[GENERAL]
# Unique name
name = MOL

# Description (use %% to escape %)
info = My Molecule 20%% mixture

# Path to OFF file
objectfile = /path/to/MOLWAT20.off

# Water model used
watermodel = TIP3P

# Box unit name in OFF file
boxunit = MOLWAT20

[PROBES]
# Map probe names to residue@atom
# PROBENAME = RESNAME@ATOMNAMES
WAT = WAT@O
C1 = MOL@C1       # Track carbon 1
O1 = MOL@O1       # Track oxygen

[TYPES]
# Assign chemical types to probes
C1 = Hyd
O1 = Don,Acc
WAT = Wat
```

Add to database:

```bash
pymdmix create solvent -f my_solvent.cfg
```

---

## Solvent Database

### Location

- **System database**: `$PYMDMIX_ROOT/data/solvents/`
- **User database**: `~/.mdmix/` (created automatically if system is read-only)

### JSON Format

Modern pyMDMix uses JSON for solvent definitions:

```json
{
  "name": "ETA",
  "info": "Ethanol 20% mixture",
  "boxunit": "ETAWAT20",
  "objectfile": "ETAWAT20.off",
  "watermodel": "TIP3P",
  "probes": {
    "WAT": {"residue": "WAT", "atoms": ["O"]},
    "CT": {"residue": "ETA", "atoms": ["C1"]},
    "OH": {"residue": "ETA", "atoms": ["O1"]}
  },
  "types": {
    "CT": ["Hyd"],
    "OH": ["Don", "Acc"],
    "WAT": ["Wat"]
  }
}
```

---

## Best Practices

1. **Start with ETA** - Ethanol is versatile and well-characterized
2. **Use multiple solvents** - Different probes reveal different interactions
3. **Match probe types** - Use consistent type names across solvents for combined analysis
4. **Check equilibration** - Ensure mixture box is properly equilibrated before use

---

## Next Steps

- [Creating Replicas](replicas.md) - Use solvents in your simulations
- [Analysis Overview](../analysis/overview.md) - Analyze solvent density
