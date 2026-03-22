# System Preparation

A pyMDMix project can contain several **systems** - macromolecules ready to be simulated. A **replica** is a system with solvent added and all MD inputs prepared for specific simulation conditions.

## Two Ways to Add a System

1. **Amber Object File (OFF)** - Preferred method. Manually prepare your system as if simulating in vacuo.
2. **PDB File** - Experimental. Use a clean PDB with correct protonation states and Amber-compatible residue names.

---

## Method 1: Using an Amber Object File

This is the recommended approach. Prepare your system using AmberTools (tLeap/xLeap) as if simulating in vacuo.

### Requirements

Your system should be:
- Completely parameterized with your chosen force field
- Correctly protonated
- Capped (ACE/NME) if needed
- Disulfide bridges connected (CYS → CYX)
- Flipped residues corrected

**Important**: Do NOT add solvent or counter-ions. pyMDMix handles this automatically.

### Example: Protein with a Cofactor

For a system containing a NADH cofactor:

```bash
# In tLeap
source leaprc.protein.ff14SB

# Load cofactor parameters
NADH = loadAmberPrep nadh.prep
NADHparams = loadAmberParams nadh.frcmod

# Load protein structure
system = loadPdb 1p44_protonated.pdb

# Check for errors
check system

# Save object file (include all units and parameters)
saveOff system inha.off
saveOff NADH inha.off
saveOff NADHparams inha.off
```

### Verify Your Object File

Before proceeding, test that topology generation works:

```bash
# In tLeap
loadOff inha.off
saveAmberParm system inha.prmtop inha.inpcrd
```

If this fails, fix parameterization errors before continuing.

---

## Method 2: Using a PDB File

> ⚠️ **Experimental Feature**

This method attempts automatic preparation but requires a very clean PDB:

- All residue names must conform to Amber standards
  - `HIE`, `HID`, `HIP` for histidines
  - `ASP`, `ASH` for aspartic acid
  - `CYX` for disulfide cysteines
- No missing residues
- Non-standard residues must be parameterized separately

---

## Adding the System to pyMDMix

### System Configuration File

Create a configuration file (`system.cfg`) to describe your system:

```ini
# pyMDMix System Configuration File

[SYSTEM]
# Unique name for this system
NAME = MyProtein

# Path to amber object file (preferred)
OFF = /path/to/myprotein.off

# OR path to PDB file (experimental)
# PDB = /path/to/myprotein.pdb

# Unit name in the object file (default: first unit found)
# UNAME = system

# Non-standard residues that are part of the solute (not solvent)
# Comma-separated list
# EXTRARES = NADH, HEM

# Force fields or frcmod files needed
# Default: leaprc.ff99SB (use ff14SB or ff19SB for modern work)
# EXTRAFF = leaprc.protein.ff14SB, nadh.frcmod
```

### Configuration Options

| Option | Required | Description |
|--------|----------|-------------|
| `NAME` | Yes | Unique identifier for the system |
| `OFF` | Yes* | Path to Amber object file |
| `PDB` | Yes* | Path to PDB file (*alternative to OFF) |
| `UNAME` | No | Unit name inside the OFF file |
| `EXTRARES` | No | Non-standard residues in the solute |
| `EXTRAFF` | No | Additional force field files |

### Add to Project

Navigate to your project folder and add the system:

```bash
cd myproject
pymdmix add system -f system.cfg
```

Verify it was added:

```bash
pymdmix info project
```

---

## Special Cases

### Multi-chain Proteins

Multi-chain proteins work normally, but pay attention to:
- Imaging during trajectory alignment (may need manual adjustment)
- Chain-specific restraints if needed

### Proteins with Ligands

If your protein contains a bound ligand:

1. Parameterize the ligand (e.g., with `antechamber`)
2. Include ligand unit and parameters in the OFF file
3. List the ligand residue name in `EXTRARES`

```ini
[SYSTEM]
NAME = ProteinLigand
OFF = complex.off
EXTRARES = LIG
EXTRAFF = ligand.frcmod
```

### Membrane Proteins

Membrane proteins require special handling not currently automated in pyMDMix. Consider:
- Pre-embedding in a membrane bilayer
- Using appropriate lipid force fields
- Adjusting solvation procedures

---

## Next Steps

- [Creating Projects](projects.md) - Project organization
- [Creating Replicas](replicas.md) - Set up MD runs
- [Solvent Mixtures](solvents.md) - Choose your probe mixture
