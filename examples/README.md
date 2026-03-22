# pyMDMix Example Configuration Files

This directory contains example configuration files for pyMDMix.

## Configuration Files

| File | Purpose | Usage |
|------|---------|-------|
| `project.toml` | Project settings | `mdmix create project -f project.toml` |
| `system.toml` | System definition | `mdmix add system -f system.toml` |
| `replica.toml` | Replica settings | `mdmix add replica -f replica.toml` |
| `mdsettings.toml` | MD parameters | Referenced by replicas |
| `analysis.toml` | Analysis settings | `mdmix analyze --config analysis.toml` |
| `custom_solvent.toml` | Custom solvent | `mdmix create solvent -f custom_solvent.toml` |

## Quick Start

```bash
# 1. Create a project
mdmix create project -n my_protein

# 2. Prepare your protein
mdmix setup prepare protein.pdb -o protein_clean.pdb

# 3. Copy and edit system.toml
cp examples/system.toml my_protein/
# Edit: set pdb = "protein_clean.pdb"

# 4. Add the system
cd my_protein
mdmix add system -f system.toml

# 5. Add replicas
mdmix add replica --system 1abc --solvent ETA --count 3

# 6. Generate simulation inputs
mdmix setup amber --replica 1abc_ETA_1

# 7. Run simulations (external)
# ...

# 8. Analyze
mdmix analyze align all
mdmix analyze density bysolvent -s ETA
mdmix analyze hotspots
```

## File Format

All configuration files use TOML format. See https://toml.io for syntax.

### Key Sections

**project.toml:**
- `[project]` - Name, description
- `[defaults]` - Default simulation parameters
- `[paths]` - Directory structure

**system.toml:**
- `[system]` - Name, PDB, solvent
- `[preparation]` - Cleaning options
- `[solvation]` - Box parameters

**mdsettings.toml:**
- `[mdsettings]` - All MD parameters

**replica.toml:**
- `[replica]` - Name, solvent, system
- `[md]` - Override MD settings
- `[analysis]` - Analysis parameters

**analysis.toml:**
- `[alignment]` - Trajectory alignment
- `[density]` - Density grid calculation
- `[energy]` - Free energy conversion
- `[hotspots]` - Hotspot detection

## See Also

- [User Guide](../docs/user-guide/)
- [Tutorials](../docs/tutorials/)
- [CLI Reference](../docs/api/cli-reference.md)
