# pyMDMix Example Configuration Files

This directory contains example configuration files for pyMDMix.

## Configuration Files

| File | Purpose | Status |
|------|---------|--------|
| `system.toml` | System definition | ✅ `pymdmix add system -f system.cfg` |
| `replica.toml` | Replica settings | ✅ `pymdmix add replica -f replica.cfg` |
| `mdsettings.toml` | MD parameters | ✅ Referenced by replicas |
| `analysis.toml` | Analysis settings | ✅ `pymdmix analyze --config analysis.toml` |
| `project.toml` | Project defaults | 🚧 Planned - not yet implemented |
| `custom_solvent.toml` | Custom solvent | 🚧 Planned |

> **Note:** System and replica configs currently use CFG (INI) format.
> TOML examples show the planned unified format.

## Quick Start

```bash
# 1. Create a project
pymdmix create project -n my_protein

# 2. Prepare your protein
pymdmix setup prepare protein.pdb -o protein_clean.pdb

# 3. Create system config (system.cfg)
cat > system.cfg << EOF
[SYSTEM]
NAME = myprotein
PDB = protein_clean.pdb
EOF

# 4. Add the system
pymdmix add system -f system.cfg

# 5. Create replica config (replica.cfg)
cat > replica.cfg << EOF
[REPLICA]
SYSTEM = myprotein
SOLVENT = ETA
NANOS = 50
RESTRMODE = HA
RESTRFORCE = 0.01
EOF

# 6. Add replicas
pymdmix add replica -f replica.cfg --count 3

# 7. Run simulations (external - Amber/OpenMM)
# ...

# 8. Analyze
pymdmix analyze align all
pymdmix analyze density all
pymdmix analyze energy all
pymdmix analyze hotspots all
```

## File Formats

### Currently Supported: CFG (INI-style)

**system.cfg:**
```ini
[SYSTEM]
NAME = myprotein
PDB = protein_clean.pdb
OFF = protein.off        # Alternative: Amber OFF file
EXTRARES = LIG,CYX       # Non-standard residues
```

**replica.cfg:**
```ini
[REPLICA]
SYSTEM = myprotein
SOLVENT = ETA
NANOS = 50
RESTRMODE = HA           # FREE | BB | HA | CUSTOM
RESTRFORCE = 0.01        # kcal/mol·Å²
RESTRAINMASK = :1-100    # Only for CUSTOM mode
```

### Planned: TOML Format

The `.toml` files in this directory show the planned unified format.
See each file for structure examples.

## See Also

- [User Guide](../docs/user-guide/)
- [Tutorials](../docs/tutorials/)
