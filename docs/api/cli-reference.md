# CLI Reference

Complete reference for the pyMDMix command-line interface.

## Usage

```bash
pymdmix [OPTIONS] COMMAND [ARGS]...

# Or via Python module
python -m pymdmix.cli [OPTIONS] COMMAND [ARGS]...
```

## Global Options

| Option | Description |
|--------|-------------|
| `--version` | Show version and exit |
| `-v, --verbose` | Enable verbose output |
| `--help` | Show help message |

---

## Commands Overview

| Command | Description |
|---------|-------------|
| `create` | Create projects, replicas, or solvents |
| `add` | Add systems, replicas, or groups to project |
| `info` | Display information about project components |
| `analyze` | Run analysis on simulation data |
| `queue` | Generate and manage queue scripts |
| `plot` | Generate plots from analysis data |
| `tools` | Utility tools for grids and structures |

---

## create

Create new pyMDMix entities.

### create project

```bash
pymdmix create project NAME [OPTIONS]
```

| Option | Description |
|--------|-------------|
| `NAME` | Project directory name |
| `-f, --force` | Overwrite existing project |

### create solvent

```bash
pymdmix create solvent -f CONFIG [OPTIONS]
```

| Option | Description |
|--------|-------------|
| `-f, --file` | Solvent configuration file |
| `--validate` | Validate without adding |

---

## add

Add components to an existing project.

### add system

```bash
pymdmix add system -f CONFIG [OPTIONS]
```

| Option | Description |
|--------|-------------|
| `-f, --file` | System configuration file |

### add replica

```bash
pymdmix add replica -f CONFIG [OPTIONS]
```

| Option | Description |
|--------|-------------|
| `-f, --file` | Replica configuration file |
| `--count N` | Create N replicas |

### add group

```bash
pymdmix add group -n NAME -s REPLICAS [OPTIONS]
```

| Option | Description |
|--------|-------------|
| `-n, --name` | Group name |
| `-s, --selection` | Replica names to include |

---

## info

Display information about project components.

### info project

```bash
pymdmix info project
```

### info systems

```bash
pymdmix info systems [--detailed]
```

### info replicas

```bash
pymdmix info replicas [--detailed]
```

### info solvents

```bash
pymdmix info solvents [--detailed]
```

### info analysis

```bash
pymdmix info analysis REPLICA_NAME
```

---

## analyze

Run analysis on simulation trajectories.

### Common Selection Syntax

```bash
pymdmix analyze COMMAND SELECTION [OPTIONS]

# Selection types:
all                           # All replicas
bysolvent -s SOLV1 SOLV2      # By solvent name
byname -s REP1 REP2           # By replica name
group -s GROUPNAME            # By group
```

### Common Options

| Option | Description |
|--------|-------------|
| `-N, --nanoselect` | Nanosecond range (e.g., 1:10) |
| `-C, --ncpus` | Number of CPUs |

### analyze align

```bash
pymdmix analyze align SELECTION [OPTIONS]
```

| Option | Description |
|--------|-------------|
| `--mask` | Alignment mask (residue range) |
| `--ref` | Reference PDB file |

### analyze density

```bash
pymdmix analyze density SELECTION [OPTIONS]
```

| Option | Description |
|--------|-------------|
| `--spacing` | Grid spacing (Å) |
| `--buffer` | Grid buffer (Å) |

### analyze energy

```bash
pymdmix analyze energy SELECTION [OPTIONS]
```

| Option | Description |
|--------|-------------|
| `--temp` | Temperature (K) |
| `--average` | Average across replicas |

### analyze hotspots

```bash
pymdmix analyze hotspots SELECTION [OPTIONS]
```

| Option | Description |
|--------|-------------|
| `--threshold` | Energy cutoff (kcal/mol) |
| `--min-size` | Minimum cluster size |

### analyze residence

```bash
pymdmix analyze residence SELECTION [OPTIONS]
```

| Option | Description |
|--------|-------------|
| `--cutoff` | Distance cutoff (Å) |
| `--min-time` | Minimum residence (ps) |

---

## queue

Manage queue submission scripts.

### queue generate

```bash
pymdmix queue generate SELECTION [OPTIONS]
```

| Option | Description |
|--------|-------------|
| `--config` | Queue configuration file |
| `--template` | Custom template file |
| `--array` | Generate array job |

### queue submit

```bash
pymdmix queue submit SELECTION [OPTIONS]
```

| Option | Description |
|--------|-------------|
| `--dry-run` | Show commands without executing |
| `--stage` | Submission stage (eq, prod) |

---

## plot

Generate visualization plots.

### plot rmsd

```bash
pymdmix plot rmsd SELECTION [OPTIONS]
```

| Option | Description |
|--------|-------------|
| `-o, --output` | Output file |
| `--format` | Output format (png, pdf, svg) |

### plot energy

```bash
pymdmix plot energy SELECTION [OPTIONS]
```

### plot density

```bash
pymdmix plot density SELECTION [OPTIONS]
```

---

## tools

Utility tools for data manipulation.

### tools grid-info

```bash
pymdmix tools grid-info FILE
```

### tools grid-math

```bash
pymdmix tools grid-math OPERATION FILES... -o OUTPUT
```

Operations: `add`, `average`, `min`, `max`, `scale`

### tools convert

```bash
pymdmix tools convert INPUT -o OUTPUT
```

### tools combine-hotspots

```bash
pymdmix tools combine-hotspots FILES... -o OUTPUT [--distance D]
```

### tools extract-frames

```bash
pymdmix tools extract-frames TRAJECTORY [OPTIONS]
```

---

## Configuration Files

### System Configuration

```ini
[SYSTEM]
NAME = string       # Required
OFF = path          # Required (or PDB)
PDB = path          # Alternative to OFF
UNAME = string      # Unit name in OFF
EXTRARES = list     # Non-standard residues
EXTRAFF = list      # Force field files
```

### Replica Configuration

```ini
[REPLICA]
SYSTEM = string     # Required
SOLVENT = string    # Required
NANOS = int         # Production length
RESTRMODE = string  # FREE, BB, HA, CUSTOM
RESTRAINMASK = mask # If RESTRMODE=CUSTOM
RESTRFORCE = float  # Force constant
MDSETTINGS = path   # Custom MD settings
```

### Solvent Configuration

```ini
[GENERAL]
name = string
info = string
objectfile = path
boxunit = string
watermodel = string

[PROBES]
PROBENAME = RESNAME@ATOMS

[TYPES]
PROBENAME = TYPE1,TYPE2
```

---

## Environment Variables

| Variable | Description |
|----------|-------------|
| `AMBERHOME` | Amber installation path |
| `MDMIX_DATA` | Custom data directory |
| `MDMIX_DEBUG` | Enable debug output |

---

## Exit Codes

| Code | Meaning |
|------|---------|
| 0 | Success |
| 1 | General error |
| 2 | Invalid arguments |
| 3 | File not found |
| 4 | Analysis error |
