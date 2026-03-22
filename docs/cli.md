# CLI Reference

Complete command-line reference for pyMDMix.

## Global Options

```bash
mdmix [OPTIONS] COMMAND
```

| Option | Description |
|--------|-------------|
| `-v, --verbose` | Enable verbose output |
| `--debug` | Enable debug output |
| `--version` | Show version |
| `--help` | Show help |

---

## create

Create new resources.

### create project

```bash
mdmix create project [OPTIONS]
```

| Option | Description | Default |
|--------|-------------|---------|
| `-n, --name` | Project name | `mdmix_project` |
| `-f, --config` | Configuration file | - |
| `-d, --directory` | Project directory | `./<name>` |
| `--force` | Overwrite existing | False |

### create solvent

```bash
mdmix create solvent -f <config> [OPTIONS]
```

| Option | Description |
|--------|-------------|
| `-f, --file` | Solvent config file (required) |
| `--validate` | Validate only, don't add |

---

## add

Add resources to a project.

### add system

```bash
mdmix add system -f <config> [OPTIONS]
```

| Option | Description | Default |
|--------|-------------|---------|
| `-f, --file` | System config file (required) | - |
| `-p, --project` | Project directory | `.` |

### add replica

```bash
mdmix add replica -f <config> [OPTIONS]
```

| Option | Description | Default |
|--------|-------------|---------|
| `-f, --file` | Replica config file (required) | - |
| `-p, --project` | Project directory | `.` |
| `--count` | Number of replicas | 1 |

### add group

```bash
mdmix add group -n <name> -s <replica> [OPTIONS]
```

| Option | Description |
|--------|-------------|
| `-n, --name` | Group name (required) |
| `-s, --selection` | Replica names (repeatable) |
| `-p, --project` | Project directory |

---

## setup

Prepare structures and inputs.

### setup prepare

```bash
mdmix setup prepare <structure.pdb> [OPTIONS]
```

| Option | Description | Default |
|--------|-------------|---------|
| `-o, --output` | Output PDB file | `<input>_clean.pdb` |
| `--cap/--no-cap` | Add ACE/NME caps | True |
| `--disulfide/--no-disulfide` | Detect disulfides | True |
| `--remove-water/--keep-water` | Remove waters | True |

### setup solvate

```bash
mdmix setup solvate <structure.pdb> -s <solvent> [OPTIONS]
```

| Option | Description | Default |
|--------|-------------|---------|
| `-s, --solvent` | Solvent name (required) | - |
| `-b, --buffer` | Box buffer (Å) | 14.0 |
| `-o, --output` | Output prefix | - |

### setup amber

```bash
mdmix setup amber --replica <name> [OPTIONS]
```

Generate Amber input files for a replica.

### setup openmm

```bash
mdmix setup openmm --replica <name> [OPTIONS]
```

Generate OpenMM input files for a replica.

---

## analyze

Run analysis actions.

### analyze align

```bash
mdmix analyze align [OPTIONS] [REPLICAS]
```

| Option | Description | Default |
|--------|-------------|---------|
| `--mask` | Alignment mask | Auto (backbone) |
| `--reference` | Reference PDB | Replica ref |
| `all` | Align all replicas | - |

### analyze density

```bash
mdmix analyze density [OPTIONS]
```

| Option | Description | Default |
|--------|-------------|---------|
| `--replica` | Replica name | - |
| `--probes` | Probe names | All |
| `--spacing` | Grid spacing (Å) | 0.5 |
| `bysolvent -s <name>` | Analyze by solvent | - |

### analyze energy

```bash
mdmix analyze energy [OPTIONS]
```

Convert density grids to free energy.

| Option | Description | Default |
|--------|-------------|---------|
| `--temperature` | Temperature (K) | 300.0 |
| `--input` | Input density grid | - |
| `--output` | Output energy grid | - |

### analyze hotspots

```bash
mdmix analyze hotspots [OPTIONS]
```

| Option | Description | Default |
|--------|-------------|---------|
| `--energy-cutoff` | Energy threshold (kcal/mol) | -0.5 |
| `--cluster-distance` | Clustering distance (Å) | 2.0 |
| `--min-points` | Minimum cluster points | 3 |
| `--output-dir` | Output directory | `.` |

### analyze residence

```bash
mdmix analyze residence [OPTIONS]
```

Calculate residence times at hotspots.

---

## info

Display information.

### info project

```bash
mdmix info project [-p <path>]
```

Show project summary.

### info systems

```bash
mdmix info systems [-p <path>]
```

List systems in project.

### info replicas

```bash
mdmix info replicas [-p <path>] [--system <name>]
```

List replicas with status.

### info solvents

```bash
mdmix info solvents [--name <name>]
```

List available solvents or show details.

### info analysis

```bash
mdmix info analysis [-p <path>]
```

Show analysis status.

---

## Examples

```bash
# Full workflow
mdmix create project -n kinase_study
mdmix setup prepare kinase.pdb -o kinase_clean.pdb
mdmix add system -f system.toml
mdmix add replica --system kinase --solvent ETA --count 3

# After running simulations...
mdmix analyze align all
mdmix analyze density bysolvent -s ETA
mdmix analyze hotspots --energy-cutoff -0.7

# Check status
mdmix info replicas
mdmix info analysis
```
