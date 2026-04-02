# Projects and Systems

A **project** is the top-level container that organises your MDMix work.  It
holds one or more **systems** (proteins) and multiple **replicas**
(independent simulations per solvent).

---

## Primary Workflow: Single Project Config File

The recommended way to create a project is with a single `project.cfg` file
that combines the `[SYSTEM]` and `[MDSETTINGS]` sections — exactly as in the
original pyMDMix.

### 1. Get the template

```bash
pymdmix create template            # → project.cfg in current directory
pymdmix create template -o sim.cfg # → custom path
```

### 2. Edit the template

```ini
[SYSTEM]
NAME    = MyProtein
OFF     = myprotein.off   # Amber Object File

[MDSETTINGS]
SOLVENTS = ETA, MAM, WAT  # solvent mixture names
NREPL    = 3              # replicas per solvent
NANOS    = 20             # production nanoseconds
TEMP     = 300            # temperature (K)
RESTR    = FREE           # no restraints
```

### 3. Create the project

```bash
pymdmix create project -n myproject -f project.cfg
```

Output:

```
Creating project 'myproject' in /path/to/myproject
✓ Project 'myproject' created
  System:   MyProtein
  Solvents: ETA, MAM, WAT
  Replicas: 9
    - MyProtein_ETA_1
    - MyProtein_ETA_2
    - MyProtein_ETA_3
    - MyProtein_MAM_1
    ...
```

---

## Project Structure

```
myproject/
├── .mdmix/                    # Project metadata
│   └── project.json
├── inputs/                    # Config files
│   └── project.cfg
├── MyProtein_ETA_1/          # Replica folders
│   ├── md/                   # MD outputs
│   ├── align/                # Aligned trajectories
│   ├── grids/                # Density grids
│   └── ...
├── MyProtein_ETA_2/
└── MyProtein_MAM_1/
```

---

## Step-by-step Alternative

If you prefer to build the project incrementally:

```bash
# Create empty project
pymdmix create project -n myproject

# Add system
cat > system.cfg << EOF
[SYSTEM]
NAME = MyProtein
OFF  = myprotein.off
EOF
pymdmix add system -f system.cfg

# Add replicas from a settings file
cat > settings.cfg << EOF
[MDSETTINGS]
SOLVENTS = ETA, MAM
NREPL    = 3
NANOS    = 20
RESTR    = FREE
EOF
pymdmix add replica -f settings.cfg
```

---

## Project Information

```bash
pymdmix info project     # project summary
pymdmix info replicas    # detailed replica listing
pymdmix info systems     # registered systems
pymdmix info solvents    # available solvents
```

---

## Managing Replicas

### Adding replicas (MDSETTINGS format — primary)

```bash
# Add 3 replicas for ETA and MAM from a settings file
cat > more_replicas.cfg << EOF
[MDSETTINGS]
SOLVENTS = ETA, MAM
NREPL    = 3
NANOS    = 20
EOF
pymdmix add replica -f more_replicas.cfg
```

### Adding replicas (REPLICA format — legacy)

```bash
cat > replica.cfg << EOF
[REPLICA]
SYSTEM    = MyProtein
SOLVENT   = ETA
NANOS     = 20
RESTRMODE = FREE
EOF
pymdmix add replica -f replica.cfg --count 3
```

### Replica Naming

Replicas are automatically named: `{SYSTEM}_{SOLVENT}_{NUMBER}`

Example: `MyProtein_ETA_1`, `MyProtein_ETA_2`, `MyProtein_MAM_1`

---

## Replica Groups

Groups allow you to organise replicas for batch operations.

```bash
# Create a group for all ethanol replicas
pymdmix add group -n ethanol_runs -s MyProtein_ETA_1 MyProtein_ETA_2

# Use in analysis
pymdmix analyze density group -s ethanol_runs
```

---

## Python API

```python
from pymdmix.io.parsers import parse_project_config
from pymdmix.project import Project

# --- Primary: parse and inspect a project config file ----------------------
cfg = parse_project_config("project.cfg")
print(cfg.system.name)        # → 'MyProtein'
print(len(cfg.settings))      # → 9  (3 solvents × 3 replicas)
for s in cfg.settings:
    print(f"  {s.solvent}: {s.nanos} ns, {s.temperature} K")

# --- Load and work with an existing project ---------------------------------
project = Project.load("myproject")

for replica in project.replicas:
    print(f"{replica.name}: {replica.state.name}")

# Get replica by name
replica = project.get_replica("MyProtein_ETA_1")

# Filter by solvent
eta_reps = project.get_replicas_by_solvent("ETA")
```

---

## Best Practices

1. **Use a single project.cfg** — the easiest way to reproduce a setup
2. **Multiple replicas** — run 3–5 replicas per solvent for statistics
3. **Descriptive system names** — use meaningful identifiers
4. **Use groups** — organise replicas for batch processing
5. **Back up project.json** — it stores the full project state

---

## Next Steps

- [Quick Start](../quickstart.md) — step-by-step tutorial
- [MD Settings](md-settings.md) — simulation parameter reference
- [Replicas](replicas.md) — replica lifecycle and queue submission
