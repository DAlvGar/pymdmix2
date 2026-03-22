# Projects and Systems

A **project** is a container that organizes your MDMix work. It can contain multiple **systems** (proteins) and multiple **replicas** (independent simulations) for each system.

## Project Structure

```
myproject/
├── .mdmix/                    # Project metadata
│   └── project.json
├── systems/                   # System definitions
│   └── MyProtein.json
├── MyProtein_ETA_1/          # Replica folders
│   ├── md/                   # MD outputs
│   ├── align/                # Aligned trajectories
│   ├── grids/                # Density grids
│   └── ...
├── MyProtein_ETA_2/
└── MyProtein_MAM_1/
```

---

## Creating a Project

```bash
# Create a new project
pymdmix create project myproject

# Or initialize in current directory
mkdir myproject && cd myproject
pymdmix create project .
```

## Project Information

```bash
# View project summary
pymdmix info project

# Detailed replica listing
pymdmix info replicas

# System information
pymdmix info systems
```

---

## Managing Systems

### Adding a System

Create a system configuration file and add it:

```bash
# Create system.cfg (see System Preparation)
pymdmix add system -f system.cfg
```

### Listing Systems

```bash
pymdmix info systems
```

### Removing a System

To remove a system, delete its configuration file and any associated replica folders:

```bash
rm systems/MyProtein.json
rm -rf MyProtein_*  # Be careful!
```

---

## Managing Replicas

### Creating Replicas

```bash
# From configuration file
pymdmix add replica -f replica.cfg

# Multiple replicas
pymdmix add replica -f replica.cfg --count 3
```

### Listing Replicas

```bash
pymdmix info replicas
```

### Replica Naming

Replicas are automatically named: `{SYSTEM}_{SOLVENT}_{NUMBER}`

Example: `MyProtein_ETA_1`, `MyProtein_ETA_2`, `MyProtein_MAM_1`

---

## Replica Groups

Groups allow you to organize replicas for batch operations.

### Creating Groups

```bash
# Create a group for all ethanol replicas
pymdmix add group -n ethanol_runs -s MyProtein_ETA_1 MyProtein_ETA_2
```

### Using Groups in Analysis

```bash
# Analyze only replicas in a group
pymdmix analyze density group -s ethanol_runs
```

---

## Working with Multiple Projects

Each project is independent. Switch between projects by changing directories:

```bash
cd /path/to/project1
pymdmix info project

cd /path/to/project2
pymdmix info project
```

---

## Python API

```python
from pymdmix.project import Project, Config

# Load existing project
project = Project.load("myproject")

# Access replicas
for replica in project.replicas:
    print(f"{replica.name}: {replica.status}")

# Get specific replica
replica = project.get_replica("MyProtein_ETA_1")

# Project configuration
config = Config.load()
print(config.amber_home)
```

---

## Best Practices

1. **One project per study** - Keep related simulations together
2. **Descriptive system names** - Use meaningful identifiers
3. **Multiple replicas** - Run 3-5 replicas per condition for statistics
4. **Use groups** - Organize replicas for batch processing
5. **Backup your project** - Especially the `.mdmix/` folder and replica configurations

---

## Next Steps

- [Creating Replicas](replicas.md) - Configure simulation replicas
- [MD Settings](md-settings.md) - Customize simulation parameters
