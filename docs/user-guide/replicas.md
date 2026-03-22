# Creating Replicas

A **replica** is an independent MD simulation of a system solvated in a specific mixture. Each replica contains:

- Solvated topology and coordinates
- MD input files (minimization, equilibration, production)
- Queue submission scripts
- Reference structures for analysis

## Replica Configuration

### Basic Configuration

Create a replica configuration file (`replica.cfg`):

```ini
[REPLICA]
# System to use (must be added to project first)
SYSTEM = MyProtein

# Solvent mixture (from solvent database)
SOLVENT = ETA

# Production simulation length in nanoseconds
NANOS = 20
```

### Full Configuration Options

```ini
[REPLICA]
# REQUIRED
SYSTEM = MyProtein
SOLVENT = ETA

# SIMULATION LENGTH
NANOS = 20          # Production nanoseconds
NPROD = 20          # Number of production steps (default: same as NANOS)

# RESTRAINTS
RESTRMODE = FREE    # FREE, BB (backbone), HA (heavy atoms), CUSTOM
RESTRAINMASK =      # Custom mask if RESTRMODE = CUSTOM
RESTRFORCE = 5.0    # Force constant (kcal/mol/Å²)

# MD SETTINGS (override defaults)
MDSETTINGS = /path/to/mdsettings.cfg

# ALIGNMENT
ALIGNMASK = 1-200   # Residue mask for trajectory alignment

# MULTIPLE REPLICAS
COUNT = 3           # Create this many replicas
```

---

## Creating Replicas

### Single Replica

```bash
pymdmix add replica -f replica.cfg
```

### Multiple Replicas

```bash
# Using COUNT in config
pymdmix add replica -f replica.cfg

# Or via command line
pymdmix add replica -f replica.cfg --count 3
```

This creates `MyProtein_ETA_1`, `MyProtein_ETA_2`, `MyProtein_ETA_3`.

---

## Restraint Modes

| Mode | Description | Use Case |
|------|-------------|----------|
| `FREE` | No restraints | Most common, allows full flexibility |
| `BB` | Backbone restrained | Study side-chain dynamics |
| `HA` | Heavy atoms restrained | Minimal structural deviation |
| `CUSTOM` | Custom atom mask | Specific requirements |

### Custom Restraints

```ini
[REPLICA]
SYSTEM = MyProtein
SOLVENT = ETA
RESTRMODE = CUSTOM
RESTRAINMASK = :1-50@CA,C,N  # Restrain backbone of residues 1-50
RESTRFORCE = 2.0
```

---

## Replica Folder Structure

After creation, each replica folder contains:

```
MyProtein_ETA_1/
├── MyProtein_ETA_1.prmtop        # Solvated topology
├── MyProtein_ETA_1.inpcrd        # Starting coordinates
├── MyProtein_ETA_1_ref.pdb       # Reference structure
├── COMMANDS.sh                   # Execution script
├── queue.sh                      # Queue submission script
├── min/                          # Minimization inputs
│   ├── min1.in
│   └── min2.in
├── eq/                           # Equilibration inputs
│   ├── eq1.in
│   └── eq2.in
├── md/                           # Production inputs
│   ├── md1.in
│   ├── md2.in
│   └── ...
└── align/                        # Analysis outputs (created later)
```

---

## Running Simulations

### View Commands

```bash
cd MyProtein_ETA_1/
cat COMMANDS.sh
```

### Local Execution

```bash
# Run directly (for testing)
bash COMMANDS.sh
```

### HPC Submission

Edit the queue script for your cluster:

```bash
# Edit queue.sh for your system (SLURM, PBS, SGE, etc.)
vim queue.sh

# Submit
sbatch queue.sh   # SLURM
qsub queue.sh     # PBS/SGE
```

### Monitoring

```bash
# Check output files
ls md/*.out

# View RMSD (after alignment)
pymdmix plot rmsd byname -s MyProtein_ETA_1
```

---

## Post-Simulation

After simulations complete:

1. **Copy trajectories** back to the project folder
2. **Verify file names** match expectations (md1.nc, md2.nc, etc.)
3. **Check outputs** for errors

```bash
# Verify expected files exist
ls MyProtein_ETA_1/md/*.nc

# Should see: md1.nc, md2.nc, ... md20.nc (for 20 ns)
```

---

## Python API

```python
from pymdmix.project import Project, Replica

# Load project
project = Project.load("myproject")

# Create replica programmatically
replica = Replica(
    system="MyProtein",
    solvent="ETA",
    nanos=20,
    restr_mode="FREE"
)

project.add_replica(replica)
project.save()

# Access replica info
replica = project.get_replica("MyProtein_ETA_1")
print(replica.topology_path)
print(replica.status)
```

---

## Troubleshooting

### LEaP Errors

If solvation fails, check:
- System OFF file is valid
- Solvent boxunit exists
- Force field compatibility

### Missing Parameters

Ensure all non-standard residues have parameters in the OFF file.

### Box Too Small

If protein is too close to box edges, increase the buffer:

```ini
[REPLICA]
BOXBUFFER = 12.0   # Angstroms (default: 10)
```

---

## Next Steps

- [MD Settings](md-settings.md) - Customize simulation parameters
- [Queue Scripts](queue-scripts.md) - HPC submission
- [Analysis Overview](../analysis/overview.md) - Analyze your simulations
