# MD Settings

MD settings control simulation parameters: temperature, pressure, timestep, restraints, and output frequency.

## Default Settings

pyMDMix uses sensible defaults for organic mixture simulations:

| Parameter | Default | Description |
|-----------|---------|-------------|
| Temperature | 300 K | Production temperature |
| Pressure | 1 atm | NPT ensemble |
| Timestep | 2 fs | Integration step |
| Cutoff | 10 Å | Nonbonded cutoff |
| PME | Yes | Electrostatics |
| SHAKE | Yes | Constrain H-bonds |
| Output | 1 ns | Trajectory file per ns |

---

## Customizing Settings

### Settings Configuration File

Create `mdsettings.cfg`:

```ini
[MD]
# Temperature in Kelvin
TEMP = 300

# Pressure in atm (0 = NVT ensemble)
PRES = 1.0

# Timestep in femtoseconds
TIMESTEP = 2.0

# Nonbonded cutoff in Angstroms
CUTOFF = 10.0

# Trajectory output frequency (steps)
TRAJFREQ = 5000

# Energy output frequency (steps)
OUTFREQ = 500

# Restart frequency (steps)
RESTFREQ = 5000

[MINIMIZATION]
# Minimization steps
MINSTEPS = 5000

# Minimization method (sd, cg)
MINMETHOD = sd

[EQUILIBRATION]
# Equilibration length in ps
EQLENGTH = 500

# Restraint force constant during equilibration
EQRESTRAINT = 5.0
```

### Apply to Replica

```ini
[REPLICA]
SYSTEM = MyProtein
SOLVENT = ETA
NANOS = 20
MDSETTINGS = /path/to/mdsettings.cfg
```

---

## Restraint Options

### Restraint Modes

| Mode | Description |
|------|-------------|
| `FREE` | No restraints (most common) |
| `BB` | Backbone atoms restrained |
| `HA` | Heavy atoms restrained |
| `CUSTOM` | User-defined mask |

### Restraint Configuration

```ini
[REPLICA]
SYSTEM = MyProtein
SOLVENT = ETA
RESTRMODE = CUSTOM
RESTRAINMASK = :1-50@CA,C,N,O
RESTRFORCE = 2.0  # kcal/mol/Å²
```

### Staged Restraint Release

For equilibration with gradual restraint release:

```ini
[EQUILIBRATION]
# Restraint values for each eq stage
RESTRAINTS = 5.0, 2.0, 1.0, 0.5, 0.0
# Length of each stage in ps
EQSTAGES = 100, 100, 100, 100, 100
```

---

## Production Settings

### Simulation Length

```ini
[REPLICA]
NANOS = 20        # Total production nanoseconds
NPROD = 20        # Number of production files (1 per ns default)
```

### Output Frequency

```ini
[MD]
# Write trajectory every N steps
TRAJFREQ = 5000   # = 10 ps at 2 fs timestep

# Write energies every N steps
OUTFREQ = 500     # = 1 ps

# Write restart every N steps
RESTFREQ = 5000   # = 10 ps
```

---

## Engine-Specific Settings

### Amber

```ini
[AMBER]
# Use pmemd.cuda for GPU
EXECUTABLE = pmemd.cuda

# Hydrogen mass repartitioning (allows 4 fs timestep)
HMASS = yes

# If HMASS = yes, timestep can be increased
TIMESTEP = 4.0
```

### OpenMM

```ini
[OPENMM]
# Platform selection
PLATFORM = CUDA

# Precision
PRECISION = mixed

# GPU device
DEVICEID = 0
```

### GROMACS

```ini
[GROMACS]
# Number of MPI ranks
NRANKS = 4

# Number of OpenMP threads
NTHREADS = 8
```

---

## Python API

```python
from pymdmix.project import Config, MDSettings

# Load default settings
settings = MDSettings.default()

# Modify
settings.temperature = 310
settings.timestep = 4.0
settings.hmass = True

# Save
settings.save("my_settings.cfg")

# Use with replica
from pymdmix.project import Replica
replica = Replica(
    system="MyProtein",
    solvent="ETA",
    md_settings=settings
)
```

---

## Common Modifications

### High Temperature

```ini
[MD]
TEMP = 350  # For enhanced sampling
```

### Longer Simulations

```ini
[REPLICA]
NANOS = 100
NPROD = 100  # 100 x 1ns files
```

### GPU Acceleration

```ini
[AMBER]
EXECUTABLE = pmemd.cuda
HMASS = yes
TIMESTEP = 4.0
```

---

## Best Practices

1. **Use defaults first** - They work for most cases
2. **H-mass repartitioning** - Enables 4 fs timestep safely
3. **Match equilibration restraints** - To replica restraint mode
4. **Monitor temperature/pressure** - Check equilibration succeeded
5. **Save restarts frequently** - For recovery from crashes

---

## Next Steps

- [Queue Scripts](queue-scripts.md) - HPC submission
- [Creating Replicas](replicas.md) - Apply settings
