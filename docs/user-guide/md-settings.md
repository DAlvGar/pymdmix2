# MD Settings

MD settings control simulation parameters: temperature, production length, timestep, restraints, and output frequency.

---

## Default Settings

pyMDMix2 ships with sensible defaults in `data/defaults/md-settings.cfg`:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `nanos` | 20 | Production nanoseconds |
| `temp` | 300 K | Simulation temperature |
| `restrMode` | FREE | No positional restraints |
| `restrForce` | 0.0 | Restraint force constant (kcal/mol·Å²) |
| `md_timestep` | 2.0 fs | Integration timestep |
| `trajfrequency` | 500 | Trajectory write frequency (steps) |
| `prod_steps` | 500 000 | Steps per production file (= 1 ns) |
| `npt_eq_steps` | 500 000 | NPT equilibration steps |

---

## Primary Workflow: [MDSETTINGS] in project.cfg

The recommended way to configure MD settings is through the combined project
configuration file.  Get an annotated template with:

```bash
pymdmix create template
```

Then set parameters in the `[MDSETTINGS]` section:

```ini
[MDSETTINGS]
SOLVENTS    = ETA, MAM, WAT
NREPL       = 3
NANOS       = 20
TEMP        = 300
RESTR       = FREE

# Override defaults
TRAJFREQUENCY = 1000      # snapshot every 2 ps at 2 fs timestep
PROD_STEPS    = 500000    # 1 ns per trajectory file
NPT_EQ_STEPS  = 500000    # 1 ns NPT equilibration
MD_TIMESTEP   = 2.0       # femtoseconds
```

Create the project in one step:

```bash
pymdmix create project -n myproject -f project.cfg
```

---

## Standalone Settings File

You can also use a standalone settings file to add replicas to an existing project:

```ini
[MDSETTINGS]
SOLVENTS = ETA
NREPL    = 2
NANOS    = 40
TEMP     = 300
RESTR    = HA
FORCE    = 0.5
```

```bash
pymdmix add replica -f settings.cfg
```

---

## Restraint Options

### Restraint Modes

| Mode | Description |
|------|-------------|
| `FREE` | No restraints (most common) |
| `BB` | Backbone atoms restrained |
| `HA` | Heavy atoms restrained |

### Example: backbone restraints for flexible loops

```ini
[MDSETTINGS]
SOLVENTS  = ETA
NREPL     = 3
NANOS     = 20
RESTR     = BB
FORCE     = 0.5
RESTRMASK = :1-50    # residues 1-50 only
ALIGNMASK = :1-50    # same region for trajectory alignment
```

### Per-solvent or per-replica overrides

```ini
[MDSETTINGS]
SOLVENTS = ETA, WAT
NREPL    = 3
RESTR    = HA, WAT//FREE    # HA for ETA, FREE for WAT
NANOS    = 20, ETA/1/40    # first ETA replica: 40 ns; rest: 20 ns
```

---

## Per-section Configuration

Use multiple `[MDSETTINGS]` sections with different parameters:

```ini
[MDSETTINGS1]
SOLVENTS = ETA
NREPL    = 3
NANOS    = 20
RESTR    = FREE

[MDSETTINGS2]
SOLVENTS = WAT
NREPL    = 1
NANOS    = 10
RESTR    = FREE
```

---

## Python API

```python
from pymdmix.io.parsers import parse_settings_config_file
from pymdmix.project.settings import MDSettings

# Parse a settings file → list of MDSettings objects
settings = parse_settings_config_file("settings.cfg")
for s in settings:
    print(f"{s.solvent}: {s.nanos} ns, {s.temperature} K, {s.restraint_mode}")

# Or build programmatically
s = MDSettings(
    solvent="ETA",
    nanos=20,
    temperature=300.0,
    restraint_mode="FREE",
    trajectory_frequency=500,
)
```

---

## Overridable Defaults

Any key present in `data/defaults/md-settings.cfg` can be overridden in the
`[MDSETTINGS]` section.  The parser uses fuzzy matching so minor spelling
differences are tolerated.

| Config key | MDSettings field | Default |
|------------|-----------------|---------|
| `TRAJFREQUENCY` | `trajectory_frequency` | 500 |
| `PROD_STEPS` | `production_steps` | 500 000 |
| `NPT_EQ_STEPS` | `npt_equilibration_steps` | 500 000 |
| `MD_TIMESTEP` | `timestep` | 2.0 |
| `MINSTEPS` | `minimization_steps` | 5 000 |
| `HEATING_STEPS` | `heating_steps` | 100 000 |
| `MDPROGRAM` | `md_program` | `AMBER` |
| `MDFOLDER` | `md_folder` | `md` |

---

## Best Practices

1. **Use defaults first** — they work for most cases
2. **Match equilibration and production restraints** — avoid sudden changes
3. **Monitor temperature/pressure** — check equilibration succeeded
4. **Use per-section configs** — for mixed restraint strategies

---

## Next Steps

- [Queue Scripts](queue-scripts.md) — HPC submission
- [Creating Replicas](replicas.md) — apply settings to replicas
- [User Guide: Projects](projects.md) — full project workflow
