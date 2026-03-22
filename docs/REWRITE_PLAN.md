# pyMDMix Rewrite Plan

**Decision:** Full rewrite rather than migration.
**Reason:** Biskit integration is too deep (private method access, 15+ PDBModel methods, subclassing). A compatibility layer would essentially reimplement Biskit.

---

## Architecture Principles

### 1. Use Libraries Idiomatically
- **parmed** for structure manipulation (don't wrap it)
- **MDAnalysis** for trajectory reading (multi-format support)
- **NumPy** for grids and numerical work
- Standard library for everything else

### 2. Modern Python
- Python 3.10+ (match statement, type unions)
- Type hints on all public APIs
- Dataclasses for data containers
- Pathlib for paths
- f-strings everywhere

### 3. Design for Flexibility
- Abstract trajectory interface (Amber, Gromacs, NAMD, generic)
- Plugin system for analysis actions
- Clear separation: setup vs simulation vs analysis

### 4. Keep What Works
- Core algorithms (density, residence, hotspots)
- Solvent library concept
- Project/Replica organization
- Grid format (DX)

---

## Module Structure

```
pymdmix/                    # lowercase, PEP8
├── __init__.py
├── cli.py                  # Click-based CLI
├── core/
│   ├── __init__.py
│   ├── structure.py        # Structure handling (parmed utilities)
│   ├── trajectory.py       # Abstract trajectory interface
│   ├── grid.py             # 3D grid class
│   └── solvent.py          # Solvent definitions
├── setup/
│   ├── __init__.py
│   ├── prepare.py          # Structure preparation (capping, protonation)
│   ├── solvate.py          # Solvation with mixtures
│   └── topology.py         # Topology generation
├── engines/
│   ├── __init__.py
│   ├── base.py             # Abstract MD engine
│   ├── amber.py            # Amber implementation
│   ├── openmm.py           # OpenMM implementation
│   └── namd.py             # NAMD implementation
├── analysis/
│   ├── __init__.py
│   ├── base.py             # Abstract analysis action
│   ├── density.py          # Density grid calculation
│   ├── residence.py        # Residence time analysis
│   └── hotspots.py         # Binding hotspot detection
├── io/
│   ├── __init__.py
│   ├── readers.py          # Trajectory readers (MDAnalysis-based)
│   ├── writers.py          # Output writers
│   └── formats.py          # File format handlers
├── project/
│   ├── __init__.py
│   ├── project.py          # Project management
│   ├── replica.py          # Replica handling
│   └── config.py           # Configuration parsing
└── data/
    ├── solvents/           # Solvent definitions
    ├── templates/          # Queue templates
    └── forcefields/        # FF patches if needed
```

---

## Core Classes

### Grid (core/grid.py)

```python
@dataclass
class Grid:
    """3D grid for storing volumetric data."""
    data: NDArray[np.float64]
    origin: tuple[float, float, float]
    spacing: float

    @classmethod
    def from_structure(cls, coords: NDArray, spacing: float = 0.5, padding: float = 5.0) -> Grid:
        """Create grid encompassing structure with padding."""
        ...

    def add_count(self, coord: NDArray) -> None:
        """Add count at coordinate position."""
        ...

    def to_density(self, n_frames: int) -> Grid:
        """Convert counts to density."""
        ...

    def write_dx(self, path: Path) -> None:
        """Write in OpenDX format."""
        ...
```

### Trajectory (core/trajectory.py)

```python
class TrajectoryReader(Protocol):
    """Protocol for trajectory readers."""

    @property
    def n_frames(self) -> int: ...
    def __iter__(self) -> Iterator[Frame]: ...

class MDAnalysisReader:
    """Universal reader using MDAnalysis (Amber, Gromacs, NAMD, etc.)"""
    ...

class AmberNetCDFReader:
    """Direct Amber NetCDF reader (minimal dependency fallback)"""
    ...

def open_trajectory(topology: Path, trajectory: Path, backend: str = 'auto') -> TrajectoryReader:
    """Open trajectory with appropriate reader."""
    ...
```

### Analysis Plugin System (analysis/base.py)

```python
class Action(ABC):
    """Abstract base for analysis actions."""
    name: str

    @abstractmethod
    def run(self, replica: Replica, trajectory: TrajectoryReader, **kwargs) -> ActionResult:
        ...

@register_action("density")
class DensityAction(Action):
    """Calculate probe density grids."""
    ...
```

---

## Timeline

| Week | Focus | Deliverable |
|------|-------|-------------|
| 1 | Core modules | Working trajectory → density pipeline |
| 2 | Setup + engines | Full Amber workflow |
| 3 | Polish + docs | Releasable v0.3.0 |

---

## Phase 1: Core Infrastructure (This Week)

### Tasks
1. [ ] Set up new `pymdmix/` package structure
2. [ ] Implement `core/grid.py` (DX read/write, density counting)
3. [ ] Implement `core/trajectory.py` (MDAnalysis + Amber fallback)
4. [ ] Implement `core/structure.py` (parmed utilities)
5. [ ] Implement `core/solvent.py` (JSON-based definitions)
6. [ ] Unit tests for all core modules
7. [ ] Basic CLI skeleton with Click

### Success Criteria
- Can read Amber trajectory
- Can calculate density grid
- Can write DX file
- All tests pass

---

## Dependencies

```toml
[project]
dependencies = [
    "numpy>=1.21",
    "scipy>=1.7",
    "parmed>=4.0",
    "click>=8.0",
]

[project.optional-dependencies]
full = [
    "MDAnalysis>=2.0",
    "matplotlib>=3.5",
]
```
