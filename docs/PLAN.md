# pyMDMix 2.0 - Migration Plan

**Goal**: Port pyMDMix from Python 2.7 to Python 3.10+, removing Biskit dependency.

**Strategy**: Wrap ParmEd (not extract Biskit). ParmEd is actively maintained, Python 3 native, and already handles all the PDB/topology operations we need. Our `core/structure.py` provides a clean abstraction layer (~350 lines) that replaces Biskit's `PDBModel` (~3000 lines).

**Key decision** (2026-03-22): We evaluated extracting Biskit classes vs wrapping ParmEd. ParmEd wins because:
1. Zero maintenance burden (someone else maintains it)
2. Already integrated with AmberTools ecosystem
3. Handles Amber prmtop/rst7 natively (Biskit couldn't)
4. Our thin wrapper is easier to test and extend

---

## Repository Structure

```
~/clawd/pymdmix2/          ← CANONICAL REPO (all new development)
~/clawd/pyMDmix-port/      ← REFERENCE ONLY (legacy Python 2.7 code)
    └── pyMDMix/           ← Original source to read from
```

---

## Biskit Dependency Removal - COMPREHENSIVE AUDIT

### Decision: Wrap ParmEd + MDAnalysis ✅

**Why NOT extract Biskit:**
- Biskit is ~15K lines of Python 2.7 code
- `PDBModel` alone is ~3K lines with complex internal state
- Would require extensive porting, testing, maintenance
- Creates another codebase we'd have to maintain forever

**Why ParmEd + MDAnalysis works:**
- ParmEd: PDB/topology handling, atom data, saving (v4.3.1, actively maintained)
- MDAnalysis: Trajectory reading, alignment/fitting, RMSD (v2.x, industry standard)
- Both are Python 3 native, well-tested, widely used in computational chemistry
- Our thin wrapper (~500 lines) is easier to test and extend

---

### COMPLETE Biskit Method Audit

**Legend:** ✅ = implemented | 🔲 = TODO | ❌ = not needed

#### 1. PDBModel - Dictionary Access (column-based data)

Legacy pattern: `model['residue_name']`, `model['serial_number']`, `model['temperature_factor']`

| Biskit Column | Usage | pymdmix2 | Status |
|---------------|-------|----------|--------|
| `['residue_name']` | Get residue names per atom | `[a.residue.name for a in struct.atoms]` | ✅ |
| `['residue_number']` | Get residue numbers | `[a.residue.number for a in struct.atoms]` | ✅ |
| `['serial_number']` | Atom serial numbers | `[a.idx for a in struct.atoms]` | ✅ |
| `['name']` | Atom names | `[a.name for a in struct.atoms]` | ✅ |
| `['chain_id']` | Chain IDs | `[a.residue.chain for a in struct.atoms]` | ✅ |
| `['temperature_factor']` | B-factors (for restraints) | `[a.bfactor for a in struct.atoms]` | ✅ |
| `['element']` | Element symbols | `[a.element_name for a in struct.atoms]` | ✅ |

**Note:** ParmEd stores these as atom/residue attributes, not dict columns. Access pattern differs but data is available.

#### 2. PDBModel - Masking Methods

| Biskit Method | Usage | pymdmix2 Replacement | Status |
|---------------|-------|----------------------|--------|
| `maskProtein()` | Protein atoms | `get_protein_mask(struct)` | ✅ |
| `maskNA()` | Nucleic acids | `get_nucleic_mask(struct)` | ✅ |
| `maskH2O()` | Waters | `get_water_mask(struct)` | ✅ |
| `maskH()` | Hydrogen atoms | `get_hydrogen_mask(struct)` | 🔲 |
| `maskHeavy()` | Non-hydrogen | `get_heavy_atom_mask(struct)` | ✅ |
| `maskBB()` | Backbone (N,CA,C,O) | `get_backbone_mask(struct)` | 🔲 |
| `maskCA()` | Alpha carbons | `get_atom_mask(struct, 'CA')` | ✅ |
| `maskFrom(col, vals)` | Generic selection | `get_residue_mask()` / `get_atom_mask()` | ✅ |

#### 3. PDBModel - Transformation & Fitting (⚠️ CRITICAL)

| Biskit Method | Usage | pymdmix2 Replacement | Status |
|---------------|-------|----------------------|--------|
| `magicFit(ref, mask)` | Superimpose to reference | `align_structures(mobile, ref, mask)` | 🔲 |
| `transform(matrix)` | Apply rotation/translation | MDAnalysis `translate()` / `rotate()` | 🔲 |
| `center()` | Get geometric center | `struct.coordinates.mean(axis=0)` | ✅ |

**Implementation for `magicFit`:**
```python
# Use MDAnalysis for alignment (add to structure.py)
from MDAnalysis.analysis import align

def align_structures(
    mobile: parmed.Structure,
    reference: parmed.Structure, 
    mask: np.ndarray | None = None
) -> tuple[parmed.Structure, float]:
    """
    Superimpose mobile onto reference (Kabsch algorithm).
    
    Returns aligned structure and RMSD.
    Uses MDAnalysis internally for robust fitting.
    """
    # Convert to MDAnalysis Universe for alignment
    # Apply transformation back to parmed Structure
    ...
```

#### 4. PDBModel - Selection & Subsetting

| Biskit Method | Usage | pymdmix2 Replacement | Status |
|---------------|-------|----------------------|--------|
| `compress(mask)` | Select atoms by bool mask | `select_atoms(struct, mask)` | ✅ |
| `take(indices)` | Select atoms by indices | `struct[indices]` (parmed native) | ✅ |
| `clone()` | Deep copy | `struct.copy(parmed.Structure)` | ✅ |
| `takeChains(indices)` | Select specific chains | `select_chains(struct, chain_ids)` | 🔲 |
| `atomRange()` | All atom indices | `range(len(struct.atoms))` | ✅ |

#### 5. PDBModel - Chain Operations

| Biskit Method | Usage | pymdmix2 Replacement | Status |
|---------------|-------|----------------------|--------|
| `lenChains()` | Number of chains | `len(get_chain_ids(struct))` | ✅ |
| `chainIndex()` | Start index of each chain | `get_chain_start_indices(struct)` | 🔲 |
| `takeChains(ids)` | Extract chains | `select_chains(struct, ids)` | 🔲 |

#### 6. PDBModel - Profile Methods

| Biskit Method | Usage | pymdmix2 Replacement | Status |
|---------------|-------|----------------------|--------|
| `atom2resProfile(col)` | Per-atom → per-residue | `atom_to_residue_property(struct, prop)` | 🔲 |
| `res2atomProfile(col)` | Per-residue → per-atom | `residue_to_atom_property(struct, prop)` | 🔲 |

**Implementation:**
```python
def atom_to_residue_property(struct, prop_name: str) -> list:
    """Get one value per residue from atom property."""
    return [getattr(res.atoms[0], prop_name) for res in struct.residues]

def residue_to_atom_property(struct, values: list) -> np.ndarray:
    """Expand per-residue values to per-atom array."""
    result = np.zeros(len(struct.atoms))
    for i, res in enumerate(struct.residues):
        for atom in res.atoms:
            result[atom.idx] = values[i]
    return result
```

#### 7. PDBModel - I/O

| Biskit Method | Usage | pymdmix2 Replacement | Status |
|---------------|-------|----------------------|--------|
| `writePdb(path, amber=1)` | Save PDB (Amber format) | `save_structure(struct, path)` | ✅ |
| `source` | Original file path | `struct.title` or track separately | ✅ |
| `PDB_KEYS` | Valid column names | Not needed (parmed has attrs) | ❌ |

#### 8. PDBModel - Properties

| Biskit Property | Usage | pymdmix2 Replacement | Status |
|-----------------|-------|----------------------|--------|
| `xyz` | Coordinates (N,3) | `struct.coordinates` | ✅ |
| `lenAtoms()` | Number of atoms | `len(struct.atoms)` | ✅ |
| `lenResidues()` | Number of residues | `len(struct.residues)` | ✅ |

#### 9. AmberParmBuilder (ACE/NME capping, SS bonds)

| Biskit Method | Usage | pymdmix2 Replacement | Status |
|---------------|-------|----------------------|--------|
| `__ssBonds(model, cutoff)` | Detect disulfides | `find_disulfides(struct, cutoff)` | ✅ |
| `__cys2cyx(model, ss)` | Rename CYS→CYX | `rename_cys_to_cyx(struct, ss)` | ✅ |
| `__fLines(template, ss)` | Generate leap SS commands | `generate_ss_leap_commands(ss_pairs)` | 🔲 |
| ACE/NME capping | Add terminal caps | tleap-based in `engines/amber.py` | ✅ |

#### 10. Biskit.tools

| Function | Usage | pymdmix2 Replacement | Status |
|----------|-------|----------------------|--------|
| `absfile(path)` | Absolute path | `Path(path).resolve()` | ✅ |
| `projectRoot()` | Find project root | Not needed (explicit paths) | ❌ |
| `dump/load` | Pickle I/O | `pickle.dump/load` or avoid | ⚠️ |

#### 11. Biskit.SettingsParser

| Feature | pymdmix2 | Status |
|---------|----------|--------|
| Config file parsing | `utils/settings_parser.py` | ✅ |
| Environment detection | Simplified in settings | ✅ |

#### 12. Other Biskit modules

| Module | pymdmix2 | Status |
|--------|----------|--------|
| `Biskit.test` | `pytest` | ✅ |
| `Biskit.LogFile` | `logging` (stdlib) | ✅ |
| `Biskit.mathUtils` | `numpy` | ✅ |
| `Biskit.PDBParseFile` | `parmed.load_file()` | ✅ |

---

### Summary: What's Missing

**Must implement (blocking):**
1. `get_hydrogen_mask()` - trivial, filter by atomic_number == 1
2. `get_backbone_mask()` - filter by atom name in {N, CA, C, O}
3. `align_structures()` - use MDAnalysis.analysis.align
4. `select_chains()` - filter by chain ID

**Should implement (helpful):**
5. `get_chain_start_indices()` - for chain iteration
6. `atom_to_residue_property()` / `residue_to_atom_property()` - profile conversion
7. `generate_ss_leap_commands()` - for Amber prep

**Don't need:**
- `PDB_KEYS` - parmed uses attributes, not dict
- `projectRoot()` - we use explicit paths
- Complex Biskit internal state management

---

## Module Migration Status

### ✅ COMPLETE (32 modules)

| Legacy | New Location | Notes |
|--------|--------------|-------|
| `Align.py` | `analysis/align.py` | - |
| `Amber.py` | `engines/amber.py` | LeapSession, AmberChecker, AmberWriter added |
| `AutoPrepare.py` | `setup/prepare.py` | PDB2PQR interface |
| `Browser.py` | `project/browser.py` | Filesystem navigation |
| `containers.py` | `core/containers.py` | Atom, Residue, Probe |
| `Energy.py` | `analysis/energy.py` | Density→ΔG conversion |
| `Executor.py` | `engines/executor.py` | Parallel job execution |
| `GridData.py` | `core/grid.py` | DX/XPLOR I/O, radial ops |
| `GridsManager.py` | `core/grid.py` | Merged into Grid class |
| `NamdDCDParser.py` | `io/dcd_parser.py` | DCD trajectory reader |
| `NAMD.py` | `engines/namd.py` | NAMD input generation |
| `OFFManager.py` | `io/off_manager.py` | Amber OFF file handling |
| `OpenMM.py` | `engines/openmm.py` | OpenMM interface |
| `OpenMMUtils.py` | `engines/openmm.py` | Merged |
| `Parsers.py` | `io/parsers.py` | Config file parsing |
| `Plotter.py` | `io/plotting.py` | Visualization |
| `Projects.py` | `project/project.py` | Project container |
| `QueueWriting.py` | `engines/queue.py` | SLURM/PBS/SGE scripts |
| `Replicas.py` | `project/replica.py` | Replica management |
| `SettingsParser.py` | `utils/settings_parser.py` | Config management |
| `settings.py` | (simplified) | No Biskit, env detection |
| `Solvents.py` | `core/solvent.py` | JSON-based solvent library |
| `Structures.py` | `core/structure.py` | parmed-based structure ops |
| `Systems.py` | `core/system.py` | System class |
| `tools.py` | `utils/tools.py` | Utility functions |
| `Trajectory.py` | `core/trajectory.py` | MDAnalysis-based |
| - | `cli.py` | New Click-based CLI |
| - | `engines/gromacs.py` | New GROMACS support |
| - | `setup/solvate.py` | Solvation workflow |
| - | `analysis/density.py` | Density calculation |
| - | `analysis/residence.py` | Residence time analysis |
| - | `analysis/hotspots.py` | Basic hotspot detection |

### ❌ NOT YET PORTED (4 modules)

| Legacy | Status | Priority |
|--------|--------|----------|
| `MDSettings.py` | Not started | P1 |
| `PDB.py` | Not started | P2 |
| `Analysis.py` | Partial (ActionsManager missing) | P1 |
| `HotSpotsManager.py` | Partial (clustering missing) | P1 |

---

## Implementation Tasks

### Phase 1: Complete Core (P0) ✅ DONE
- [x] Grid with DX/XPLOR I/O
- [x] Replica with folder creation, status checks
- [x] Amber with LeapSession, AmberChecker
- [x] CLI with all commands wired

### Phase 2: Missing Modules (P1)

#### 2.1 MDSettings Class
**Source**: `pyMDmix-port/pyMDMix/MDSettings.py`
**Target**: `pymdmix2/pymdmix/project/settings.py`

```python
@dataclass
class MDSettings:
    """MD simulation parameters."""
    solvent: str
    nanos: int = 20
    temperature: float = 300.0
    timestep: float = 2.0  # fs
    restraint_mode: str = "FREE"  # FREE, BB, HA, CUSTOM
    restraint_force: float = 0.0
    restraint_mask: str = ""
    align_mask: str = ""
    # ... parsing from config files
```

#### 2.2 ActionsManager
**Source**: `pyMDmix-port/pyMDMix/Analysis.py`
**Target**: `pymdmix2/pymdmix/analysis/manager.py`

```python
class ActionsManager:
    """Run analysis actions on multiple replicas in parallel."""
    def __init__(self, ncpus: int = 1):
        ...
    def add_replicas(self, replicas: list[Replica]):
        ...
    def add_actions(self, actions: list[Action]):
        ...
    def run(self):
        # ThreadPoolExecutor-based parallel execution
```

#### 2.3 HotSpotSet & Clustering
**Source**: `pyMDmix-port/pyMDMix/HotSpotsManager.py`
**Target**: `pymdmix2/pymdmix/analysis/hotspots.py` (expand)

```python
@dataclass
class HotSpotSet:
    """Collection of hotspots with clustering."""
    hotspots: list[Hotspot]
    
    def cluster(self, cutoff: float = 2.5) -> list[HotSpotCluster]:
        # scipy.cluster.hierarchy
        ...
    
    def prune_by_energy(self, max_energy: float):
        ...
    
    def prune_by_volume(self, min_volume: float):
        ...
```

### Phase 3: Complete structure.py (P2)

#### 3.1 Missing Mask Functions

```python
# Add to core/structure.py:

def get_hydrogen_mask(struct: parmed.Structure) -> NDArray[np.bool_]:
    """Get boolean mask for hydrogen atoms."""
    return np.array([a.atomic_number == 1 for a in struct.atoms])

def get_backbone_mask(struct: parmed.Structure) -> NDArray[np.bool_]:
    """Get boolean mask for protein backbone atoms (N, CA, C, O)."""
    bb_names = {'N', 'CA', 'C', 'O'}
    protein = get_protein_mask(struct)
    names = get_atom_mask(struct, bb_names)
    return protein & names
```

#### 3.2 Chain Operations

```python
def get_chain_start_indices(struct: parmed.Structure) -> list[int]:
    """Get atom index where each chain starts."""
    indices = [0]
    current_chain = struct.atoms[0].residue.chain if struct.atoms else None
    for i, atom in enumerate(struct.atoms):
        if atom.residue.chain != current_chain:
            indices.append(i)
            current_chain = atom.residue.chain
    return indices

def select_chains(struct: parmed.Structure, chain_ids: Sequence[str]) -> parmed.Structure:
    """Select atoms belonging to specified chains."""
    chain_set = set(chain_ids)
    mask = np.array([a.residue.chain in chain_set for a in struct.atoms])
    return select_atoms(struct, mask)
```

#### 3.3 Structure Alignment (uses MDAnalysis)

```python
def align_structures(
    mobile: parmed.Structure,
    reference: parmed.Structure,
    mask: NDArray[np.bool_] | None = None,
) -> tuple[parmed.Structure, float]:
    """
    Superimpose mobile structure onto reference (Kabsch algorithm).
    
    Parameters
    ----------
    mobile : parmed.Structure
        Structure to align (will be modified in place)
    reference : parmed.Structure  
        Reference structure
    mask : NDArray[np.bool_] | None
        Atoms to use for fitting. If None, use all atoms.
        
    Returns
    -------
    tuple[parmed.Structure, float]
        Aligned structure and RMSD
        
    Notes
    -----
    Uses MDAnalysis internally for robust Kabsch superposition.
    Equivalent to Biskit's magicFit().
    """
    try:
        import MDAnalysis as mda
        from MDAnalysis.analysis import align as mda_align
    except ImportError:
        raise ImportError("MDAnalysis required for alignment: pip install MDAnalysis")
    
    # Create temporary MDAnalysis universes
    # ... implementation using mda_align.alignto()
```

#### 3.4 Profile Conversion

```python
def atom_to_residue_values(struct: parmed.Structure, atom_values: Sequence) -> list:
    """
    Convert per-atom values to per-residue (takes first atom of each residue).
    Equivalent to Biskit's atom2resProfile.
    """
    result = []
    for res in struct.residues:
        first_atom_idx = res.atoms[0].idx
        result.append(atom_values[first_atom_idx])
    return result

def residue_to_atom_values(struct: parmed.Structure, res_values: Sequence) -> NDArray:
    """
    Expand per-residue values to per-atom array.
    Equivalent to Biskit's res2atomProfile.
    """
    result = np.zeros(len(struct.atoms))
    for i, res in enumerate(struct.residues):
        for atom in res.atoms:
            result[atom.idx] = res_values[i]
    return result
```

#### 3.5 Probe/Density Helpers (for analysis)

```python
def detect_solvent_type(
    struct: parmed.Structure, 
    known_solvents: dict[str, list[str]]
) -> str | None:
    """Auto-detect solvent mixture from residue composition."""
    # Get organic solvent residues (not water, not ions, not protein/NA)
    organic_mask = ~(get_water_mask(struct) | get_ion_mask(struct) | get_solute_mask(struct))
    resnames = set(a.residue.name for i, a in enumerate(struct.atoms) if organic_mask[i])
    
    for solvent_name, residue_names in known_solvents.items():
        if resnames == set(residue_names):
            return solvent_name
    return None

def get_probe_coords(
    struct: parmed.Structure, 
    residue_name: str, 
    atom_names: list[str]
) -> NDArray[np.float64]:
    """Extract coordinates of probe atoms for density analysis."""
    mask = get_residue_mask(struct, residue_name) & get_atom_mask(struct, atom_names)
    return struct.coordinates[mask]

def iter_residue_coords(
    struct: parmed.Structure,
    residue_name: str,
) -> Iterator[NDArray[np.float64]]:
    """Iterate over coordinates of each residue with given name."""
    for res in struct.residues:
        if res.name == residue_name:
            indices = [a.idx for a in res.atoms]
            yield struct.coordinates[indices]

def get_residue_com_coords(
    struct: parmed.Structure,
    residue_name: str,
) -> NDArray[np.float64]:
    """Get center-of-mass coordinates for each residue of given type."""
    coms = []
    for coords in iter_residue_coords(struct, residue_name):
        coms.append(coords.mean(axis=0))
    return np.array(coms)
```

#### 3.6 Leap Command Generation (for Amber prep)

```python
# Add to engines/amber.py:

def generate_ss_leap_commands(
    unit_name: str,
    ss_pairs: list[tuple[int, int]]
) -> list[str]:
    """
    Generate tleap commands for disulfide bonds.
    
    Equivalent to Biskit AmberParmBuilder.__fLines(ss_bond, ss)
    """
    commands = []
    for res1, res2 in ss_pairs:
        commands.append(f"bond {unit_name}.{res1}.SG {unit_name}.{res2}.SG")
    return commands
```

---

## Testing

Current: **569 tests passing** (updated 2026-03-22)

Each new module needs:
- Unit tests in `tests/test_<module>.py`
- Integration test if it interacts with files

Run: `cd ~/clawd/pymdmix2 && ~/.local/bin/pytest tests/ -q`

---

## File Locations

| Purpose | Location |
|---------|----------|
| New code | `~/clawd/pymdmix2/pymdmix/` |
| Tests | `~/clawd/pymdmix2/tests/` |
| This plan | `~/clawd/pymdmix2/docs/PLAN.md` |
| Legacy reference | `~/clawd/pyMDmix-port/pyMDMix/` |

---

## Quick Reference

```bash
# Run tests
cd ~/clawd/pymdmix2 && ~/.local/bin/pytest tests/ -q

# Check legacy code
grep -n "def something" ~/clawd/pyMDmix-port/pyMDMix/SomeFile.py

# Line counts
wc -l ~/clawd/pymdmix2/pymdmix/**/*.py | tail -5
```

---

## Progress Tracking

Update this section as work progresses:

| Date | Work Done |
|------|-----------|
| 2026-03-20 | Initial scaffolding, core modules |
| 2026-03-21 | Solvent library, CLI, tests |
| 2026-03-22 AM | Expanded replica, amber, grid; repo consolidation; Biskit audit |
| 2026-03-22 PM | **Biskit replacement COMPLETE**: 18 new functions in structure.py + amber.py |
| Next | MDSettings, ActionsManager, HotSpot clustering |

---

## Biskit Replacement Checklist

**Status: COMPLETE** ✅ (2026-03-22 13:55)
**Tests: 569 passing**

### structure.py Functions

**Masks:**
- [x] `get_protein_mask()`
- [x] `get_water_mask()`
- [x] `get_nucleic_mask()`
- [x] `get_ion_mask()`
- [x] `get_solute_mask()`
- [x] `get_solvent_mask()`
- [x] `get_heavy_atom_mask()`
- [x] `get_residue_mask()`
- [x] `get_atom_mask()`
- [x] `get_hydrogen_mask()` ✅ NEW
- [x] `get_backbone_mask()` ✅ NEW

**Selection:**
- [x] `select_atoms()`
- [x] `load_structure()`
- [x] `save_structure()`
- [x] `select_chains()` ✅ NEW

**Chains:**
- [x] `get_chain_ids()`
- [x] `get_chain_start_indices()` ✅ NEW

**Alignment (critical for trajectory processing):**
- [x] `align_structures()` ✅ NEW (pure numpy Kabsch, no MDAnalysis needed)
- [x] `align_to_reference()` ✅ NEW (convenience wrapper)
- [x] `compute_rmsd()` ✅ NEW

**Profiles:**
- [x] `atom_to_residue_values()` ✅ NEW
- [x] `residue_to_atom_values()` ✅ NEW

**Disulfides:**
- [x] `find_disulfides()`
- [x] `rename_cys_to_cyx()`

**Structure info:**
- [x] `count_residues()`
- [x] `center_structure()`
- [x] `get_residue_names()`

**Density/probe helpers:**
- [x] `detect_solvent_type()` ✅ NEW
- [x] `get_probe_coords()` ✅ NEW
- [x] `iter_residue_coords()` ✅ NEW
- [x] `get_residue_com_coords()` ✅ NEW

### engines/amber.py Functions

- [x] LeapSession (tleap wrapper)
- [x] AmberChecker
- [x] AmberWriter
- [x] `generate_ss_leap_commands()` ✅ NEW

### Summary

All 18 planned Biskit replacement functions have been implemented and tested.
The alignment uses pure numpy (Kabsch algorithm) - no MDAnalysis dependency required for this.
