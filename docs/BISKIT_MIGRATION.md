# Biskit Migration Assessment

## TL;DR: Migration is NOT trivial

The codebase has **deep entanglement** with Biskit, not just simple imports.

---

## Critical Dependencies

### 1. Class Inheritance (HARD)

```python
# PDB.py - Direct subclass
class SolvatedPDB(bi.PDBModel):
    """Custom PDB handling for solvated systems"""
    
# AutoPrepare.py - Inherits Amber tools  
class AmberPDBCleaner(bi.AmberParmBuilder):
    """Structure cleaning with ACE/NME capping"""
```

**Problem**: These classes inherit behavior, not just data structures.

### 2. PDBModel Methods Used (33 unique calls)

| Method | Count | Description | parmed Equivalent |
|--------|-------|-------------|-------------------|
| `maskFrom()` | 11 | Boolean mask from attribute | Custom function needed |
| `compress()` | 10 | Select atoms by mask | `struct[mask]` or loop |
| `writePdb()` | 9 | Write PDB file | `struct.save()` ✅ |
| `concat()` | 6 | Concatenate models | `struct + struct2` or manual |
| `take()` | 4 | Take subset by indices | Index selection |
| `maskProtein()` | 4 | Mask protein atoms | `get_protein_mask()` ✅ |
| `maskHeavy()` | 3 | Mask heavy atoms | Custom function needed |
| `maskBB()` | 3 | Mask backbone | Custom function needed |
| `maskNA()` | 2 | Mask nucleic acids | `get_nucleic_mask()` ✅ |
| `maskH()` | 2 | Mask hydrogens | Custom function needed |
| `atom2resProfile()` | - | Atom→residue aggregation | Custom function needed |
| `res2atomProfile()` | - | Residue→atom expansion | Custom function needed |
| `lenResidues()` | - | Count residues | `len(struct.residues)` ✅ |
| `takeChains()` | - | Extract chains | Custom selection |
| `resModels()` | - | Get residue models | Custom iteration |
| `magicFit()` | - | Superposition | MDAnalysis/scipy |
| `mergeChains()` | - | Merge chains | Custom function |
| `remove()` | - | Remove atoms | `struct.strip()` or filter |

### 3. Dictionary-like Access

```python
# Biskit PDBModel provides:
pdb['residue_name']      # Array of all residue names
pdb['serial_number']     # Array of atom numbers  
pdb.xyz                  # Coordinates array

# parmed Structure provides:
[a.residue.name for a in struct.atoms]  # List, needs conversion
struct.coordinates                       # Works similar ✅
```

### 4. AmberParmBuilder Dependencies

The `AmberPDBCleaner` class uses Biskit's Amber preparation tools:
- ACE/NME capping (complex chain manipulation)
- Chain extraction and merging
- Atom renumbering
- Terminal handling

---

## Migration Strategy Options

### Option A: Wrapper Class (RECOMMENDED)

Create a `SolvatedStructure` class that wraps parmed.Structure and provides Biskit-like API:

```python
class SolvatedStructure:
    """parmed.Structure with Biskit-like convenience methods."""
    
    def __init__(self, source: str | Path | parmed.Structure):
        if isinstance(source, parmed.Structure):
            self._struct = source
        else:
            self._struct = parmed.load_file(str(source))
        self._cache = {}  # Cache expensive operations
    
    @property
    def xyz(self) -> np.ndarray:
        return self._struct.coordinates
    
    def __getitem__(self, key: str) -> np.ndarray:
        """Biskit-style attribute access."""
        if key == 'residue_name':
            return np.array([a.residue.name for a in self._struct.atoms])
        elif key == 'serial_number':
            return np.array([a.idx + 1 for a in self._struct.atoms])
        # ... etc
    
    def mask_from(self, attr: str, values: list) -> np.ndarray:
        """Create boolean mask where attribute matches values."""
        data = self[attr]
        return np.isin(data, values)
    
    def compress(self, mask: np.ndarray) -> SolvatedStructure:
        """Return new structure with only masked atoms."""
        # parmed doesn't have direct mask support, need to build selection
        indices = np.where(mask)[0]
        new_struct = self._struct[indices]  # If supported
        return SolvatedStructure(new_struct)
    
    def mask_protein(self) -> np.ndarray:
        return get_protein_mask(self._struct)
    
    def mask_heavy(self) -> np.ndarray:
        return np.array([a.element != 1 for a in self._struct.atoms])
    
    def mask_backbone(self) -> np.ndarray:
        bb_names = {'CA', 'C', 'N', 'O'}
        return np.array([a.name in bb_names for a in self._struct.atoms])
    
    def write_pdb(self, path: str | Path) -> None:
        self._struct.save(str(path), overwrite=True)
    
    # ... implement remaining methods
```

**Effort**: ~500 lines, 2-3 days

### Option B: Pure Functions (CURRENT APPROACH)

Keep `parmed.Structure` as-is, add helper functions:

```python
# Already in core/structure.py
def get_protein_mask(struct): ...
def get_water_mask(struct): ...
def get_residue_mask(struct, resnames): ...
```

**Problem**: Doesn't handle inherited classes or complex operations.

### Option C: MDAnalysis Universe

MDAnalysis has richer selection language:

```python
import MDAnalysis as mda
u = mda.Universe("protein.pdb")
protein = u.select_atoms("protein")
heavy = u.select_atoms("not name H*")
backbone = u.select_atoms("backbone")
```

**Problem**: Different data model, would require larger refactor.

---

## Implementation Priority

### Phase 1: Core Masks (✅ DONE)
- `get_protein_mask()`
- `get_water_mask()`
- `get_nucleic_mask()`
- `get_solute_mask()`
- `get_solvent_mask()`
- `get_residue_mask()`

### Phase 2: SolvatedStructure Wrapper (TODO)
Replace `SolvatedPDB` with:
```
pymdmix2/pymdmix/core/solvated_structure.py
```

Key methods needed:
- [ ] `__getitem__` for attribute arrays
- [ ] `mask_from()` 
- [ ] `compress()` / `take()`
- [ ] `concat()`
- [ ] `mask_heavy()`, `mask_backbone()`, `mask_hydrogen()`
- [ ] `atom2res_profile()`, `res2atom_profile()`
- [ ] Coordinate access (`xyz` property)
- [ ] Probe coordinate extraction

### Phase 3: ACE/NME Capping (TODO)
Replace `AmberPDBCleaner` with:
- Use parmed or pdb4amber for capping
- Or implement chain manipulation directly

### Phase 4: Remaining Biskit.tools (TODO)
Functions in `tools.py` that came from Biskit:
- `projectRoot()` → `project_root()` ✅
- `dataRoot()` → `data_root()` ✅
- Pickle utilities → standard pickle ✅
- File utilities → pathlib ✅

---

## Risk Assessment

| Component | Risk | Mitigation |
|-----------|------|------------|
| SolvatedPDB | HIGH | Create SolvatedStructure wrapper |
| AmberPDBCleaner | MEDIUM | Use pdb4amber or parmed's tools |
| PDBModel masks | LOW | Already have helper functions |
| Coordinate access | LOW | parmed.Structure.coordinates works |
| Chain operations | MEDIUM | Need custom implementation |
| Superposition | LOW | Use MDAnalysis or scipy |

---

## Recommendation

1. **Create `SolvatedStructure` class** wrapping parmed - this gives us the Biskit-like API without Biskit
2. **Don't port AmberPDBCleaner** - use `pdb4amber` CLI or parmed's `AmberOFFLibrary` instead
3. **Test incrementally** - port one use case at a time with tests

This is **2-3 days of focused work**, not a simple find-replace.
