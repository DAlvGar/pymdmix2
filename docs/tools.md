# Additional Tools

pyMDMix includes utility tools for grid manipulation, format conversion, and analysis helpers.

## Grid Tools

### Grid Information

```bash
# Display grid properties
pymdmix tools grid-info grids/ETA_CT_DG.dx
```

Output:
```
File: ETA_CT_DG.dx
Dimensions: 100 x 100 x 100
Spacing: 0.5 Å
Origin: (-25.0, -25.0, -25.0)
Min value: -2.31
Max value: 5.67
```

### Grid Math

```bash
# Add grids
pymdmix tools grid-math add grid1.dx grid2.dx -o sum.dx

# Average grids
pymdmix tools grid-math average grid1.dx grid2.dx grid3.dx -o avg.dx

# Minimum (most favorable)
pymdmix tools grid-math min grid1.dx grid2.dx -o min.dx

# Scale grid values
pymdmix tools grid-math scale grid.dx -f 0.5 -o scaled.dx
```

### Grid Conversion

```bash
# Convert DX to CCP4/MRC
pymdmix tools convert grids/ETA_CT_DG.dx -o ETA_CT_DG.mrc

# Convert to Gaussian cube
pymdmix tools convert grids/ETA_CT_DG.dx -o ETA_CT_DG.cube
```

### Grid Resampling

```bash
# Change grid resolution
pymdmix tools resample grid.dx --spacing 0.25 -o fine_grid.dx

# Crop to region
pymdmix tools crop grid.dx --center 10,20,30 --size 20 -o cropped.dx
```

---

## Hotspot Tools

### Combine Hotspots

Merge hotspots from multiple replicas/solvents:

```bash
pymdmix tools combine-hotspots \
    replica1/hotspots/hotspots.pdb \
    replica2/hotspots/hotspots.pdb \
    -o combined.pdb \
    --distance 3.0  # Merge within 3 Å
```

### Filter Hotspots

```bash
# Keep only strong hotspots
pymdmix tools filter-hotspots hotspots.pdb \
    --min-energy -1.5 \
    -o strong_hotspots.pdb

# Keep only specific probe types
pymdmix tools filter-hotspots hotspots.pdb \
    --probes HYD,DON \
    -o filtered.pdb
```

### Hotspot Clustering

```bash
# Cluster nearby hotspots
pymdmix tools cluster-hotspots hotspots.pdb \
    --cutoff 5.0 \
    -o clustered.pdb
```

---

## Trajectory Tools

### Extract Frames

```bash
# Extract specific frames
pymdmix tools extract-frames trajectory.nc \
    --frames 1,100,200 \
    -o frames.pdb

# Extract every Nth frame
pymdmix tools extract-frames trajectory.nc \
    --stride 10 \
    -o subsampled.nc
```

### Concatenate Trajectories

```bash
pymdmix tools concat-traj md1.nc md2.nc md3.nc -o combined.nc
```

### Convert Formats

```bash
# NetCDF to DCD
pymdmix tools convert-traj trajectory.nc -o trajectory.dcd

# DCD to NetCDF
pymdmix tools convert-traj trajectory.dcd \
    --topology system.prmtop \
    -o trajectory.nc
```

---

## Structure Tools

### Generate Reference

```bash
# Create reference PDB from topology
pymdmix tools make-ref system.prmtop system.inpcrd -o reference.pdb
```

### Extract Protein

```bash
# Extract protein only from solvated system
pymdmix tools extract-solute system.prmtop system.inpcrd \
    -o protein.pdb
```

---

## Density to Energy Tools

### Batch Conversion

```bash
# Convert all density grids in a folder
pymdmix tools density-to-energy grids/ \
    --temp 300 \
    --suffix _DG
```

### Custom Normalization

```bash
# Use custom bulk density
pymdmix tools density-to-energy grid_density.dx \
    --bulk-density 0.033 \
    -o grid_DG.dx
```

---

## Solvent Database Tools

### List Solvents

```bash
pymdmix info solvents
pymdmix info solvents --detailed
```

### Export Solvent

```bash
# Export solvent definition to JSON
pymdmix tools export-solvent ETA -o eta_solvent.json
```

### Validate Solvent

```bash
# Check solvent configuration
pymdmix tools validate-solvent my_solvent.cfg
```

---

## Python API

```python
from pymdmix.core import Grid
from pymdmix.io import grids

# Load and manipulate grids
grid1 = Grid.from_dx("grid1.dx")
grid2 = Grid.from_dx("grid2.dx")

# Average
avg = Grid.average([grid1, grid2])
avg.write_dx("averaged.dx")

# Minimum
combined = Grid.minimum([grid1, grid2])

# Convert formats
grids.dx_to_mrc("input.dx", "output.mrc")
grids.dx_to_cube("input.dx", "output.cube", structure="protein.pdb")
```

---

## Scripting Examples

### Batch Process All Replicas

```bash
#!/bin/bash
for replica in */; do
    if [ -d "$replica/grids" ]; then
        echo "Processing $replica..."
        pymdmix tools grid-info "$replica/grids/"*.dx
    fi
done
```

### Compare Replica Energies

```python
#!/usr/bin/env python
from pymdmix.core import Grid
from pathlib import Path

replicas = list(Path(".").glob("*_ETA_*/grids/ETA_CT_DG.dx"))

for rep in replicas:
    grid = Grid.from_dx(rep)
    print(f"{rep.parent.parent.name}: min={grid.data.min():.2f} kcal/mol")
```

---

## Getting Help

```bash
# List all tools
pymdmix tools --help

# Help for specific tool
pymdmix tools grid-math --help
pymdmix tools combine-hotspots --help
```
