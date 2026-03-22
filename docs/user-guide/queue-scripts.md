# Queue Scripts

pyMDMix generates queue submission scripts for HPC clusters. Customize these for your computing environment.

## Supported Queue Systems

- **SLURM** - Most common modern scheduler
- **PBS/Torque** - Traditional HPC scheduler
- **SGE** - Sun Grid Engine
- **LSF** - IBM Spectrum LSF

---

## Generating Queue Scripts

### During Replica Creation

Queue scripts are generated automatically when creating replicas:

```bash
pymdmix add replica -f replica.cfg
```

Creates `queue.sh` in each replica folder.

### Regenerate Scripts

```bash
# Regenerate for all replicas
pymdmix queue generate all

# For specific replicas
pymdmix queue generate byname -s MyProtein_ETA_1
```

---

## Queue Configuration

### Project-Level Settings

Create `queue.cfg` in your project:

```ini
[QUEUE]
# Queue system: slurm, pbs, sge, lsf
SYSTEM = slurm

# Partition/queue name
PARTITION = gpu

# Account for billing
ACCOUNT = myproject

# Wall time (HH:MM:SS)
WALLTIME = 24:00:00

# Number of nodes
NODES = 1

# CPUs per node
CPUS = 4

# GPUs per node
GPUS = 1

# Memory per node
MEMORY = 16G

# Email notifications
EMAIL = user@example.com
MAILTYPE = END,FAIL
```

### Apply to Replicas

```bash
pymdmix queue generate all --config queue.cfg
```

---

## Queue System Templates

### SLURM

```bash
#!/bin/bash
#SBATCH --job-name=MyProtein_ETA_1
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:1
#SBATCH --time=24:00:00
#SBATCH --mem=16G

module load amber/22
cd $SLURM_SUBMIT_DIR

# Run simulation
bash COMMANDS.sh
```

### PBS

```bash
#!/bin/bash
#PBS -N MyProtein_ETA_1
#PBS -q gpu
#PBS -l nodes=1:ppn=4:gpus=1
#PBS -l walltime=24:00:00
#PBS -l mem=16gb

module load amber/22
cd $PBS_O_WORKDIR

bash COMMANDS.sh
```

### SGE

```bash
#!/bin/bash
#$ -N MyProtein_ETA_1
#$ -q gpu.q
#$ -pe smp 4
#$ -l h_rt=24:00:00
#$ -l gpu=1

module load amber/22
cd $SGE_O_WORKDIR

bash COMMANDS.sh
```

---

## Custom Templates

### Create Template

Save as `my_template.sh`:

```bash
#!/bin/bash
#SBATCH --job-name={{JOB_NAME}}
#SBATCH --partition={{PARTITION}}
#SBATCH --time={{WALLTIME}}
#SBATCH --gres=gpu:{{GPUS}}

# Custom module loads
module purge
module load cuda/11.8
module load amber/22-gpu

# Custom environment
export CUDA_VISIBLE_DEVICES=0

cd {{WORKDIR}}
bash COMMANDS.sh
```

### Use Template

```bash
pymdmix queue generate all --template my_template.sh
```

### Template Variables

| Variable | Description |
|----------|-------------|
| `{{JOB_NAME}}` | Replica name |
| `{{WORKDIR}}` | Replica directory |
| `{{PARTITION}}` | Queue/partition name |
| `{{WALLTIME}}` | Wall time limit |
| `{{NODES}}` | Number of nodes |
| `{{CPUS}}` | CPUs per node |
| `{{GPUS}}` | GPUs per node |
| `{{MEMORY}}` | Memory allocation |

---

## Submitting Jobs

### Single Replica

```bash
cd MyProtein_ETA_1/
sbatch queue.sh
```

### All Replicas

```bash
# Generate submission script
pymdmix queue submit all --dry-run  # Preview
pymdmix queue submit all            # Submit all
```

### With Dependencies

```bash
# Submit equilibration first, then production
pymdmix queue submit all --stage eq
pymdmix queue submit all --stage prod --dependency afterok
```

---

## Monitoring Jobs

```bash
# SLURM
squeue -u $USER

# PBS
qstat -u $USER

# Check output
tail -f MyProtein_ETA_1/slurm-*.out
```

---

## Array Jobs

For many replicas, use array jobs:

```bash
# Generate array job script
pymdmix queue generate all --array

# Submit single script for all replicas
sbatch replica_array.sh
```

---

## Troubleshooting

### Job Fails Immediately

- Check module loads are correct
- Verify paths exist
- Check resource requests

### Out of Time

- Increase `WALLTIME`
- Use checkpointing/restarts
- Split into multiple jobs

### GPU Not Found

- Verify `--gres=gpu:1` or equivalent
- Check CUDA modules loaded
- Set `CUDA_VISIBLE_DEVICES`

---

## Next Steps

- [Creating Replicas](replicas.md) - Set up simulations
- [Analysis Overview](../analysis/overview.md) - After simulations complete
