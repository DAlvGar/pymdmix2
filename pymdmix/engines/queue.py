"""
Queue Script Generation
=======================

Generate job submission scripts for HPC clusters.

Supported systems:
- SLURM
- SGE (Sun Grid Engine)
- PBS/Torque
- Local (bash script)

Examples
--------
>>> from pymdmix.engines import QueueConfig, generate_queue_script
>>>
>>> config = QueueConfig(
...     system="slurm",
...     partition="gpu",
...     n_gpus=1,
...     time_hours=24,
... )
>>>
>>> script = generate_queue_script(
...     config=config,
...     job_name="md_production",
...     commands=["pmemd.cuda -O -i prod.in ..."],
... )
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field, fields
from pathlib import Path
from string import Template

log = logging.getLogger(__name__)


@dataclass
class QueueConfig:
    """
    Queue system configuration.

    Attributes
    ----------
    system : str
        Queue system: slurm, sge, pbs, local
    partition : str
        Partition/queue name
    n_nodes : int
        Number of nodes
    n_cpus : int
        CPUs per node
    n_gpus : int
        GPUs per node
    memory_gb : int
        Memory per node in GB
    time_hours : int
        Wall time in hours
    account : str
        Account/project for billing
    email : str
        Email for notifications
    """

    system: str = "slurm"
    partition: str = "gpu"
    n_nodes: int = 1
    n_cpus: int = 4
    n_gpus: int = 1
    memory_gb: int = 16
    time_hours: int = 24
    account: str | None = None
    email: str | None = None
    extra_directives: list[str] = field(default_factory=list)
    modules: list[str] = field(default_factory=list)
    environment: dict[str, str] = field(default_factory=dict)

    @classmethod
    def from_file(cls, path: Path | str) -> QueueConfig:
        """
        Load QueueConfig from a TOML or JSON file.

        Parameters
        ----------
        path : Path | str
            Config file path (.toml or .json)

        Returns
        -------
        QueueConfig
            Loaded configuration

        Examples
        --------
        >>> cfg = QueueConfig.from_file("queue.toml")
        """
        import json

        path = Path(path)
        if not path.exists():
            raise FileNotFoundError(f"Queue config file not found: {path}")

        if path.suffix == ".toml":
            try:
                import tomllib
            except ImportError:
                import tomli as tomllib  # type: ignore[no-redef]
            with open(path, "rb") as f:
                data = tomllib.load(f)
        elif path.suffix == ".json":
            with open(path) as f:
                data = json.load(f)
        else:
            raise ValueError(
                f"Unsupported queue config format: {path.suffix!r}. Use .toml or .json"
            )

        valid = {f.name for f in fields(cls)}
        filtered = {k: v for k, v in data.items() if k in valid}
        return cls(**filtered)


# SLURM template
SLURM_TEMPLATE = Template("""#!/bin/bash
#SBATCH --job-name=${job_name}
#SBATCH --partition=${partition}
#SBATCH --nodes=${n_nodes}
#SBATCH --ntasks-per-node=${n_cpus}
#SBATCH --gres=gpu:${n_gpus}
#SBATCH --mem=${memory_gb}G
#SBATCH --time=${time_hours}:00:00
#SBATCH --output=${job_name}_%j.out
#SBATCH --error=${job_name}_%j.err
${account_line}${email_line}${extra_directives}

# Environment setup
${modules}
${environment}

# Job information
echo "Job started on $$(hostname) at $$(date)"
echo "Working directory: $$(pwd)"

# Commands
${commands}

echo "Job finished at $$(date)"
""")

# SGE template
SGE_TEMPLATE = Template("""#!/bin/bash
#$$$$ -N ${job_name}
#$$$$ -q ${partition}
#$$$$ -pe smp ${n_cpus}
#$$$$ -l h_rt=${time_hours}:00:00
#$$$$ -l h_vmem=${memory_gb}G
#$$$$ -cwd
#$$$$ -o ${job_name}.o$$$$JOB_ID
#$$$$ -e ${job_name}.e$$$$JOB_ID
${extra_directives}

# Environment setup
${modules}
${environment}

# Commands
${commands}
""")

# PBS template
PBS_TEMPLATE = Template("""#!/bin/bash
#PBS -N ${job_name}
#PBS -q ${partition}
#PBS -l nodes=${n_nodes}:ppn=${n_cpus}:gpus=${n_gpus}
#PBS -l walltime=${time_hours}:00:00
#PBS -l mem=${memory_gb}gb
#PBS -o ${job_name}.o
#PBS -e ${job_name}.e
${extra_directives}

cd $$$$PBS_O_WORKDIR

# Environment setup
${modules}
${environment}

# Commands
${commands}
""")

# Local template (simple bash)
LOCAL_TEMPLATE = Template("""#!/bin/bash
# Local execution script for ${job_name}

# Environment setup
${modules}
${environment}

echo "Starting ${job_name} at $$(date)"

# Commands
${commands}

echo "Finished ${job_name} at $$(date)"
""")


def generate_queue_script(
    config: QueueConfig,
    job_name: str,
    commands: list[str],
    work_dir: Path | None = None,
) -> str:
    """
    Generate a queue submission script.

    Parameters
    ----------
    config : QueueConfig
        Queue configuration
    job_name : str
        Job name
    commands : list[str]
        Commands to execute
    work_dir : Path
        Working directory

    Returns
    -------
    str
        Script content
    """
    # Select template
    templates = {
        "slurm": SLURM_TEMPLATE,
        "sge": SGE_TEMPLATE,
        "pbs": PBS_TEMPLATE,
        "local": LOCAL_TEMPLATE,
    }

    template = templates.get(config.system.lower())
    if template is None:
        raise ValueError(f"Unknown queue system: {config.system}")

    # Format modules
    if config.modules:
        modules = "\n".join(f"module load {m}" for m in config.modules)
    else:
        modules = "# No modules to load"

    # Format environment
    if config.environment:
        env_lines = [f"export {k}={v}" for k, v in config.environment.items()]
        environment = "\n".join(env_lines)
    else:
        environment = ""

    # Format extra directives
    extra = ""
    if config.extra_directives:
        if config.system == "slurm":
            extra = "\n".join(f"#SBATCH {d}" for d in config.extra_directives)
        elif config.system == "sge":
            extra = "\n".join(f"#$ {d}" for d in config.extra_directives)
        elif config.system == "pbs":
            extra = "\n".join(f"#PBS {d}" for d in config.extra_directives)

    # Account line (SLURM specific)
    account_line = ""
    if config.account and config.system == "slurm":
        account_line = f"#SBATCH --account={config.account}\n"

    # Email line (SLURM specific)
    email_line = ""
    if config.email and config.system == "slurm":
        email_line = f"#SBATCH --mail-user={config.email}\n#SBATCH --mail-type=END,FAIL\n"

    # Format commands
    commands_str = "\n".join(commands)

    # Substitute
    script = template.substitute(
        job_name=job_name,
        partition=config.partition,
        n_nodes=config.n_nodes,
        n_cpus=config.n_cpus,
        n_gpus=config.n_gpus,
        memory_gb=config.memory_gb,
        time_hours=config.time_hours,
        modules=modules,
        environment=environment,
        commands=commands_str,
        extra_directives=extra,
        account_line=account_line,
        email_line=email_line,
    )

    return script


def generate_mdmix_production_script(
    config: QueueConfig,
    job_name: str,
    topology: Path,
    coordinates: Path,
    n_runs: int = 10,
    steps_per_run: int = 500000,
    amber_exe: str = "pmemd.cuda",
) -> str:
    """
    Generate a complete MDMix production script.

    Runs multiple sequential production runs with automatic
    restart handling.

    Parameters
    ----------
    config : QueueConfig
        Queue configuration
    job_name : str
        Job name
    topology : Path
        Topology file
    coordinates : Path
        Initial coordinates
    n_runs : int
        Number of production runs
    steps_per_run : int
        Steps per run
    amber_exe : str
        Amber executable

    Returns
    -------
    str
        Script content
    """
    from pymdmix.engines.amber import AmberEngine

    engine = AmberEngine(exe=amber_exe)

    # Generate production input
    mdin_content = engine.production_input(nsteps=steps_per_run)

    # Build commands
    commands = [
        "# Write production input",
        "cat > prod.in << 'EOF'",
        mdin_content.strip(),
        "EOF",
        "",
        "# Production runs",
        f"TOPOLOGY={topology}",
        f"COORDS={coordinates}",
        "",
        f"for i in $(seq 1 {n_runs}); do",
        '    echo "Starting run $i"',
        "    if [ $i -eq 1 ]; then",
        "        INPUT=$COORDS",
        "    else",
        "        PREV=$((i-1))",
        "        INPUT=md${PREV}.rst7",
        "    fi",
        f"    {amber_exe} -O \\",
        "        -i prod.in \\",
        "        -o md${i}.out \\",
        "        -p $TOPOLOGY \\",
        "        -c $INPUT \\",
        "        -r md${i}.rst7 \\",
        "        -x md${i}.nc \\",
        "        -ref $INPUT \\",
        "        -inf md${i}.mdinfo",
        '    echo "Finished run $i"',
        "done",
    ]

    return generate_queue_script(config, job_name, commands)
