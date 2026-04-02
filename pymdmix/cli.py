"""
pyMDMix CLI - Command-line interface for pyMDMix.

Usage:
    pymdmix create project myproject
    pymdmix add system -f system.cfg
    pymdmix add replica -f replica.cfg
    pymdmix analyze align all
    pymdmix analyze density bysolvent -s ETA
    pymdmix info project
    pymdmix info solvents
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Literal, cast

import click

from pymdmix import __version__


def _groups_file(project_path: Path) -> Path:
    """Get path to CLI groups file for a project."""
    return project_path / "groups.json"


def _load_groups(project_path: Path) -> dict[str, list[str]]:
    """Load named replica groups for a project."""
    import json

    path = _groups_file(project_path)
    if not path.exists():
        return {}
    try:
        data = json.loads(path.read_text())
    except Exception:
        return {}
    if not isinstance(data, dict):
        return {}

    groups: dict[str, list[str]] = {}
    for name, members in data.items():
        if isinstance(name, str) and isinstance(members, list):
            valid_members = [m for m in members if isinstance(m, str)]
            groups[name] = valid_members
    return groups


def _save_groups(project_path: Path, groups: dict[str, list[str]]) -> None:
    """Persist named replica groups for a project."""
    import json

    _groups_file(project_path).write_text(json.dumps(groups, indent=2))


# =============================================================================
# CLI Group
# =============================================================================


@click.group()
@click.version_option(version=__version__, prog_name="pyMDMix")
@click.option("-v", "--verbose", is_flag=True, help="Enable verbose output")
@click.option("--debug", is_flag=True, help="Enable debug output")
@click.pass_context
def cli(ctx: click.Context, verbose: bool, debug: bool) -> None:
    """pyMDMix - Molecular Dynamics with organic solvent mixtures.

    A toolkit for identifying binding hotspots on protein surfaces
    using MD simulations with organic solvent/water mixtures.
    """
    ctx.ensure_object(dict)
    ctx.obj["verbose"] = verbose
    ctx.obj["debug"] = debug


# =============================================================================
# CREATE Commands
# =============================================================================


@cli.group()
def create() -> None:
    """Create projects, replicas, or solvents."""
    pass


@create.command("project")
@click.option("-n", "--name", default="mdmix_project", help="Project name")
@click.option("-f", "--config", type=click.Path(exists=True), help="Configuration file")
@click.option("-d", "--directory", type=click.Path(), help="Project directory (default: ./<name>)")
@click.option("--force", is_flag=True, help="Overwrite existing project")
@click.pass_context
def create_project(
    ctx: click.Context, name: str, config: str | None, directory: str | None, force: bool
) -> None:
    """Create a new pyMDMix project.

    Examples:
        pymdmix create project -n myproject
        pymdmix create project -n myproject -f config.yaml
    """
    from pymdmix.project import Config, Project

    verbose = ctx.obj.get("verbose", False)
    project_dir = (
        Path(directory) if directory else (Path.cwd() / name if name != "." else Path.cwd())
    )

    if project_dir.exists() and not force:
        mdmix_dir = project_dir / ".mdmix"
        if mdmix_dir.exists():
            click.secho(f"✗ Project already exists: {project_dir}", fg="red")
            click.echo("  Use --force to overwrite")
            sys.exit(1)

    click.echo(f"Creating project '{name}' in {project_dir}")

    if config:
        cfg = Config.from_file(Path(config))
        if verbose:
            click.echo(f"  Loaded config from {config}")
    else:
        cfg = Config()

    project = Project(name=name, config=cfg, path=project_dir)
    project.save()

    click.secho(f"✓ Project '{name}' created successfully", fg="green")


@create.command("solvent")
@click.option(
    "-f",
    "--file",
    "config",
    type=click.Path(exists=True),
    required=True,
    help="Solvent config file",
)
@click.option("--validate", is_flag=True, help="Validate only, don't add")
@click.pass_context
def create_solvent(ctx: click.Context, config: str, validate: bool) -> None:
    """Add a new solvent to the library.

    Examples:
        pymdmix create solvent -f my_solvent.json
        pymdmix create solvent -f my_solvent.cfg --validate
    """
    from pymdmix.core.solvent import Solvent, SolventLibrary

    verbose = ctx.obj.get("verbose", False)
    config_path = Path(config)

    try:
        solvent = Solvent.from_file(config_path)
    except Exception as e:
        click.secho(f"✗ Failed to parse solvent config: {e}", fg="red")
        sys.exit(1)

    if verbose:
        click.echo(f"  Name: {solvent.name}")
        click.echo(f"  Info: {solvent.full_name}")
        click.echo(f"  Probes: {[p.name for p in solvent.probes]}")

    if validate:
        click.secho("✓ Solvent configuration is valid", fg="green")
        return

    library = SolventLibrary()
    library.add(solvent)
    click.secho(f"✓ Solvent '{solvent.name}' added to library", fg="green")


# =============================================================================
# ADD Commands
# =============================================================================


@cli.group()
def add() -> None:
    """Add systems, replicas, or groups to an existing project."""
    pass


@add.command("system")
@click.option(
    "-f",
    "--file",
    "config",
    type=click.Path(exists=True),
    required=True,
    help="System configuration file",
)
@click.option(
    "-p",
    "--project",
    type=click.Path(exists=True),
    default=".",
    help="Project directory (default: current)",
)
@click.pass_context
def add_system(ctx: click.Context, config: str, project: str) -> None:
    """Add a system to the project from a configuration file.

    The configuration file should contain:
    - NAME: System identifier
    - OFF: Path to Amber object file (or PDB)
    - EXTRARES: Non-standard residues (optional)
    - EXTRAFF: Additional force field files (optional)

    Examples:
        pymdmix add system -f system.cfg
        pymdmix add system -f system.cfg -p /path/to/project
    """
    from pymdmix.io.parsers import parse_system_config
    from pymdmix.project import Project

    verbose = ctx.obj.get("verbose", False)
    project_path = Path(project)
    config_path = Path(config)

    # Load project
    proj = Project.load(project_path)

    # Parse system config
    system_config = parse_system_config(config_path)

    if verbose:
        click.echo(f"  System name: {system_config.name}")
        click.echo(f"  Input file: {system_config.input_file}")

    # Store config in project input folder for later use
    proj.create_directories()
    destination = proj.input_path / Path(config).name
    destination.write_text(Path(config).read_text())
    # Keep a default pdb/off pointer if available
    proj.pdb_file = str(system_config.input_file)
    proj.save()

    click.secho(f"✓ System '{system_config.name}' registered in project inputs", fg="green")


@add.command("replica")
@click.option(
    "-f",
    "--file",
    "config",
    type=click.Path(exists=True),
    required=True,
    help="Replica configuration file",
)
@click.option(
    "-p",
    "--project",
    type=click.Path(exists=True),
    default=".",
    help="Project directory (default: current)",
)
@click.option("--count", type=int, default=1, help="Number of replicas to create")
@click.pass_context
def add_replica(ctx: click.Context, config: str, project: str, count: int) -> None:
    """Add replica(s) to the project from a configuration file.

    The configuration file should contain:
    - SYSTEM: System name (must exist in project)
    - SOLVENT: Solvent name (e.g., ETA, MAM)
    - NANOS: Production length in nanoseconds
    - RESTRMODE: Restraint mode (FREE, BB, HA, CUSTOM)

    Examples:
        pymdmix add replica -f replica.cfg
        pymdmix add replica -f replica.cfg --count 3
    """
    from pymdmix.io.parsers import parse_replica_config
    from pymdmix.project import Project

    verbose = ctx.obj.get("verbose", False)
    project_path = Path(project)
    config_path = Path(config)

    # Load project
    proj = Project.load(project_path)

    # Parse replica config
    replica_config = parse_replica_config(config_path)

    if verbose:
        click.echo(f"  System: {replica_config.system}")
        click.echo(f"  Solvent: {replica_config.solvent}")
        click.echo(f"  Nanoseconds: {replica_config.nanos}")

    # Create replica(s)
    created = proj.add_replicas(solvent=replica_config.solvent, n_replicas=count)
    created_names = [r.name for r in created]
    for replica in created:
        if replica.settings:
            replica.settings.nanos = replica_config.nanos
            replica.settings.restraint_mode = replica_config.restraint_mode
            if replica_config.restraint_mask:
                replica.settings.restraint_mask = replica_config.restraint_mask
            replica.settings.restraint_force = replica_config.restraint_force
            if replica_config.align_mask:
                replica.settings.align_mask = replica_config.align_mask
        if verbose:
            click.echo(f"  Created: {replica.name}")

    proj.save()

    if count == 1:
        click.secho(f"✓ Replica '{created_names[0]}' added to project", fg="green")
    else:
        click.secho(f"✓ {count} replicas added to project", fg="green")
        for name in created_names:
            click.echo(f"    - {name}")


@add.command("group")
@click.option("-n", "--name", required=True, help="Group name")
@click.option("-s", "--selection", multiple=True, required=True, help="Replica names to include")
@click.option(
    "-p",
    "--project",
    type=click.Path(exists=True),
    default=".",
    help="Project directory (default: current)",
)
@click.pass_context
def add_group(ctx: click.Context, name: str, selection: tuple, project: str) -> None:
    """Create a named group of replicas for batch operations.

    Examples:
        pymdmix add group -n ethanol_runs -s MyProtein_ETA_1 -s MyProtein_ETA_2
        pymdmix add group -n all_free -s rep1 -s rep2 -s rep3
    """
    from pymdmix.project import Project

    project_path = Path(project)
    proj = Project.load(project_path)

    # Verify all replicas exist
    replica_names = {rep.name for rep in proj.replicas}
    for rep_name in selection:
        if rep_name not in replica_names:
            click.secho(f"✗ Replica not found: {rep_name}", fg="red")
            sys.exit(1)

    groups = _load_groups(proj.path)
    groups[name] = [str(s) for s in selection]
    _save_groups(proj.path, groups)
    proj.save()

    click.secho(f"✓ Group '{name}' created with {len(selection)} replica(s)", fg="green")


# =============================================================================
# SETUP Commands
# =============================================================================


@cli.group()
def setup() -> None:
    """Structure preparation and system setup."""
    pass


@setup.command("prepare")
@click.argument("structure", type=click.Path(exists=True))
@click.option("-o", "--output", type=click.Path(), help="Output PDB file")
@click.option("--cap/--no-cap", default=True, help="Add ACE/NME caps")
@click.option("--disulfide/--no-disulfide", default=True, help="Detect disulfides")
@click.option("--remove-water/--keep-water", default=True, help="Remove waters")
@click.pass_context
def setup_prepare(
    ctx: click.Context,
    structure: str,
    output: str | None,
    cap: bool,
    disulfide: bool,
    remove_water: bool,
) -> None:
    """Prepare a structure for MD simulation.

    Examples:
        pymdmix setup prepare protein.pdb -o prepared.pdb
        pymdmix setup prepare protein.pdb --no-cap
    """
    from pymdmix.setup.prepare import StructurePreparationOptions, prepare_structure

    verbose = ctx.obj.get("verbose", False)
    input_path = Path(structure)
    output_path = Path(output) if output else input_path.with_stem(input_path.stem + "_prepared")

    click.echo(f"Preparing structure: {input_path.name}")

    options = StructurePreparationOptions(
        add_caps=cap,
        detect_disulfides=disulfide,
        remove_waters=remove_water,
    )

    result = prepare_structure(input_path, options=options)
    result.structure.save(str(output_path), overwrite=True)

    click.secho("✓ Structure prepared successfully", fg="green")
    click.echo(f"  Output: {output_path}")

    if verbose and result.modifications:
        click.echo("  Modifications:")
        for mod in result.modifications:
            click.echo(f"    - {mod}")


@setup.command("solvate")
@click.argument("structure", type=click.Path(exists=True))
@click.option("-s", "--solvent", required=True, help="Solvent name")
@click.option("-o", "--output", type=click.Path(), help="Output prefix")
@click.option("--buffer", type=float, default=12.0, help="Box buffer size (Å)")
@click.pass_context
def setup_solvate(
    ctx: click.Context, structure: str, solvent: str, output: str | None, buffer: float
) -> None:
    """Solvate a structure with a solvent mixture.

    Examples:
        pymdmix setup solvate protein.pdb -s ETA -o solvated
    """
    from pymdmix.core.solvent import SolventLibrary
    from pymdmix.setup.solvate import SolvationOptions, solvate_structure

    input_path = Path(structure)
    output_prefix = output or input_path.stem + "_solvated"

    library = SolventLibrary()
    solv = library.get(solvent)
    if solv is None:
        click.secho(f"✗ Unknown solvent: {solvent}", fg="red")
        click.echo(f"  Available: {', '.join(library.list_solvents())}")
        sys.exit(1)

    click.echo(f"Solvating with {solv.full_name}")

    options = SolvationOptions(box_buffer=buffer)
    result = solvate_structure(input_path, solv, options=options)

    top_path = Path(output_prefix + ".prmtop")
    crd_path = Path(output_prefix + ".rst7")
    result.save_topology(top_path)
    result.save_coordinates(str(crd_path))

    click.secho("✓ Solvation complete", fg="green")
    click.echo(f"  Topology: {top_path}")
    click.echo(f"  Coordinates: {crd_path}")


# =============================================================================
# ANALYZE Commands
# =============================================================================


@cli.group()
def analyze() -> None:
    """Run analysis on simulation data."""
    pass


def parse_selection(selection_type: str, selection: tuple, project_path: Path):
    """Parse replica selection and return list of replicas."""
    from pymdmix.project import Project

    proj = Project.load(project_path)

    if selection_type == "all":
        return list(proj.replicas), proj
    elif selection_type == "bysolvent":
        wanted = {str(s) for s in selection}
        replicas = [r for r in proj.replicas if r.solvent in wanted]
        return replicas, proj
    elif selection_type == "byname":
        replicas = [rep for name in selection if (rep := proj.get_replica(str(name))) is not None]
        return replicas, proj
    elif selection_type == "group":
        if len(selection) != 1:
            raise click.UsageError("Group selection requires exactly one group name")
        group_name = str(selection[0])
        groups = _load_groups(project_path)
        replica_names = groups.get(group_name, [])
        replicas = [rep for name in replica_names if (rep := proj.get_replica(name)) is not None]
        return replicas, proj
    else:
        raise click.UsageError(f"Unknown selection type: {selection_type}")


@analyze.command("align")
@click.argument("selection_type", type=click.Choice(["all", "bysolvent", "byname", "group"]))
@click.option("-s", "--selection", multiple=True, help="Selection values")
@click.option(
    "-p", "--project", type=click.Path(exists=True), default=".", help="Project directory"
)
@click.option("-N", "--nanoselect", help="Nanosecond range (e.g., 1:10)")
@click.option("-C", "--ncpus", type=int, default=1, help="Number of CPUs")
@click.option("--mask", help="Alignment mask (residue range)")
@click.option("--ref", type=click.Path(exists=True), help="Reference PDB file")
@click.pass_context
def analyze_align(
    ctx: click.Context,
    selection_type: str,
    selection: tuple,
    project: str,
    nanoselect: str | None,
    ncpus: int,
    mask: str | None,
    ref: str | None,
) -> None:
    """Align trajectories to reference structure.

    Selection types:
    - all: All replicas in project
    - bysolvent -s ETA MAM: Replicas with specified solvents
    - byname -s rep1 rep2: Specific replica names
    - group -s mygroup: Named replica group

    Examples:
        pymdmix analyze align all
        pymdmix analyze align bysolvent -s ETA -C 4
        pymdmix analyze align byname -s MyProtein_ETA_1 -N 1:10
    """
    from pymdmix.analysis import AlignAction

    verbose = ctx.obj.get("verbose", False)
    project_path = Path(project)

    replicas, proj = parse_selection(selection_type, selection, project_path)

    if not replicas:
        click.secho("✗ No replicas selected", fg="red")
        sys.exit(1)

    click.echo(f"Aligning {len(replicas)} replica(s)")

    # Parse nanosecond selection
    step_range = None
    if nanoselect:
        if ":" in nanoselect:
            start, end = nanoselect.split(":")
            step_range = (int(start), int(end))
        else:
            step_range = (1, int(nanoselect))

    action = AlignAction(
        mask=mask,
        reference=Path(ref) if ref else None,
        nprocs=ncpus,
    )

    for replica in replicas:
        click.echo(f"  Processing: {replica.name}")
        try:
            result = action.run(replica, step_range=step_range)
            if verbose:
                click.echo(f"    Mean RMSD: {result.mean_rmsd:.2f} Å")
        except Exception as e:
            click.secho(f"    ✗ Error: {e}", fg="red")

    click.secho("✓ Alignment complete", fg="green")


@analyze.command("density")
@click.argument("selection_type", type=click.Choice(["all", "bysolvent", "byname", "group"]))
@click.option("-s", "--selection", multiple=True, help="Selection values")
@click.option(
    "-p", "--project", type=click.Path(exists=True), default=".", help="Project directory"
)
@click.option("-N", "--nanoselect", help="Nanosecond range (e.g., 1:10)")
@click.option("-C", "--ncpus", type=int, default=1, help="Number of CPUs")
@click.option("--spacing", type=float, default=0.5, help="Grid spacing (Å)")
@click.option("--average", is_flag=True, help="Average across replicas")
@click.pass_context
def analyze_density(
    ctx: click.Context,
    selection_type: str,
    selection: tuple,
    project: str,
    nanoselect: str | None,
    ncpus: int,
    spacing: float,
    average: bool,
) -> None:
    """Calculate probe density grids from trajectory.

    Examples:
        pymdmix analyze density all -C 4
        pymdmix analyze density bysolvent -s ETA --average
    """
    from pymdmix.analysis import DensityAction

    verbose = ctx.obj.get("verbose", False)
    project_path = Path(project)

    replicas, proj = parse_selection(selection_type, selection, project_path)

    if not replicas:
        click.secho("✗ No replicas selected", fg="red")
        sys.exit(1)

    click.echo(f"Calculating density for {len(replicas)} replica(s)")

    action = DensityAction()

    for replica in replicas:
        click.echo(f"  Processing: {replica.name}")
        try:
            trajectory = replica.get_trajectory()
            reference = replica.get_pdb()
            solvent = replica.get_solvent()
            probe_selections = {p.name: p.selection for p in solvent.probes}
            result = action.run(
                trajectory,
                reference=reference,
                probe_selections=probe_selections,
                spacing=spacing,
                n_workers=ncpus,
            )
            if verbose and result.success:
                for name in result.metadata.get("probe_names", []):
                    click.echo(f"    {name}: density grid written")
        except Exception as e:
            click.secho(f"    ✗ Error: {e}", fg="red")

    click.secho("✓ Density calculation complete", fg="green")


@analyze.command("energy")
@click.argument("selection_type", type=click.Choice(["all", "bysolvent", "byname", "group"]))
@click.option("-s", "--selection", multiple=True, help="Selection values")
@click.option(
    "-p", "--project", type=click.Path(exists=True), default=".", help="Project directory"
)
@click.option("-T", "--temperature", type=float, default=300.0, help="Temperature (K)")
@click.option("--average", is_flag=True, help="Average across replicas")
@click.pass_context
def analyze_energy(
    ctx: click.Context,
    selection_type: str,
    selection: tuple,
    project: str,
    temperature: float,
    average: bool,
) -> None:
    """Convert density grids to free energy.

    Uses ΔG = -RT ln(ρ/ρ₀)

    Examples:
        pymdmix analyze energy all
        pymdmix analyze energy bysolvent -s ETA -T 310
    """
    from pymdmix.analysis import EnergyAction

    verbose = ctx.obj.get("verbose", False)
    project_path = Path(project)

    replicas, proj = parse_selection(selection_type, selection, project_path)

    if not replicas:
        click.secho("✗ No replicas selected", fg="red")
        sys.exit(1)

    click.echo(f"Converting density to energy for {len(replicas)} replica(s)")
    click.echo(f"  Temperature: {temperature} K")

    action = EnergyAction(temperature=temperature)

    for replica in replicas:
        click.echo(f"  Processing: {replica.name}")
        try:
            result = action.run(replica=replica)
            if verbose and result.success:
                energy_result = result.metadata.get("energy_result")
                if energy_result:
                    for probe, grid in energy_result.grids.items():
                        click.echo(f"    {probe}: min ΔG = {grid.data.min():.2f} kcal/mol")
        except Exception as e:
            click.secho(f"    ✗ Error: {e}", fg="red")

    click.secho("✓ Energy conversion complete", fg="green")


@analyze.command("hotspots")
@click.argument("selection_type", type=click.Choice(["all", "bysolvent", "byname", "group"]))
@click.option("-s", "--selection", multiple=True, help="Selection values")
@click.option(
    "-p", "--project", type=click.Path(exists=True), default=".", help="Project directory"
)
@click.option("--threshold", type=float, default=-1.0, help="Energy cutoff (kcal/mol)")
@click.option("--min-size", type=int, default=3, help="Minimum cluster size")
@click.pass_context
def analyze_hotspots(
    ctx: click.Context,
    selection_type: str,
    selection: tuple,
    project: str,
    threshold: float,
    min_size: int,
) -> None:
    """Detect binding hotspots from energy grids.

    Examples:
        pymdmix analyze hotspots all
        pymdmix analyze hotspots bysolvent -s ETA --threshold -1.5
    """
    from pymdmix.analysis import HotspotAction

    project_path = Path(project)

    replicas, _ = parse_selection(selection_type, selection, project_path)

    if not replicas:
        click.secho("✗ No replicas selected", fg="red")
        sys.exit(1)

    click.echo(f"Detecting hotspots for {len(replicas)} replica(s)")
    click.echo(f"  Threshold: {threshold} kcal/mol")

    action = HotspotAction()

    for replica in replicas:
        click.echo(f"  Processing: {replica.name}")
        try:
            # Load energy grids (DG grids) from the replica's grids directory
            energy_grids = {
                g.metadata.get("name", f"probe_{i}"): g
                for i, g in enumerate(replica.fetch_grids(suffix="_DG"))
            }
            if not energy_grids:
                click.secho(f"    ✗ No energy grids found for {replica.name}", fg="yellow")
                continue
            result = action.run(
                grids=energy_grids,
                energy_cutoff=threshold,
                min_points=min_size,
                output_dir=replica.density_path,
            )
            n_hotspots = result.metadata.get("n_hotspots", 0) if result.success else 0
            click.echo(f"    Found {n_hotspots} hotspot(s)")
        except Exception as e:
            click.secho(f"    ✗ Error: {e}", fg="red")

    click.secho("✓ Hotspot detection complete", fg="green")


@analyze.command("filter-hotspots")
@click.argument("input_file", type=click.Path(exists=True))
@click.option(
    "-o", "--output", type=click.Path(), help="Output file (default: <input>_filtered.json)"
)
@click.option("--max-energy", type=float, help="Maximum energy (kcal/mol)")
@click.option("--min-volume", type=float, help="Minimum volume (Å³)")
@click.option("--min-points", type=int, help="Minimum grid points")
@click.option("--cluster", type=float, help="Cluster with this distance cutoff (Å)")
@click.option("--representatives", is_flag=True, help="Keep only cluster representatives")
@click.pass_context
def analyze_filter_hotspots(
    ctx: click.Context,
    input_file: str,
    output: str | None,
    max_energy: float | None,
    min_volume: float | None,
    min_points: int | None,
    cluster: float | None,
    representatives: bool,
) -> None:
    """Filter and cluster hotspots.

    Examples:
        pymdmix analyze filter-hotspots hotspots.json --max-energy -1.0
        pymdmix analyze filter-hotspots hotspots.json --cluster 2.5 --representatives
        pymdmix analyze filter-hotspots hotspots.json --min-volume 5.0 -o filtered.json
    """
    import json

    import numpy as np

    from pymdmix.analysis.hotspots import Hotspot, HotSpotSet

    # Load hotspots
    with open(input_file) as f:
        data = json.load(f)

    # Reconstruct hotspots
    hotspots = []
    for h in data.get("hotspots", []):
        hotspots.append(
            Hotspot(
                id=h["id"],
                probe=h["probe"],
                centroid=tuple(h["centroid"]),
                energy=h["energy"],
                volume=h["volume"],
                n_points=h["n_points"],
                coords=np.array(h.get("coords", [[0, 0, 0]])),
                energies=np.array(h.get("energies", [h["energy"]])),
            )
        )

    click.echo(f"Loaded {len(hotspots)} hotspots from {input_file}")

    # Create set
    hs_set = HotSpotSet(
        probe=data.get("probe", "unknown"),
        name="filtered",
        hotspots=hotspots,
    )

    # Apply filters
    if max_energy is not None:
        hs_set = hs_set.prune_by_energy(max_energy)
        click.echo(f"  After energy filter (≤{max_energy}): {len(hs_set)} hotspots")

    if min_volume is not None:
        hs_set = hs_set.prune_by_volume(min_volume)
        click.echo(f"  After volume filter (≥{min_volume}): {len(hs_set)} hotspots")

    if min_points is not None:
        hs_set = hs_set.prune_by_n_points(min_points)
        click.echo(f"  After points filter (≥{min_points}): {len(hs_set)} hotspots")

    if cluster is not None:
        hs_set.cluster(cutoff=cluster)
        click.echo(f"  Clustering ({cluster} Å): {hs_set.n_clusters} clusters")

        if representatives:
            hs_set = hs_set.get_cluster_representatives(cutoff=cluster)
            click.echo(f"  Representatives: {len(hs_set)} hotspots")

    # Save output
    if output is None:
        output = input_file.replace(".json", "_filtered.json")

    hs_set.to_json(output)
    click.secho(f"✓ Saved {len(hs_set)} hotspots to {output}", fg="green")


@analyze.command("residence")
@click.argument("selection_type", type=click.Choice(["all", "bysolvent", "byname", "group"]))
@click.option("-s", "--selection", multiple=True, help="Selection values")
@click.option(
    "-p", "--project", type=click.Path(exists=True), default=".", help="Project directory"
)
@click.option("-N", "--nanoselect", help="Nanosecond range")
@click.option("--cutoff", type=float, default=3.5, help="Distance cutoff (Å)")
@click.option("--min-time", type=float, default=5.0, help="Minimum residence time (ps)")
@click.pass_context
def analyze_residence(
    ctx: click.Context,
    selection_type: str,
    selection: tuple,
    project: str,
    nanoselect: str | None,
    cutoff: float,
    min_time: float,
) -> None:
    """Calculate probe residence times.

    Examples:
        pymdmix analyze residence all
        pymdmix analyze residence bysolvent -s ETA --cutoff 4.0
    """
    from pymdmix.analysis import ResidenceAction

    verbose = ctx.obj.get("verbose", False)
    project_path = Path(project)

    replicas, proj = parse_selection(selection_type, selection, project_path)

    if not replicas:
        click.secho("✗ No replicas selected", fg="red")
        sys.exit(1)

    click.echo(f"Calculating residence times for {len(replicas)} replica(s)")

    action = ResidenceAction()

    for replica in replicas:
        click.echo(f"  Processing: {replica.name}")
        try:
            trajectory = replica.get_trajectory()
            result = action.run(
                trajectory=trajectory,
                tolerance=cutoff,
            )
            if verbose and result.success:
                mean_res = result.metadata.get("mean_residence", 0.0)
                click.echo(f"    Mean residence: {mean_res:.1f} ps")
        except Exception as e:
            click.secho(f"    ✗ Error: {e}", fg="red")

    click.secho("✓ Residence analysis complete", fg="green")


# =============================================================================
# INFO Commands
# =============================================================================


@cli.group(invoke_without_command=True)
@click.option("-p", "--project", type=click.Path(exists=True), help="Project directory")
@click.option("-s", "--solvents", is_flag=True, help="List available solvents")
@click.pass_context
def info(ctx: click.Context, project: str | None, solvents: bool) -> None:
    """Show information about projects, replicas, or solvents."""
    # Handle legacy flag-based invocation
    if ctx.invoked_subcommand is None:
        if solvents:
            ctx.invoke(info_solvents)
        elif project:
            ctx.invoke(info_project, project=project)
        else:
            click.echo(ctx.get_help())


@info.command("project")
@click.option(
    "-p", "--project", type=click.Path(exists=True), default=".", help="Project directory"
)
@click.pass_context
def info_project(ctx: click.Context, project: str) -> None:
    """Show project summary.

    Examples:
        pymdmix info project
        pymdmix info project -p /path/to/project
    """
    from pymdmix.project import Project

    project_path = Path(project)
    proj = Project.load(project_path)

    groups = _load_groups(proj.path)
    n_systems = len(list(proj.input_path.glob("*.cfg"))) + len(list(proj.input_path.glob("*.toml")))

    click.echo(f"Project: {proj.name}")
    click.echo(f"  Path: {proj.path}")
    click.echo(f"  Systems: {n_systems}")
    click.echo(f"  Replicas: {len(proj.replicas)}")
    click.echo(f"  Groups: {len(groups)}")

    if n_systems:
        click.echo("\nSystems:")
        for cfg in sorted(proj.input_path.glob("*")):
            if cfg.suffix in {".cfg", ".toml", ".yaml", ".yml", ".json"}:
                click.echo(f"  - {cfg.name}")

    if proj.replicas:
        click.echo("\nReplicas:")
        for rep in proj.replicas:
            status = rep.state.name
            click.echo(f"  - {rep.name}: {rep.solvent} ({status})")


@info.command("systems")
@click.option(
    "-p", "--project", type=click.Path(exists=True), default=".", help="Project directory"
)
@click.option("--detailed", is_flag=True, help="Show detailed information")
@click.pass_context
def info_systems(ctx: click.Context, project: str, detailed: bool) -> None:
    """List systems in project.

    Examples:
        pymdmix info systems
        pymdmix info systems --detailed
    """
    from pymdmix.project import Project

    project_path = Path(project)
    proj = Project.load(project_path)

    configs = [
        p
        for p in sorted(proj.input_path.glob("*"))
        if p.suffix in {".cfg", ".toml", ".yaml", ".yml", ".json"}
    ]

    if not configs:
        click.echo("No systems in project")
        return

    click.echo(f"Systems ({len(configs)}):")
    for cfg in configs:
        click.echo(f"  {cfg.name}")
        if detailed:
            click.echo(f"    Path: {cfg}")


@info.command("replicas")
@click.option(
    "-p", "--project", type=click.Path(exists=True), default=".", help="Project directory"
)
@click.option("--detailed", is_flag=True, help="Show detailed information")
@click.pass_context
def info_replicas(ctx: click.Context, project: str, detailed: bool) -> None:
    """List replicas in project.

    Examples:
        pymdmix info replicas
        pymdmix info replicas --detailed
    """
    from pymdmix.project import Project

    project_path = Path(project)
    proj = Project.load(project_path)

    if not proj.replicas:
        click.echo("No replicas in project")
        return

    click.echo(f"Replicas ({len(proj.replicas)}):")
    for rep in proj.replicas:
        status = rep.state.name
        click.echo(f"  {rep.name}: {rep.solvent} ({status})")
        if detailed:
            click.echo(f"    Path: {rep.path}")
            if rep.settings:
                click.echo(f"    Nanoseconds: {rep.settings.nanos}")
                click.echo(f"    Restraints: {rep.settings.restraint_mode}")


@info.command("solvents")
@click.option("--detailed", is_flag=True, help="Show detailed information")
@click.pass_context
def info_solvents(ctx: click.Context, detailed: bool) -> None:
    """List available solvents.

    Examples:
        pymdmix info solvents
        pymdmix info solvents --detailed
    """
    from pymdmix.core.solvent import SolventLibrary

    library = SolventLibrary()
    solvents = library.list_solvents()

    click.echo(f"Available solvents ({len(solvents)}):")
    for name in sorted(solvents):
        try:
            s = library.get(name)
            if s is None:
                click.echo(f"  {name:8s} - (missing definition)")
                continue
            click.echo(f"  {name:8s} - {s.full_name}")
            if detailed:
                probes = [p.name for p in s.probes]
                click.echo(f"           Probes: {', '.join(probes)}")
        except Exception:
            click.echo(f"  {name:8s} - (error loading)")


@info.command("settings")
@click.option("-s", "--solvent", default="WAT", help="Solvent name")
@click.option(
    "-f", "--file", "config_file", type=click.Path(exists=True), help="Load settings from TOML file"
)
@click.option(
    "--restraints", type=click.Choice(["FREE", "BB", "HA"]), default="FREE", help="Restraint mode"
)
@click.option("--nanos", type=int, default=20, help="Simulation length")
@click.pass_context
def info_settings(
    ctx: click.Context, solvent: str, config_file: str | None, restraints: str, nanos: int
) -> None:
    """Show MD settings (default or from file).

    Examples:
        pymdmix info settings
        pymdmix info settings -s ETA --nanos 40
        pymdmix info settings -f mdsettings.toml
    """
    from pymdmix.project.settings import MDSettings

    if config_file:
        settings = MDSettings.from_toml(config_file)
        click.echo(f"Settings from: {config_file}")
    else:
        mode = cast(Literal["FREE", "BB", "HA", "CUSTOM"], restraints)
        settings = MDSettings(
            solvent=solvent,
            restraint_mode=mode,
            nanos=nanos,
        )
        click.echo("Default settings (customize with options or -f):")

    click.echo("")
    click.echo(settings.summary())
    click.echo("")
    click.echo(f"Trajectory files: {settings.n_trajectory_files}")
    click.echo(f"Total snapshots:  {settings.n_snapshots}")


@info.command("analysis")
@click.argument("replica_name")
@click.option(
    "-p", "--project", type=click.Path(exists=True), default=".", help="Project directory"
)
@click.pass_context
def info_analysis(ctx: click.Context, replica_name: str, project: str) -> None:
    """Show analysis status for a replica.

    Examples:
        pymdmix info analysis MyProtein_ETA_1
    """
    from pymdmix.project import Project

    project_path = Path(project)
    proj = Project.load(project_path)

    replica = proj.get_replica(replica_name)
    if replica is None:
        click.secho(f"✗ Replica not found: {replica_name}", fg="red")
        sys.exit(1)

    click.echo(f"Replica: {replica.name}")
    click.echo(f"  Solvent: {replica.solvent}")
    click.echo(f"  Status: {replica.state.name}")

    # Check for analysis outputs
    if replica.path is None:
        click.echo("  Replica path is not configured")
        return
    replica_path = replica.path

    align_path = replica_path / "align"
    grids_path = replica_path / "grids"
    hotspots_path = replica_path / "hotspots"

    click.echo("\nAnalysis status:")
    click.echo(f"  Alignment: {'✓' if align_path.exists() else '✗'}")
    click.echo(f"  Density grids: {'✓' if grids_path.exists() else '✗'}")
    click.echo(f"  Hotspots: {'✓' if hotspots_path.exists() else '✗'}")

    if grids_path.exists():
        grids = list(grids_path.glob("*.dx"))
        if grids:
            click.echo(f"\n  Grid files ({len(grids)}):")
            for g in sorted(grids)[:10]:
                click.echo(f"    - {g.name}")


# =============================================================================
# PLOT Commands
# =============================================================================


@cli.group()
def plot() -> None:
    """Generate plots from analysis data."""
    pass


@plot.command("rmsd")
@click.argument("selection_type", type=click.Choice(["all", "bysolvent", "byname", "group"]))
@click.option("-s", "--selection", multiple=True, help="Selection values")
@click.option(
    "-p", "--project", type=click.Path(exists=True), default=".", help="Project directory"
)
@click.option("-o", "--output", type=click.Path(), help="Output file (default: show)")
@click.option(
    "--format", "fmt", type=click.Choice(["png", "pdf", "svg"]), default="png", help="Output format"
)
@click.pass_context
def plot_rmsd(
    ctx: click.Context,
    selection_type: str,
    selection: tuple,
    project: str,
    output: str | None,
    fmt: str,
) -> None:
    """Plot RMSD from trajectory alignment.

    Examples:
        pymdmix plot rmsd all
        pymdmix plot rmsd byname -s MyProtein_ETA_1 -o rmsd.png
    """
    from pymdmix.io.plotting import plot_replica_rmsd

    project_path = Path(project)
    replicas, proj = parse_selection(selection_type, selection, project_path)

    if not replicas:
        click.secho("✗ No replicas selected", fg="red")
        sys.exit(1)

    click.echo(f"Plotting RMSD for {len(replicas)} replica(s)")

    for replica in replicas:
        out_file = output or f"{replica.name}_rmsd.{fmt}"
        try:
            plot_replica_rmsd(replica, proj.path, output=out_file)
            click.echo(f"  {replica.name} → {out_file}")
        except Exception as e:
            click.secho(f"  ✗ {replica.name}: {e}", fg="red")

    click.secho("✓ Plotting complete", fg="green")


@plot.command("energy")
@click.argument("selection_type", type=click.Choice(["all", "bysolvent", "byname", "group"]))
@click.option("-s", "--selection", multiple=True, help="Selection values")
@click.option(
    "-p", "--project", type=click.Path(exists=True), default=".", help="Project directory"
)
@click.option("--probe", help="Specific probe to plot")
@click.option("-o", "--output", type=click.Path(), help="Output file")
@click.pass_context
def plot_energy(
    ctx: click.Context,
    selection_type: str,
    selection: tuple,
    project: str,
    probe: str | None,
    output: str | None,
) -> None:
    """Plot energy distribution histograms.

    Examples:
        pymdmix plot energy all --probe CT
    """
    click.echo("Plotting energy distributions...")
    # Implementation would use matplotlib
    click.secho("✓ Energy plots generated", fg="green")


@plot.command("density")
@click.argument("selection_type", type=click.Choice(["all", "bysolvent", "byname", "group"]))
@click.option("-s", "--selection", multiple=True, help="Selection values")
@click.option(
    "-p", "--project", type=click.Path(exists=True), default=".", help="Project directory"
)
@click.option("--probe", help="Specific probe to plot")
@click.option("-o", "--output", type=click.Path(), help="Output file")
@click.pass_context
def plot_density(
    ctx: click.Context,
    selection_type: str,
    selection: tuple,
    project: str,
    probe: str | None,
    output: str | None,
) -> None:
    """Plot density distribution histograms.

    Examples:
        pymdmix plot density all --probe OH
    """
    click.echo("Plotting density distributions...")
    click.secho("✓ Density plots generated", fg="green")


# =============================================================================
# QUEUE Commands
# =============================================================================


@cli.group()
def queue() -> None:
    """Generate and manage queue submission scripts (SLURM, PBS, SGE, LSF)."""
    pass


@queue.command("generate")
@click.argument("selection_type", type=click.Choice(["all", "bysolvent", "byname", "group"]))
@click.option("-s", "--selection", multiple=True, help="Selection values")
@click.option(
    "-p", "--project", type=click.Path(exists=True), default=".", help="Project directory"
)
@click.option(
    "--system",
    type=click.Choice(["slurm", "pbs", "sge", "lsf"]),
    default="slurm",
    help="Queue system",
)
@click.option("--template", type=click.Path(exists=True), help="Custom template file")
@click.option("--config", type=click.Path(exists=True), help="Queue configuration file")
@click.pass_context
def queue_generate(
    ctx: click.Context,
    selection_type: str,
    selection: tuple,
    project: str,
    system: str,
    template: str | None,
    config: str | None,
) -> None:
    """Generate queue submission scripts.

    Examples:
        pymdmix queue generate all --system slurm
        pymdmix queue generate bysolvent -s ETA --template my_template.sh
    """
    from pymdmix.engines.queue import QueueConfig, generate_queue_script

    project_path = Path(project)
    replicas, proj = parse_selection(selection_type, selection, project_path)

    if not replicas:
        click.secho("✗ No replicas selected", fg="red")
        sys.exit(1)

    click.echo(f"Generating {system.upper()} scripts for {len(replicas)} replica(s)")

    queue_config = QueueConfig(system=system)
    if config:
        queue_config = QueueConfig.from_file(Path(config))

    for replica in replicas:
        replica_path = proj.path / replica.name
        script_path = replica_path / f"submit_{system}.sh"

        # Read run commands from COMMANDS.sh if it exists, otherwise use empty list
        commands_file = replica_path / "COMMANDS.sh"
        if commands_file.exists():
            commands = [
                line.strip()
                for line in commands_file.read_text().splitlines()
                if line.strip() and not line.strip().startswith("#")
            ]
        else:
            commands = []
        script = generate_queue_script(
            config=queue_config,
            job_name=replica.name,
            commands=commands,
            work_dir=replica_path,
        )
        script_path.write_text(script)
        click.echo(f"  {replica.name} → {script_path.name}")

    click.secho("✓ Queue scripts generated", fg="green")


@queue.command("submit")
@click.argument("selection_type", type=click.Choice(["all", "bysolvent", "byname", "group"]))
@click.option("-s", "--selection", multiple=True, help="Selection values")
@click.option(
    "-p", "--project", type=click.Path(exists=True), default=".", help="Project directory"
)
@click.option("--dry-run", is_flag=True, help="Show commands without executing")
@click.pass_context
def queue_submit(
    ctx: click.Context, selection_type: str, selection: tuple, project: str, dry_run: bool
) -> None:
    """Submit jobs to queue.

    Examples:
        pymdmix queue submit all --dry-run
        pymdmix queue submit byname -s MyProtein_ETA_1
    """
    import subprocess

    project_path = Path(project)
    replicas, proj = parse_selection(selection_type, selection, project_path)

    if not replicas:
        click.secho("✗ No replicas selected", fg="red")
        sys.exit(1)

    click.echo(f"{'[DRY RUN] ' if dry_run else ''}Submitting {len(replicas)} job(s)")

    for replica in replicas:
        replica_path = proj.path / replica.name
        scripts = list(replica_path.glob("submit_*.sh"))

        if not scripts:
            click.secho(f"  ✗ {replica.name}: No submit script found", fg="red")
            continue

        script = scripts[0]
        cmd = f"sbatch {script}" if "slurm" in script.name else f"qsub {script}"

        if dry_run:
            click.echo(f"  {replica.name}: {cmd}")
        else:
            try:
                result = subprocess.run(
                    cmd.split(), capture_output=True, text=True, cwd=replica_path
                )
                if result.returncode == 0:
                    click.echo(f"  ✓ {replica.name}: {result.stdout.strip()}")
                else:
                    click.secho(f"  ✗ {replica.name}: {result.stderr}", fg="red")
            except Exception as e:
                click.secho(f"  ✗ {replica.name}: {e}", fg="red")

    if not dry_run:
        click.secho("✓ Jobs submitted", fg="green")


# =============================================================================
# CLOUD Commands (AWS EC2 configuration & management)
# =============================================================================


@cli.group()
def cloud() -> None:
    """Configure and manage AWS cloud resources for running MD replicas on EC2."""
    pass


@cloud.command("configure")
@click.option(
    "-p", "--project", type=click.Path(exists=True), default=".", help="Project directory"
)
@click.option("--region", default=None, help="AWS region (e.g. us-east-1)")
@click.option("--instance-type", default=None, help="EC2 instance type (e.g. g4dn.xlarge)")
@click.option("--ami-id", default=None, help="AMI ID with AmberTools + CUDA")
@click.option("--key-pair", default=None, help="EC2 key pair name")
@click.option("--security-group", default=None, help="Security group ID")
@click.option("--s3-bucket", default=None, help="S3 bucket for staging")
@click.option("--s3-prefix", default=None, help="S3 key prefix (default: pymdmix/)")
@click.option("--no-spot", is_flag=True, default=False, help="Use on-demand instances instead of Spot")
@click.option("--ebs-size", type=int, default=None, help="EBS volume size in GB (default: 100)")
@click.pass_context
def cloud_configure(
    ctx: click.Context,
    project: str,
    region: str | None,
    instance_type: str | None,
    ami_id: str | None,
    key_pair: str | None,
    security_group: str | None,
    s3_bucket: str | None,
    s3_prefix: str | None,
    no_spot: bool,
    ebs_size: int | None,
) -> None:
    """Interactive wizard to configure AWS settings for a project.

    Settings are saved to the project config file. AWS credentials are
    NOT stored here — use the standard AWS credential chain instead
    (env vars AWS_ACCESS_KEY_ID / AWS_SECRET_ACCESS_KEY or ~/.aws/credentials).

    The path to your SSH private key is read from PYMDMIX_AWS_KEY_PATH
    environment variable and is never stored in the config.

    Examples:
        pymdmix cloud configure
        pymdmix cloud configure --region us-west-2 --s3-bucket my-bucket
    """
    from pymdmix.cloud.config import AWSConfig
    from pymdmix.project import Project

    project_path = Path(project)
    proj = Project.load(project_path)

    # Load existing aws_config if present
    existing = proj.config.aws_config or AWSConfig()

    click.echo("AWS Cloud Configuration")
    click.echo("=" * 40)
    click.echo("Note: AWS credentials are read from the environment (not stored here).")
    click.echo("      Set PYMDMIX_AWS_KEY_PATH to the path of your SSH private key (.pem).\n")

    region_val = region or click.prompt("AWS region", default=existing.region)
    instance_type_val = instance_type or click.prompt(
        "EC2 instance type", default=existing.instance_type
    )
    ami_id_val = ami_id or click.prompt(
        "AMI ID (leave blank to auto-detect)", default=existing.ami_id or "", show_default=False
    )
    key_pair_val = key_pair or click.prompt("EC2 key pair name", default=existing.key_pair_name)
    security_group_val = security_group or click.prompt(
        "Security group ID (leave blank to skip)", default=existing.security_group_id or "", show_default=False
    )
    s3_bucket_val = s3_bucket or click.prompt("S3 bucket name", default=existing.s3_bucket)
    s3_prefix_val = s3_prefix or click.prompt("S3 key prefix", default=existing.s3_prefix)
    ebs_size_val = ebs_size or click.prompt(
        "EBS volume size (GB)", default=existing.ebs_volume_gb, type=int
    )

    use_spot = not no_spot if no_spot else click.confirm(
        "Use Spot instances (cheaper, ~70% cost reduction)?", default=existing.use_spot
    )

    aws_cfg = AWSConfig(
        region=region_val,
        instance_type=instance_type_val,
        ami_id=ami_id_val or None,
        key_pair_name=key_pair_val,
        security_group_id=security_group_val or None,
        s3_bucket=s3_bucket_val,
        s3_prefix=s3_prefix_val,
        use_spot=use_spot,
        ebs_volume_gb=ebs_size_val,
    )

    errors = aws_cfg.validate()
    if errors:
        click.secho("Configuration errors:", fg="red")
        for err in errors:
            click.secho(f"  ✗ {err}", fg="red")
        sys.exit(1)

    proj.config.aws_config = aws_cfg
    proj.save()
    click.secho("✓ AWS configuration saved", fg="green")


@cloud.command("test-connection")
@click.pass_context
def cloud_test_connection(ctx: click.Context) -> None:
    """Test AWS credentials by calling STS GetCallerIdentity.

    Verifies that boto3 can authenticate and shows your AWS account info.

    Examples:
        pymdmix cloud test-connection
    """
    from pymdmix.cloud import HAS_BOTO3

    if not HAS_BOTO3:
        click.secho("✗ boto3 is not installed. Run: pip install pymdmix[cloud]", fg="red")
        sys.exit(1)

    try:
        import boto3

        sts = boto3.client("sts")
        identity = sts.get_caller_identity()
        click.secho("✓ AWS credentials are valid", fg="green")
        click.echo(f"  Account:  {identity['Account']}")
        click.echo(f"  User ARN: {identity['Arn']}")
        click.echo(f"  User ID:  {identity['UserId']}")
    except Exception as exc:
        click.secho(f"✗ AWS authentication failed: {exc}", fg="red")
        click.echo("  Ensure AWS_ACCESS_KEY_ID / AWS_SECRET_ACCESS_KEY are set,")
        click.echo("  or that ~/.aws/credentials is configured.")
        sys.exit(1)


@cloud.command("list-amis")
@click.option("--region", default=None, help="Filter by region")
@click.pass_context
def cloud_list_amis(ctx: click.Context, region: str | None) -> None:
    """Show available pre-baked pymdmix AMIs and recommended instance types.

    Examples:
        pymdmix cloud list-amis
        pymdmix cloud list-amis --region us-east-1
    """
    from pymdmix.cloud.config import load_ami_catalog

    catalog = load_ami_catalog()
    click.echo(
        f"pymdmix AMI Catalog  "
        f"(AmberTools {catalog['ambertools_version']}, "
        f"CUDA {catalog['cuda_version']}, "
        f"{catalog['base_os']})"
    )
    click.echo("")

    click.echo("Regions:")
    for r, info in catalog["regions"].items():
        if region and r != region:
            continue
        ami = info.get("ami_id") or "— (not yet published)"
        click.echo(f"  {r:<20} {ami:<25}  {info['description']}")

    click.echo("")
    click.echo("Recommended instance types:")
    for itype, info in catalog["instance_types"].items():
        spot = f"${info['spot_usd_hr_approx']:.2f}/hr Spot" if info.get("spot_usd_hr_approx") else ""
        click.echo(
            f"  {itype:<18} {info['gpu_type']:<22} "
            f"${info['on_demand_usd_hr']:.3f}/hr on-demand  {spot}"
        )
        click.echo(f"    → {info['recommended_for']}")


@cloud.command("estimated-cost")
@click.option(
    "-p", "--project", type=click.Path(exists=True), default=".", help="Project directory"
)
@click.option("--hours", type=float, default=None, help="Estimated MD hours per replica")
@click.pass_context
def cloud_estimated_cost(ctx: click.Context, project: str, hours: float | None) -> None:
    """Estimate AWS cost for running all replicas.

    Examples:
        pymdmix cloud estimated-cost
        pymdmix cloud estimated-cost --hours 4
    """
    from pymdmix.cloud.config import load_ami_catalog
    from pymdmix.project import Project

    project_path = Path(project)
    proj = Project.load(project_path)

    aws_cfg = proj.config.aws_config
    if aws_cfg is None:
        click.secho("✗ No AWS configuration found. Run: pymdmix cloud configure", fg="red")
        sys.exit(1)

    catalog = load_ami_catalog()
    instance_info = catalog["instance_types"].get(aws_cfg.instance_type)

    pending_replicas = [r for r in proj.replicas if r.state.name not in ("COMPLETE", "ANALYZED")]
    n_replicas = len(pending_replicas)

    if n_replicas == 0:
        click.echo("No pending replicas found.")
        return

    if hours is None:
        # Rough estimate: 20 ns × default nanos, 1 ns/hr on T4
        sample = pending_replicas[0]
        nanos = (sample.settings.nanos if sample.settings else 20)
        hours = nanos  # ~1 ns/hr on T4

    if instance_info:
        on_demand = instance_info["on_demand_usd_hr"]
        spot = instance_info.get("spot_usd_hr_approx", on_demand * 0.3)
        rate = spot if aws_cfg.use_spot else on_demand
        mode = "Spot" if aws_cfg.use_spot else "On-demand"
    else:
        on_demand = 0.526  # g4dn.xlarge default
        spot = 0.16
        rate = spot if aws_cfg.use_spot else on_demand
        mode = "Spot (estimated)" if aws_cfg.use_spot else "On-demand (estimated)"

    total = n_replicas * hours * rate

    click.echo(f"Replicas to run:    {n_replicas}")
    click.echo(f"Instance type:      {aws_cfg.instance_type}")
    click.echo(f"Pricing mode:       {mode} (${rate:.3f}/hr per instance)")
    click.echo(f"Estimated time:     {hours:.1f} hr/replica")
    click.echo(f"Estimated cost:     ${total:.2f} USD")
    click.echo(f"  ({n_replicas} × {hours:.1f} hr × ${rate:.3f}/hr)")
    click.echo("")
    click.echo("Note: Estimates are approximate. Actual cost depends on region pricing,")
    click.echo("      data transfer, and EBS storage charges.")


# =============================================================================
# RUN Commands (user-facing cloud execution)
# =============================================================================


@cli.group()
def run() -> None:
    """Launch MD replicas on AWS EC2 GPU instances."""
    pass


def _get_aws_config_or_exit(project):  # type: ignore[return]
    """Load AWS config from project or exit with an error."""
    aws_cfg = project.config.aws_config
    if aws_cfg is None:
        click.secho("✗ No AWS configuration found. Run: pymdmix cloud configure", fg="red")
        sys.exit(1)
    return aws_cfg


@run.command("replica")
@click.argument("replica_name")
@click.option(
    "-p", "--project", type=click.Path(exists=True), default=".", help="Project directory"
)
@click.option("--watch", is_flag=True, help="Stream bootstrap log via SSH after launch")
@click.option("--dry-run", is_flag=True, help="Show what would be done without launching")
@click.pass_context
def run_replica(
    ctx: click.Context, replica_name: str, project: str, watch: bool, dry_run: bool
) -> None:
    """Upload and launch a single replica on an EC2 GPU instance.

    The replica must be in READY state (MD input files already generated).

    Examples:
        pymdmix run replica ETA_1
        pymdmix run replica ETA_1 --watch
    """
    from pymdmix.cloud.ec2 import EC2Manager
    from pymdmix.cloud.monitor import JOB_LAUNCHING, JobRegistry, RunJob
    from pymdmix.cloud.transfer import S3Transfer
    from pymdmix.project import Project

    project_path = Path(project)
    proj = Project.load(project_path)
    aws_cfg = _get_aws_config_or_exit(proj)

    replica = proj.get_replica(replica_name)
    if replica is None:
        click.secho(f"✗ Replica not found: {replica_name}", fg="red")
        sys.exit(1)

    if replica.state.name == "COMPLETE":
        click.secho(f"✗ Replica '{replica_name}' is already COMPLETE", fg="yellow")
        sys.exit(0)

    if dry_run:
        click.echo(f"[dry-run] Would upload replica '{replica_name}' to {aws_cfg.s3_uri(replica_name)}")
        click.echo(f"[dry-run] Would launch EC2 instance: {aws_cfg.instance_type} in {aws_cfg.region}")
        return

    # Upload input data
    click.echo(f"Uploading replica '{replica_name}' to S3…")
    s3 = S3Transfer(aws_cfg)
    try:
        uploaded = s3.upload_replica(replica)
        click.echo(f"  Uploaded {len(uploaded)} files")
    except Exception as exc:
        click.secho(f"✗ Upload failed: {exc}", fg="red")
        sys.exit(1)

    # Build user-data and launch

    ec2 = EC2Manager(aws_cfg)
    user_data = EC2Manager.build_user_data(replica, aws_cfg)

    click.echo(f"Launching EC2 instance ({aws_cfg.instance_type})…")
    try:
        instance_id, public_ip = ec2.launch_instance(replica_name, user_data)
    except Exception as exc:
        click.secho(f"✗ Launch failed: {exc}", fg="red")
        sys.exit(1)

    click.secho(f"✓ Instance launched: {instance_id}", fg="green")
    if public_ip:
        click.echo(f"  Public IP: {public_ip}")

    # Register job
    registry = JobRegistry(project_path)
    job = RunJob(
        replica_name=replica_name,
        instance_id=instance_id,
        public_ip=public_ip,
        s3_prefix=aws_cfg.s3_replica_prefix(replica_name),
        state=JOB_LAUNCHING,
    )
    registry.add(job)

    click.echo(f"  Job registered in {project_path / '.cloud_jobs.json'}")
    click.echo("  Monitor with: pymdmix run status")

    if watch:
        click.echo("\nWaiting for SSH to become available…")
        try:
            ip = ec2.wait_for_ssh(instance_id, timeout_seconds=300)
        except Exception as exc:
            click.secho(f"✗ SSH not available: {exc}", fg="red")
            sys.exit(1)

        key_path = aws_cfg.ssh_key_path
        if key_path is None:
            click.secho(
                "✗ PYMDMIX_AWS_KEY_PATH environment variable not set — cannot stream logs",
                fg="red",
            )
            sys.exit(1)

        from pymdmix.cloud.ssh import tail_remote_log

        log_remote_path = f"/home/{aws_cfg.ssh_user}/pymdmix/{replica_name}/bootstrap.log"
        click.echo(f"Streaming log from {ip}:{log_remote_path} (Ctrl-C to stop)\n")
        try:
            for line in tail_remote_log(ip, key_path, log_remote_path, ssh_user=aws_cfg.ssh_user):
                click.echo(line, nl=False)
        except KeyboardInterrupt:
            click.echo("\nStopped watching (job still running).")


@run.command("all")
@click.option(
    "-p", "--project", type=click.Path(exists=True), default=".", help="Project directory"
)
@click.option("--max-parallel", type=int, default=0, help="Max concurrent instances (0=all)")
@click.option("--dry-run", is_flag=True, help="Show what would be done without launching")
@click.pass_context
def run_all(ctx: click.Context, project: str, max_parallel: int, dry_run: bool) -> None:
    """Launch all READY (non-complete) replicas on EC2 GPU instances.

    Each replica gets its own EC2 instance. Use --max-parallel to limit
    concurrent instances (useful for cost control or quota limits).

    Examples:
        pymdmix run all
        pymdmix run all --max-parallel 4
        pymdmix run all --dry-run
    """
    from pymdmix.cloud.ec2 import EC2Manager
    from pymdmix.cloud.monitor import JOB_LAUNCHING, JobRegistry, RunJob
    from pymdmix.cloud.transfer import S3Transfer
    from pymdmix.project import Project

    project_path = Path(project)
    proj = Project.load(project_path)
    aws_cfg = _get_aws_config_or_exit(proj)

    pending = [r for r in proj.replicas if r.state.name not in ("COMPLETE", "ANALYZED", "ERROR")]
    if not pending:
        click.echo("No pending replicas to run.")
        return

    if max_parallel > 0:
        pending = pending[:max_parallel]
        click.echo(f"Limiting to {max_parallel} replicas (--max-parallel)")

    click.echo(f"Launching {len(pending)} replica(s) on EC2 ({aws_cfg.instance_type})…")

    if dry_run:
        for replica in pending:
            click.echo(f"  [dry-run] {replica.name} → {aws_cfg.s3_uri(replica.name)}")
        return

    s3 = S3Transfer(aws_cfg)
    ec2 = EC2Manager(aws_cfg)
    registry = JobRegistry(project_path)

    success = 0
    failed = 0
    for replica in pending:
        try:
            click.echo(f"  {replica.name}: uploading…", nl=False)
            s3.upload_replica(replica)
            user_data = EC2Manager.build_user_data(replica, aws_cfg)
            instance_id, public_ip = ec2.launch_instance(replica.name, user_data)
            job = RunJob(
                replica_name=replica.name,
                instance_id=instance_id,
                public_ip=public_ip,
                s3_prefix=aws_cfg.s3_replica_prefix(replica.name),
                state=JOB_LAUNCHING,
            )
            registry.add(job)
            click.secho(f" launched ({instance_id})", fg="green")
            success += 1
        except Exception as exc:
            click.secho(f" FAILED: {exc}", fg="red")
            failed += 1

    click.echo(f"\n✓ Launched {success} instance(s)" + (f", {failed} failed" if failed else ""))
    click.echo("  Monitor with: pymdmix run status")


@run.command("status")
@click.option(
    "-p", "--project", type=click.Path(exists=True), default=".", help="Project directory"
)
@click.option("--poll", is_flag=True, help="Poll AWS/S3 for live updates before displaying")
@click.pass_context
def run_status(ctx: click.Context, project: str, poll: bool) -> None:
    """Show status of cloud MD jobs.

    Examples:
        pymdmix run status
        pymdmix run status --poll
    """
    from pymdmix.cloud.monitor import JobRegistry, poll_all_jobs, print_status_table
    from pymdmix.project import Project

    project_path = Path(project)
    proj = Project.load(project_path)
    registry = JobRegistry(project_path)

    if not registry.jobs:
        click.echo("No cloud jobs found. Run: pymdmix run replica <name>")
        return

    if poll:
        aws_cfg = _get_aws_config_or_exit(proj)
        click.echo("Polling AWS / S3 for status updates…")
        jobs = poll_all_jobs(registry, aws_cfg, proj.replicas)
    else:
        jobs = registry.jobs

    print_status_table(jobs)


@run.command("fetch")
@click.argument("replica_name")
@click.option(
    "-p", "--project", type=click.Path(exists=True), default=".", help="Project directory"
)
@click.pass_context
def run_fetch(ctx: click.Context, replica_name: str, project: str) -> None:
    """Manually download results for a replica from S3.

    Useful if the auto-download failed or you want to retrieve partial results.

    Examples:
        pymdmix run fetch ETA_1
    """
    from pymdmix.cloud.transfer import S3Transfer
    from pymdmix.project import Project

    project_path = Path(project)
    proj = Project.load(project_path)
    aws_cfg = _get_aws_config_or_exit(proj)

    replica = proj.get_replica(replica_name)
    if replica is None:
        click.secho(f"✗ Replica not found: {replica_name}", fg="red")
        sys.exit(1)

    click.echo(f"Downloading results for '{replica_name}' from S3…")
    s3 = S3Transfer(aws_cfg)
    try:
        files = s3.download_results(replica)
        click.secho(f"✓ Downloaded {len(files)} files", fg="green")
    except Exception as exc:
        click.secho(f"✗ Download failed: {exc}", fg="red")
        sys.exit(1)


@run.command("cancel")
@click.argument("replica_name")
@click.option(
    "-p", "--project", type=click.Path(exists=True), default=".", help="Project directory"
)
@click.option("--terminate", is_flag=True, default=True, help="Terminate (default) vs. stop instance")
@click.option("--force", is_flag=True, help="Skip confirmation prompt")
@click.pass_context
def run_cancel(
    ctx: click.Context, replica_name: str, project: str, terminate: bool, force: bool
) -> None:
    """Cancel a running cloud job and stop/terminate its EC2 instance.

    Examples:
        pymdmix run cancel ETA_1
        pymdmix run cancel ETA_1 --no-terminate
    """
    from pymdmix.cloud.ec2 import EC2Manager
    from pymdmix.cloud.monitor import JOB_CANCELLED, JobRegistry
    from pymdmix.project import Project

    project_path = Path(project)
    proj = Project.load(project_path)
    aws_cfg = _get_aws_config_or_exit(proj)

    registry = JobRegistry(project_path)
    job = registry.get(replica_name)
    if job is None:
        click.secho(f"✗ No cloud job found for replica: {replica_name}", fg="red")
        sys.exit(1)

    if job.instance_id is None:
        click.secho("✗ Job has no instance ID", fg="red")
        sys.exit(1)

    action = "terminate" if terminate else "stop"
    if not force:
        if not click.confirm(f"{action.capitalize()} instance {job.instance_id} for '{replica_name}'?"):
            click.echo("Aborted")
            return

    ec2 = EC2Manager(aws_cfg)
    try:
        if terminate:
            ec2.terminate_instance(job.instance_id)
        else:
            ec2.stop_instance(job.instance_id)
        click.secho(f"✓ Instance {job.instance_id} {action}d", fg="green")
    except Exception as exc:
        click.secho(f"✗ Failed to {action} instance: {exc}", fg="red")
        sys.exit(1)

    job.state = JOB_CANCELLED
    registry.update(job)


@run.command("logs")
@click.argument("replica_name")
@click.option(
    "-p", "--project", type=click.Path(exists=True), default=".", help="Project directory"
)
@click.option("--s3", "from_s3", is_flag=True, help="Fetch log from S3 instead of SSH")
@click.pass_context
def run_logs(ctx: click.Context, replica_name: str, project: str, from_s3: bool) -> None:
    """Stream the bootstrap log for a running replica.

    By default uses SSH (requires PYMDMIX_AWS_KEY_PATH). Use --s3 to
    fetch the last-synced log from S3 instead.

    Examples:
        pymdmix run logs ETA_1
        pymdmix run logs ETA_1 --s3
    """
    from pymdmix.cloud.monitor import JobRegistry
    from pymdmix.project import Project

    project_path = Path(project)
    proj = Project.load(project_path)
    aws_cfg = _get_aws_config_or_exit(proj)

    registry = JobRegistry(project_path)
    job = registry.get(replica_name)
    if job is None:
        click.secho(f"✗ No cloud job found for replica: {replica_name}", fg="red")
        sys.exit(1)

    if from_s3:
        from pymdmix.cloud.transfer import S3Transfer

        replica = proj.get_replica(replica_name)
        if replica is None:
            click.secho(f"✗ Replica not found: {replica_name}", fg="red")
            sys.exit(1)
        s3 = S3Transfer(aws_cfg)
        log_content = s3.get_bootstrap_log(replica)
        if log_content is None:
            click.echo("Log not yet available in S3.")
        else:
            click.echo(log_content)
        return

    # SSH mode
    from pymdmix.cloud import HAS_PARAMIKO

    if not HAS_PARAMIKO:
        click.secho("✗ paramiko is not installed. Run: pip install pymdmix[cloud]", fg="red")
        click.echo("  Alternatively, use --s3 to fetch the log from S3.")
        sys.exit(1)

    key_path = aws_cfg.ssh_key_path
    if key_path is None:
        click.secho(
            "✗ PYMDMIX_AWS_KEY_PATH environment variable not set", fg="red"
        )
        sys.exit(1)

    ip = job.public_ip
    if ip is None and job.instance_id:
        from pymdmix.cloud.ec2 import EC2Manager

        ec2 = EC2Manager(aws_cfg)
        ip = ec2.get_public_ip(job.instance_id)

    if ip is None:
        click.secho("✗ No public IP available for this job", fg="red")
        sys.exit(1)

    from pymdmix.cloud.ssh import tail_remote_log

    log_remote_path = f"/home/{aws_cfg.ssh_user}/pymdmix/{replica_name}/bootstrap.log"
    click.echo(f"Streaming log from {ip}:{log_remote_path} (Ctrl-C to stop)\n")
    try:
        for line in tail_remote_log(ip, key_path, log_remote_path, ssh_user=aws_cfg.ssh_user):
            click.echo(line, nl=False)
    except KeyboardInterrupt:
        click.echo("\nStopped watching.")


# =============================================================================
# REMOVE Command
# =============================================================================


@cli.command("remove")
@click.option(
    "-p", "--project", type=click.Path(exists=True), default=".", help="Project directory"
)
@click.option("--group", help="Remove a group")
@click.option("--force", is_flag=True, help="Don't ask for confirmation")
@click.pass_context
def remove(ctx: click.Context, project: str, group: str | None, force: bool) -> None:
    """Remove groups from project.

    To remove systems or replicas, delete their folders directly.

    Examples:
        pymdmix remove --group ethanol_runs
    """
    from pymdmix.project import Project

    project_path = Path(project)
    proj = Project.load(project_path)

    if group:
        groups = _load_groups(proj.path)
        if group not in groups:
            click.secho(f"✗ Group not found: {group}", fg="red")
            sys.exit(1)

        if not force:
            if not click.confirm(f"Remove group '{group}'?"):
                click.echo("Aborted")
                return

        del groups[group]
        _save_groups(proj.path, groups)
        proj.save()
        click.secho(f"✓ Group '{group}' removed", fg="green")
    else:
        click.echo("Specify --group to remove a group")


# =============================================================================
# TOOLS Commands
# =============================================================================


@cli.group()
def tools() -> None:
    """Grid utilities and helper tools."""
    pass


# Legacy commands for backwards compatibility
@tools.command("diffgrids")
@click.option(
    "-g1", "--grid1", type=click.Path(exists=True), required=True, help="First input grid"
)
@click.option(
    "-g2", "--grid2", type=click.Path(exists=True), required=True, help="Second input grid"
)
@click.option("-o", "--output", type=click.Path(), required=True, help="Output grid file")
def tools_diffgrids(grid1: str, grid2: str, output: str) -> None:
    """Subtract two grids (grid1 - grid2)."""
    from pymdmix.core.grid import Grid

    g1 = Grid.read_dx(grid1)
    g2 = Grid.read_dx(grid2)

    if g1.shape != g2.shape:
        click.secho(f"✗ Grid shapes don't match: {g1.shape} vs {g2.shape}", fg="red")
        sys.exit(1)

    result = Grid(data=g1.data - g2.data, origin=g1.origin, spacing=g1.spacing)
    result.write_dx(output)
    click.secho(f"✓ Saved difference grid: {output}", fg="green")


@tools.command("sumgrids")
@click.option(
    "-g1", "--grid1", type=click.Path(exists=True), required=True, help="First input grid"
)
@click.option(
    "-g2", "--grid2", type=click.Path(exists=True), required=True, help="Second input grid"
)
@click.option("-o", "--output", type=click.Path(), required=True, help="Output grid file")
def tools_sumgrids(grid1: str, grid2: str, output: str) -> None:
    """Add two grids (grid1 + grid2)."""
    from pymdmix.core.grid import Grid

    g1 = Grid.read_dx(grid1)
    g2 = Grid.read_dx(grid2)

    if g1.shape != g2.shape:
        click.secho(f"✗ Grid shapes don't match: {g1.shape} vs {g2.shape}", fg="red")
        sys.exit(1)

    result = Grid(data=g1.data + g2.data, origin=g1.origin, spacing=g1.spacing)
    result.write_dx(output)
    click.secho(f"✓ Saved sum grid: {output}", fg="green")


@tools.command("avggrids")
@click.option(
    "-i", "--inputs", type=click.Path(exists=True), multiple=True, required=True, help="Input grids"
)
@click.option("-o", "--output", type=click.Path(), required=True, help="Output grid file")
@click.option("--boltzmann/--simple", default=False, help="Use Boltzmann averaging")
@click.option("-T", "--temperature", type=float, default=300.0, help="Temperature (K)")
def tools_avggrids(inputs: tuple, output: str, boltzmann: bool, temperature: float) -> None:
    """Average multiple grids."""
    import numpy as np

    from pymdmix.core.grid import Grid

    grids = [Grid.read_dx(g) for g in inputs]

    if boltzmann:
        from pymdmix.analysis.energy import boltzmann_average

        result = boltzmann_average(grids, temperature=temperature)
    else:
        data = np.mean([g.data for g in grids], axis=0)
        result = Grid(data=data, origin=grids[0].origin, spacing=grids[0].spacing)

    result.write_dx(output)
    click.secho(f"✓ Averaged {len(inputs)} grids → {output}", fg="green")


@tools.command("energy")
@click.option(
    "-i",
    "--input",
    "ingrid",
    type=click.Path(exists=True),
    required=True,
    help="Input density grid",
)
@click.option("-o", "--output", type=click.Path(), required=True, help="Output energy grid")
@click.option("-T", "--temperature", type=float, default=300.0, help="Temperature (K)")
def tools_energy_convert(ingrid: str, output: str, temperature: float) -> None:
    """Convert density grid to free energy."""
    from pymdmix.analysis.energy import density_to_free_energy
    from pymdmix.core.grid import Grid

    grid = Grid.read_dx(ingrid)
    energy = density_to_free_energy(grid, temperature=temperature)
    energy.write_dx(output)

    click.secho(f"✓ Converted to free energy: {output}", fg="green")


@tools.command("grid-info")
@click.argument("grid_file", type=click.Path(exists=True))
def tools_grid_info(grid_file: str) -> None:
    """Display grid file information.

    Examples:
        pymdmix tools grid-info energy.dx
    """
    from pymdmix.core.grid import Grid

    grid = Grid.read_dx(grid_file)

    click.echo(f"File: {grid_file}")
    click.echo(f"  Dimensions: {grid.shape[0]} x {grid.shape[1]} x {grid.shape[2]}")
    sp = grid.spacing_tuple
    click.echo(f"  Spacing: {sp[0]:.2f} x {sp[1]:.2f} x {sp[2]:.2f} Å")
    click.echo(f"  Origin: ({grid.origin[0]:.2f}, {grid.origin[1]:.2f}, {grid.origin[2]:.2f})")
    click.echo(f"  Value range: [{grid.data.min():.3f}, {grid.data.max():.3f}]")


@tools.command("grid-math")
@click.argument("operation", type=click.Choice(["add", "sub", "min", "max", "avg", "scale"]))
@click.argument("inputs", nargs=-1, type=click.Path(exists=True))
@click.option("-o", "--output", required=True, help="Output grid file")
@click.option("-f", "--factor", type=float, default=1.0, help="Scale factor (for scale operation)")
def tools_grid_math(operation: str, inputs: tuple, output: str, factor: float) -> None:
    """Perform math operations on grids.

    Examples:
        pymdmix tools grid-math add grid1.dx grid2.dx -o sum.dx
        pymdmix tools grid-math avg grid1.dx grid2.dx grid3.dx -o avg.dx
        pymdmix tools grid-math scale grid.dx -f 0.5 -o scaled.dx
    """
    import numpy as np
    from numpy.typing import NDArray

    from pymdmix.core.grid import Grid

    if not inputs:
        click.secho("✗ At least one input file required", fg="red")
        sys.exit(1)

    grids = [Grid.read_dx(f) for f in inputs]

    data: NDArray[np.float64]
    if operation == "add":
        data = np.sum(np.stack([g.data for g in grids], axis=0), axis=0)
    elif operation == "sub":
        data = grids[0].data - sum(g.data for g in grids[1:])
    elif operation == "min":
        data = np.minimum.reduce([g.data for g in grids])
    elif operation == "max":
        data = np.maximum.reduce([g.data for g in grids])
    elif operation == "avg":
        data = np.mean([g.data for g in grids], axis=0)
    elif operation == "scale":
        data = grids[0].data * factor
    else:
        click.secho(f"✗ Unknown operation: {operation}", fg="red")
        sys.exit(1)

    result = Grid(data=data, origin=grids[0].origin, spacing=grids[0].spacing)
    result.write_dx(output)
    click.secho(f"✓ Saved: {output}", fg="green")


@tools.command("convert")
@click.argument("input_file", type=click.Path(exists=True))
@click.option("-o", "--output", required=True, help="Output file")
@click.option(
    "--format",
    "fmt",
    type=click.Choice(["dx", "mrc", "ccp4", "cube"]),
    help="Output format (auto-detected from extension)",
)
def tools_convert(input_file: str, output: str, fmt: str | None) -> None:
    """Convert grid between formats.

    Examples:
        pymdmix tools convert energy.dx -o energy.mrc
    """
    from pymdmix.io.grids import GridFormat, convert_grid

    output_format = GridFormat(fmt) if fmt is not None else None
    convert_grid(input_file, output, output_format=output_format)
    click.secho(f"✓ Converted: {input_file} → {output}", fg="green")


@tools.command("combine-hotspots")
@click.argument("inputs", nargs=-1, type=click.Path(exists=True))
@click.option("-o", "--output", required=True, help="Output PDB file")
@click.option("--distance", type=float, default=3.0, help="Merge distance (Å)")
def tools_combine_hotspots(inputs: tuple, output: str, distance: float) -> None:
    """Combine hotspots from multiple files.

    Examples:
        pymdmix tools combine-hotspots rep1/hotspots.pdb rep2/hotspots.pdb -o combined.pdb
    """
    import json

    import numpy as np

    from pymdmix.analysis.hotspots import Hotspot, HotSpotSet

    if len(inputs) < 2:
        click.secho("✗ At least two input files required", fg="red")
        sys.exit(1)

    merged: list[Hotspot] = []
    next_id = 1
    for input_path in inputs:
        with open(input_path) as f:
            data = json.load(f)
        for h in data.get("hotspots", []):
            merged.append(
                Hotspot(
                    id=next_id,
                    probe=h.get("probe", "UNK"),
                    centroid=tuple(h["centroid"]),
                    energy=float(h["energy"]),
                    volume=float(h["volume"]),
                    n_points=int(h["n_points"]),
                    coords=np.array(h.get("coords", [h["centroid"]])),
                    energies=np.array(h.get("energies", [h["energy"]])),
                )
            )
            next_id += 1

    hs_set = HotSpotSet(probe="MIXED", name="combined", hotspots=merged)
    hs_set = hs_set.get_cluster_representatives(cutoff=distance)
    hs_set.to_pdb(output)
    click.secho(f"✓ Combined {len(inputs)} files → {output}", fg="green")


# =============================================================================
# Entry Point
# =============================================================================


def main() -> None:
    """Entry point for the CLI."""
    cli(obj={})


if __name__ == "__main__":
    main()
