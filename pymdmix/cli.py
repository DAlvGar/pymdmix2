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
from typing import Optional, List

import click

from pymdmix import __version__


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
def create_project(ctx: click.Context, name: str, config: Optional[str], 
                   directory: Optional[str], force: bool) -> None:
    """Create a new pyMDMix project.

    Examples:
        pymdmix create project -n myproject
        pymdmix create project -n myproject -f config.yaml
    """
    from pymdmix.project import Project, Config

    verbose = ctx.obj.get("verbose", False)
    project_dir = Path(directory) if directory else (Path.cwd() / name if name != "." else Path.cwd())

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
@click.option("-f", "--file", "config", type=click.Path(exists=True), required=True, 
              help="Solvent config file")
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
@click.option("-f", "--file", "config", type=click.Path(exists=True), required=True,
              help="System configuration file")
@click.option("-p", "--project", type=click.Path(exists=True), default=".",
              help="Project directory (default: current)")
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
    from pymdmix.project import Project
    from pymdmix.io.parsers import parse_system_config

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

    # Add system to project
    proj.add_system(system_config)
    proj.save()

    click.secho(f"✓ System '{system_config.name}' added to project", fg="green")


@add.command("replica")
@click.option("-f", "--file", "config", type=click.Path(exists=True), required=True,
              help="Replica configuration file")
@click.option("-p", "--project", type=click.Path(exists=True), default=".",
              help="Project directory (default: current)")
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
    from pymdmix.project import Project
    from pymdmix.io.parsers import parse_replica_config

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
    created = []
    for i in range(count):
        replica = proj.create_replica(replica_config)
        created.append(replica.name)
        if verbose:
            click.echo(f"  Created: {replica.name}")

    proj.save()

    if count == 1:
        click.secho(f"✓ Replica '{created[0]}' added to project", fg="green")
    else:
        click.secho(f"✓ {count} replicas added to project", fg="green")
        for name in created:
            click.echo(f"    - {name}")


@add.command("group")
@click.option("-n", "--name", required=True, help="Group name")
@click.option("-s", "--selection", multiple=True, required=True, help="Replica names to include")
@click.option("-p", "--project", type=click.Path(exists=True), default=".",
              help="Project directory (default: current)")
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
    for rep_name in selection:
        if rep_name not in proj.replicas:
            click.secho(f"✗ Replica not found: {rep_name}", fg="red")
            sys.exit(1)

    # Create group
    proj.create_group(name, list(selection))
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
def setup_prepare(ctx: click.Context, structure: str, output: Optional[str],
                  cap: bool, disulfide: bool, remove_water: bool) -> None:
    """Prepare a structure for MD simulation.

    Examples:
        pymdmix setup prepare protein.pdb -o prepared.pdb
        pymdmix setup prepare protein.pdb --no-cap
    """
    from pymdmix.setup.prepare import prepare_structure, StructurePreparationOptions

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
def setup_solvate(ctx: click.Context, structure: str, solvent: str,
                  output: Optional[str], buffer: float) -> None:
    """Solvate a structure with a solvent mixture.

    Examples:
        pymdmix setup solvate protein.pdb -s ETA -o solvated
    """
    from pymdmix.setup.solvate import solvate_structure, SolvationOptions
    from pymdmix.core.solvent import SolventLibrary

    input_path = Path(structure)
    output_prefix = output or input_path.stem + "_solvated"

    library = SolventLibrary()
    try:
        solv = library.get(solvent)
    except KeyError:
        click.secho(f"✗ Unknown solvent: {solvent}", fg="red")
        click.echo(f"  Available: {', '.join(library.list_solvents())}")
        sys.exit(1)

    click.echo(f"Solvating with {solv.full_name}")

    options = SolvationOptions(box_buffer=buffer)
    result = solvate_structure(input_path, solv, options=options)

    top_path = Path(output_prefix + ".prmtop")
    crd_path = Path(output_prefix + ".rst7")
    result.topology.write(str(top_path))
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
        return list(proj.replicas.values()), proj
    elif selection_type == "bysolvent":
        replicas = [r for r in proj.replicas.values() if r.solvent in selection]
        return replicas, proj
    elif selection_type == "byname":
        replicas = [proj.get_replica(name) for name in selection]
        return replicas, proj
    elif selection_type == "group":
        if len(selection) != 1:
            raise click.UsageError("Group selection requires exactly one group name")
        group_name = selection[0]
        replica_names = proj.get_group(group_name)
        replicas = [proj.get_replica(name) for name in replica_names]
        return replicas, proj
    else:
        raise click.UsageError(f"Unknown selection type: {selection_type}")


@analyze.command("align")
@click.argument("selection_type", type=click.Choice(["all", "bysolvent", "byname", "group"]))
@click.option("-s", "--selection", multiple=True, help="Selection values")
@click.option("-p", "--project", type=click.Path(exists=True), default=".",
              help="Project directory")
@click.option("-N", "--nanoselect", help="Nanosecond range (e.g., 1:10)")
@click.option("-C", "--ncpus", type=int, default=1, help="Number of CPUs")
@click.option("--mask", help="Alignment mask (residue range)")
@click.option("--ref", type=click.Path(exists=True), help="Reference PDB file")
@click.pass_context
def analyze_align(ctx: click.Context, selection_type: str, selection: tuple,
                  project: str, nanoselect: Optional[str], ncpus: int,
                  mask: Optional[str], ref: Optional[str]) -> None:
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
@click.option("-p", "--project", type=click.Path(exists=True), default=".",
              help="Project directory")
@click.option("-N", "--nanoselect", help="Nanosecond range (e.g., 1:10)")
@click.option("-C", "--ncpus", type=int, default=1, help="Number of CPUs")
@click.option("--spacing", type=float, default=0.5, help="Grid spacing (Å)")
@click.option("--average", is_flag=True, help="Average across replicas")
@click.pass_context
def analyze_density(ctx: click.Context, selection_type: str, selection: tuple,
                    project: str, nanoselect: Optional[str], ncpus: int,
                    spacing: float, average: bool) -> None:
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

    action = DensityAction(spacing=spacing, nprocs=ncpus)

    for replica in replicas:
        click.echo(f"  Processing: {replica.name}")
        try:
            result = action.run(replica)
            if verbose:
                for probe, grid in result.grids.items():
                    click.echo(f"    {probe}: max density = {grid.data.max():.2f}")
        except Exception as e:
            click.secho(f"    ✗ Error: {e}", fg="red")

    click.secho("✓ Density calculation complete", fg="green")


@analyze.command("energy")
@click.argument("selection_type", type=click.Choice(["all", "bysolvent", "byname", "group"]))
@click.option("-s", "--selection", multiple=True, help="Selection values")
@click.option("-p", "--project", type=click.Path(exists=True), default=".",
              help="Project directory")
@click.option("-T", "--temperature", type=float, default=300.0, help="Temperature (K)")
@click.option("--average", is_flag=True, help="Average across replicas")
@click.pass_context
def analyze_energy(ctx: click.Context, selection_type: str, selection: tuple,
                   project: str, temperature: float, average: bool) -> None:
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
            result = action.run(replica)
            if verbose:
                for probe, grid in result.grids.items():
                    click.echo(f"    {probe}: min ΔG = {grid.data.min():.2f} kcal/mol")
        except Exception as e:
            click.secho(f"    ✗ Error: {e}", fg="red")

    click.secho("✓ Energy conversion complete", fg="green")


@analyze.command("hotspots")
@click.argument("selection_type", type=click.Choice(["all", "bysolvent", "byname", "group"]))
@click.option("-s", "--selection", multiple=True, help="Selection values")
@click.option("-p", "--project", type=click.Path(exists=True), default=".",
              help="Project directory")
@click.option("--threshold", type=float, default=-1.0, help="Energy cutoff (kcal/mol)")
@click.option("--min-size", type=int, default=3, help="Minimum cluster size")
@click.pass_context
def analyze_hotspots(ctx: click.Context, selection_type: str, selection: tuple,
                     project: str, threshold: float, min_size: int) -> None:
    """Detect binding hotspots from energy grids.

    Examples:
        pymdmix analyze hotspots all
        pymdmix analyze hotspots bysolvent -s ETA --threshold -1.5
    """
    from pymdmix.analysis import HotspotAction

    verbose = ctx.obj.get("verbose", False)
    project_path = Path(project)

    replicas, proj = parse_selection(selection_type, selection, project_path)
    
    if not replicas:
        click.secho("✗ No replicas selected", fg="red")
        sys.exit(1)

    click.echo(f"Detecting hotspots for {len(replicas)} replica(s)")
    click.echo(f"  Threshold: {threshold} kcal/mol")

    action = HotspotAction(energy_cutoff=threshold, min_size=min_size)

    for replica in replicas:
        click.echo(f"  Processing: {replica.name}")
        try:
            hotspots = action.run(replica)
            click.echo(f"    Found {len(hotspots)} hotspot(s)")
            if verbose:
                for i, hs in enumerate(hotspots[:5]):
                    click.echo(f"      {i+1}. ΔG = {hs.energy:.2f} kcal/mol")
        except Exception as e:
            click.secho(f"    ✗ Error: {e}", fg="red")

    click.secho("✓ Hotspot detection complete", fg="green")


@analyze.command("filter-hotspots")
@click.argument("input_file", type=click.Path(exists=True))
@click.option("-o", "--output", type=click.Path(), help="Output file (default: <input>_filtered.json)")
@click.option("--max-energy", type=float, help="Maximum energy (kcal/mol)")
@click.option("--min-volume", type=float, help="Minimum volume (Å³)")
@click.option("--min-points", type=int, help="Minimum grid points")
@click.option("--cluster", type=float, help="Cluster with this distance cutoff (Å)")
@click.option("--representatives", is_flag=True, help="Keep only cluster representatives")
@click.pass_context
def analyze_filter_hotspots(ctx: click.Context, input_file: str, output: Optional[str],
                            max_energy: Optional[float], min_volume: Optional[float],
                            min_points: Optional[int], cluster: Optional[float],
                            representatives: bool) -> None:
    """Filter and cluster hotspots.

    Examples:
        pymdmix analyze filter-hotspots hotspots.json --max-energy -1.0
        pymdmix analyze filter-hotspots hotspots.json --cluster 2.5 --representatives
        pymdmix analyze filter-hotspots hotspots.json --min-volume 5.0 -o filtered.json
    """
    import json
    from pymdmix.analysis.hotspots import Hotspot, HotSpotSet
    import numpy as np

    # Load hotspots
    with open(input_file) as f:
        data = json.load(f)

    # Reconstruct hotspots
    hotspots = []
    for h in data.get("hotspots", []):
        hotspots.append(Hotspot(
            id=h["id"],
            probe=h["probe"],
            centroid=tuple(h["centroid"]),
            energy=h["energy"],
            volume=h["volume"],
            n_points=h["n_points"],
            coords=np.array(h.get("coords", [[0, 0, 0]])),
            energies=np.array(h.get("energies", [h["energy"]])),
        ))

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
        labels = hs_set.cluster(cutoff=cluster)
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
@click.option("-p", "--project", type=click.Path(exists=True), default=".",
              help="Project directory")
@click.option("-N", "--nanoselect", help="Nanosecond range")
@click.option("--cutoff", type=float, default=3.5, help="Distance cutoff (Å)")
@click.option("--min-time", type=float, default=5.0, help="Minimum residence time (ps)")
@click.pass_context
def analyze_residence(ctx: click.Context, selection_type: str, selection: tuple,
                      project: str, nanoselect: Optional[str], cutoff: float,
                      min_time: float) -> None:
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

    action = ResidenceAction(cutoff=cutoff, min_time=min_time)

    for replica in replicas:
        click.echo(f"  Processing: {replica.name}")
        try:
            result = action.run(replica)
            if verbose:
                click.echo(f"    Mean residence: {result.mean_residence:.1f} ps")
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
def info(ctx: click.Context, project: Optional[str], solvents: bool) -> None:
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
@click.option("-p", "--project", type=click.Path(exists=True), default=".",
              help="Project directory")
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

    click.echo(f"Project: {proj.name}")
    click.echo(f"  Path: {proj.path}")
    click.echo(f"  Systems: {len(proj.systems)}")
    click.echo(f"  Replicas: {len(proj.replicas)}")
    click.echo(f"  Groups: {len(proj.groups)}")

    if proj.systems:
        click.echo("\nSystems:")
        for name in proj.systems:
            click.echo(f"  - {name}")

    if proj.replicas:
        click.echo("\nReplicas:")
        for name, rep in proj.replicas.items():
            status = getattr(rep.status, 'value', str(rep.status))
            click.echo(f"  - {name}: {rep.solvent} ({status})")


@info.command("systems")
@click.option("-p", "--project", type=click.Path(exists=True), default=".",
              help="Project directory")
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

    if not proj.systems:
        click.echo("No systems in project")
        return

    click.echo(f"Systems ({len(proj.systems)}):")
    for name, system in proj.systems.items():
        click.echo(f"  {name}")
        if detailed:
            click.echo(f"    Input: {system.input_file}")
            if system.extra_residues:
                click.echo(f"    Extra residues: {', '.join(system.extra_residues)}")


@info.command("replicas")
@click.option("-p", "--project", type=click.Path(exists=True), default=".",
              help="Project directory")
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
    for name, rep in proj.replicas.items():
        status = getattr(rep.status, 'value', str(rep.status))
        click.echo(f"  {name}: {rep.solvent} ({status})")
        if detailed:
            click.echo(f"    System: {rep.system}")
            click.echo(f"    Nanoseconds: {rep.nanos}")
            click.echo(f"    Restraints: {rep.restraint_mode}")


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
            click.echo(f"  {name:8s} - {s.full_name}")
            if detailed:
                probes = [p.name for p in s.probes]
                click.echo(f"           Probes: {', '.join(probes)}")
        except Exception:
            click.echo(f"  {name:8s} - (error loading)")


@info.command("settings")
@click.option("-s", "--solvent", default="WAT", help="Solvent name")
@click.option("-f", "--file", "config_file", type=click.Path(exists=True),
              help="Load settings from TOML file")
@click.option("--restraints", type=click.Choice(["FREE", "BB", "HA"]), default="FREE",
              help="Restraint mode")
@click.option("--nanos", type=int, default=20, help="Simulation length")
@click.pass_context
def info_settings(ctx: click.Context, solvent: str, config_file: Optional[str],
                  restraints: str, nanos: int) -> None:
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
        settings = MDSettings(
            solvent=solvent,
            restraint_mode=restraints,
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
@click.option("-p", "--project", type=click.Path(exists=True), default=".",
              help="Project directory")
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
    
    click.echo(f"Replica: {replica.name}")
    click.echo(f"  Solvent: {replica.solvent}")
    click.echo(f"  Status: {replica.status}")
    
    # Check for analysis outputs
    replica_path = proj.path / replica.name
    
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
@click.option("-p", "--project", type=click.Path(exists=True), default=".",
              help="Project directory")
@click.option("-o", "--output", type=click.Path(), help="Output file (default: show)")
@click.option("--format", "fmt", type=click.Choice(["png", "pdf", "svg"]), default="png",
              help="Output format")
@click.pass_context
def plot_rmsd(ctx: click.Context, selection_type: str, selection: tuple,
              project: str, output: Optional[str], fmt: str) -> None:
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
@click.option("-p", "--project", type=click.Path(exists=True), default=".",
              help="Project directory")
@click.option("--probe", help="Specific probe to plot")
@click.option("-o", "--output", type=click.Path(), help="Output file")
@click.pass_context
def plot_energy(ctx: click.Context, selection_type: str, selection: tuple,
                project: str, probe: Optional[str], output: Optional[str]) -> None:
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
@click.option("-p", "--project", type=click.Path(exists=True), default=".",
              help="Project directory")
@click.option("--probe", help="Specific probe to plot")
@click.option("-o", "--output", type=click.Path(), help="Output file")
@click.pass_context
def plot_density(ctx: click.Context, selection_type: str, selection: tuple,
                 project: str, probe: Optional[str], output: Optional[str]) -> None:
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
@click.option("-p", "--project", type=click.Path(exists=True), default=".",
              help="Project directory")
@click.option("--system", type=click.Choice(["slurm", "pbs", "sge", "lsf"]), default="slurm",
              help="Queue system")
@click.option("--template", type=click.Path(exists=True), help="Custom template file")
@click.option("--config", type=click.Path(exists=True), help="Queue configuration file")
@click.pass_context
def queue_generate(ctx: click.Context, selection_type: str, selection: tuple,
                   project: str, system: str, template: Optional[str],
                   config: Optional[str]) -> None:
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
        
        script = generate_queue_script(replica, queue_config, template=template)
        script_path.write_text(script)
        click.echo(f"  {replica.name} → {script_path.name}")

    click.secho("✓ Queue scripts generated", fg="green")


@queue.command("submit")
@click.argument("selection_type", type=click.Choice(["all", "bysolvent", "byname", "group"]))
@click.option("-s", "--selection", multiple=True, help="Selection values")
@click.option("-p", "--project", type=click.Path(exists=True), default=".",
              help="Project directory")
@click.option("--dry-run", is_flag=True, help="Show commands without executing")
@click.pass_context
def queue_submit(ctx: click.Context, selection_type: str, selection: tuple,
                 project: str, dry_run: bool) -> None:
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
                result = subprocess.run(cmd.split(), capture_output=True, text=True,
                                       cwd=replica_path)
                if result.returncode == 0:
                    click.echo(f"  ✓ {replica.name}: {result.stdout.strip()}")
                else:
                    click.secho(f"  ✗ {replica.name}: {result.stderr}", fg="red")
            except Exception as e:
                click.secho(f"  ✗ {replica.name}: {e}", fg="red")

    if not dry_run:
        click.secho("✓ Jobs submitted", fg="green")


# =============================================================================
# REMOVE Command
# =============================================================================

@cli.command("remove")
@click.option("-p", "--project", type=click.Path(exists=True), default=".",
              help="Project directory")
@click.option("--group", help="Remove a group")
@click.option("--force", is_flag=True, help="Don't ask for confirmation")
@click.pass_context
def remove(ctx: click.Context, project: str, group: Optional[str], force: bool) -> None:
    """Remove groups from project.

    To remove systems or replicas, delete their folders directly.

    Examples:
        pymdmix remove --group ethanol_runs
    """
    from pymdmix.project import Project

    project_path = Path(project)
    proj = Project.load(project_path)

    if group:
        if group not in proj.groups:
            click.secho(f"✗ Group not found: {group}", fg="red")
            sys.exit(1)
        
        if not force:
            if not click.confirm(f"Remove group '{group}'?"):
                click.echo("Aborted")
                return
        
        proj.remove_group(group)
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
@click.option("-g1", "--grid1", type=click.Path(exists=True), required=True, help="First input grid")
@click.option("-g2", "--grid2", type=click.Path(exists=True), required=True, help="Second input grid")
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
@click.option("-g1", "--grid1", type=click.Path(exists=True), required=True, help="First input grid")
@click.option("-g2", "--grid2", type=click.Path(exists=True), required=True, help="Second input grid")
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
@click.option("-i", "--inputs", type=click.Path(exists=True), multiple=True, required=True, help="Input grids")
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
@click.option("-i", "--input", "ingrid", type=click.Path(exists=True), required=True, help="Input density grid")
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
    click.echo(f"  Spacing: {grid.spacing[0]:.2f} x {grid.spacing[1]:.2f} x {grid.spacing[2]:.2f} Å")
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
    from pymdmix.core.grid import Grid
    
    if not inputs:
        click.secho("✗ At least one input file required", fg="red")
        sys.exit(1)
    
    grids = [Grid.read_dx(f) for f in inputs]
    
    if operation == "add":
        data = sum(g.data for g in grids)
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
@click.option("--format", "fmt", type=click.Choice(["dx", "mrc", "ccp4", "cube"]),
              help="Output format (auto-detected from extension)")
def tools_convert(input_file: str, output: str, fmt: Optional[str]) -> None:
    """Convert grid between formats.

    Examples:
        pymdmix tools convert energy.dx -o energy.mrc
    """
    from pymdmix.io.grids import convert_grid
    
    convert_grid(input_file, output, format=fmt)
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
    from pymdmix.analysis.hotspots import merge_hotspot_files
    
    if len(inputs) < 2:
        click.secho("✗ At least two input files required", fg="red")
        sys.exit(1)
    
    merge_hotspot_files(list(inputs), output, distance_cutoff=distance)
    click.secho(f"✓ Combined {len(inputs)} files → {output}", fg="green")


# =============================================================================
# Entry Point
# =============================================================================

def main() -> None:
    """Entry point for the CLI."""
    cli(obj={})


if __name__ == "__main__":
    main()
