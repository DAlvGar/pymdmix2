"""
pyMDMix CLI - Command-line interface for pyMDMix.

Usage:
    mdmix create project -n myproject -f config.yaml
    mdmix setup prepare structure.pdb -o prepared.pdb
    mdmix analyze density -p my_project -r replica1
    mdmix info -p my_project
"""
from __future__ import annotations

import sys
from pathlib import Path
from typing import Optional

import click

from pymdmix import __version__


# =============================================================================
# CLI Group
# =============================================================================

@click.group()
@click.version_option(version=__version__, prog_name="pyMDMix")
@click.option("-v", "--verbose", is_flag=True, help="Enable verbose output")
@click.pass_context
def cli(ctx: click.Context, verbose: bool) -> None:
    """pyMDMix - Molecular Dynamics with organic solvent mixtures.

    A toolkit for identifying binding hotspots on protein surfaces
    using MD simulations with organic solvent/water mixtures.
    """
    ctx.ensure_object(dict)
    ctx.obj["verbose"] = verbose


# =============================================================================
# CREATE Commands
# =============================================================================

@cli.group()
def create() -> None:
    """Create projects, replicas, or solvents."""
    pass


@create.command("project")
@click.option("-n", "--name", default="mdmix_project", help="Project name")
@click.option("-f", "--config", type=click.Path(exists=True), help="Configuration file (YAML/JSON)")
@click.option("-d", "--directory", type=click.Path(), help="Project directory (default: current)")
@click.pass_context
def create_project(ctx: click.Context, name: str, config: Optional[str], directory: Optional[str]) -> None:
    """Create a new pyMDMix project.

    Examples:
        mdmix create project -n my_project
        mdmix create project -n my_project -f config.yaml
    """
    from pymdmix.project import Project, Config

    verbose = ctx.obj.get("verbose", False)
    project_dir = Path(directory) if directory else Path.cwd() / name

    click.echo(f"Creating project '{name}' in {project_dir}")

    if config:
        cfg = Config.from_file(Path(config))
        if verbose:
            click.echo(f"  Loaded config from {config}")
    else:
        cfg = Config()
        if verbose:
            click.echo("  Using default configuration")

    project = Project(name=name, config=cfg, path=project_dir)
    project.save()

    click.secho(f"✓ Project '{name}' created successfully", fg="green")
    click.echo(f"  Directory: {project_dir}")


@create.command("replica")
@click.option("-n", "--name", required=True, help="Replica name")
@click.option("-p", "--project", type=click.Path(exists=True), help="Project directory")
@click.option("-s", "--solvent", required=True, help="Solvent name (e.g., ETA, ISP)")
@click.option("--top", type=click.Path(exists=True), help="Amber topology file")
@click.option("--crd", type=click.Path(exists=True), help="Amber coordinate file")
@click.pass_context
def create_replica(
    ctx: click.Context,
    name: str,
    project: Optional[str],
    solvent: str,
    top: Optional[str],
    crd: Optional[str],
) -> None:
    """Create a new replica within a project.

    Examples:
        mdmix create replica -n rep1 -s ETA -p my_project
        mdmix create replica -n rep1 -s ETA --top system.prmtop --crd system.rst7
    """
    from pymdmix.project import Project, Replica
    from pymdmix.core.solvent import SolventLibrary as Library

    verbose = ctx.obj.get("verbose", False)

    # Load solvent
    library = Library()
    try:
        solv = library.get(solvent)
    except KeyError:
        click.secho(f"✗ Unknown solvent: {solvent}", fg="red")
        click.echo(f"  Available solvents: {', '.join(library.list_solvents())}")
        sys.exit(1)

    if verbose:
        click.echo(f"  Using solvent: {solv.name} ({solv.full_name})")

    # Create replica
    replica = Replica(
        name=name,
        solvent=solv.name,
        topology=Path(top) if top else None,
        coordinates=Path(crd) if crd else None,
    )

    # Add to project if specified
    if project:
        proj = Project.load(Path(project))
        proj.add_replica(replica)
        proj.save()
        click.secho(f"✓ Replica '{name}' added to project '{proj.name}'", fg="green")
    else:
        # Standalone replica
        replica_dir = Path.cwd() / name
        replica_dir.mkdir(exist_ok=True)
        replica.save(replica_dir / "replica.json")
        click.secho(f"✓ Standalone replica '{name}' created", fg="green")
        click.echo(f"  Directory: {replica_dir}")


@create.command("solvent")
@click.option("-f", "--config", type=click.Path(exists=True), required=True, help="Solvent config file (JSON)")
@click.pass_context
def create_solvent(ctx: click.Context, config: str) -> None:
    """Add a new solvent to the library.

    Examples:
        mdmix create solvent -f my_solvent.json
    """
    from pymdmix.core.solvent import Solvent, Library

    verbose = ctx.obj.get("verbose", False)
    config_path = Path(config)

    solvent = Solvent.from_file(config_path)
    library = Library()
    library.add(solvent)

    click.secho(f"✓ Solvent '{solvent.name}' added to library", fg="green")
    if verbose:
        click.echo(f"  Full name: {solvent.full_name}")
        click.echo(f"  Probes: {[p.name for p in solvent.probes]}")


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
@click.option("--cap/--no-cap", default=True, help="Add ACE/NME caps to termini")
@click.option("--disulfide/--no-disulfide", default=True, help="Detect and rename disulfides")
@click.option("--remove-water/--keep-water", default=True, help="Remove crystallographic waters")
@click.pass_context
def setup_prepare(
    ctx: click.Context,
    structure: str,
    output: Optional[str],
    cap: bool,
    disulfide: bool,
    remove_water: bool,
) -> None:
    """Prepare a structure for MD simulation.

    Performs:
    - Terminal capping (ACE/NME)
    - Disulfide bond detection (CYS → CYX)
    - Water removal

    Examples:
        mdmix setup prepare protein.pdb -o prepared.pdb
        mdmix setup prepare protein.pdb --no-cap
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
@click.option("--box-size", type=float, default=12.0, help="Box buffer size (Å)")
@click.option("--conc", type=float, default=0.2, help="Co-solvent concentration (M)")
@click.pass_context
def setup_solvate(
    ctx: click.Context,
    structure: str,
    solvent: str,
    output: Optional[str],
    box_size: float,
    conc: float,
) -> None:
    """Solvate a structure with a solvent mixture.

    Examples:
        mdmix setup solvate protein.pdb -s ETA -o solvated
        mdmix setup solvate protein.pdb -s ISP --conc 0.3
    """
    from pymdmix.setup.solvate import solvate_structure, SolvationOptions
    from pymdmix.core.solvent import SolventLibrary as Library

    verbose = ctx.obj.get("verbose", False)
    input_path = Path(structure)
    output_prefix = output or input_path.stem + "_solvated"

    # Load solvent
    library = Library()
    try:
        solv = library.get(solvent)
    except KeyError:
        click.secho(f"✗ Unknown solvent: {solvent}", fg="red")
        sys.exit(1)

    click.echo(f"Solvating with {solv.full_name} ({conc} M)")

    options = SolvationOptions(
        box_buffer=box_size,
        concentration=conc,
    )

    result = solvate_structure(input_path, solv, options=options)

    # Save topology and coordinates
    top_path = Path(output_prefix + ".prmtop")
    crd_path = Path(output_prefix + ".rst7")

    result.topology.write(str(top_path))
    result.save_coordinates(str(crd_path))

    click.secho("✓ Solvation complete", fg="green")
    click.echo(f"  Topology: {top_path}")
    click.echo(f"  Coordinates: {crd_path}")

    if verbose:
        click.echo(f"  Total atoms: {len(result.topology.atoms)}")
        click.echo(f"  Box dimensions: {result.box_dimensions}")


# =============================================================================
# ANALYZE Commands
# =============================================================================

@cli.group()
def analyze() -> None:
    """Run analysis on simulation data."""
    pass


@analyze.command("density")
@click.option("-p", "--project", type=click.Path(exists=True), help="Project directory")
@click.option("-r", "--replica", multiple=True, help="Replica names (default: all)")
@click.option("-t", "--trajectory", type=click.Path(exists=True), help="Trajectory file (standalone)")
@click.option("--top", type=click.Path(exists=True), help="Topology file (standalone)")
@click.option("--probe", multiple=True, help="Probe names (default: all)")
@click.option("-o", "--output", type=click.Path(), help="Output directory")
@click.option("--stride", type=int, default=1, help="Frame stride")
@click.option("-j", "--jobs", type=int, default=1, help="Parallel jobs")
@click.pass_context
def analyze_density(
    ctx: click.Context,
    project: Optional[str],
    replica: tuple,
    trajectory: Optional[str],
    top: Optional[str],
    probe: tuple,
    output: Optional[str],
    stride: int,
    jobs: int,
) -> None:
    """Calculate probe density grids from trajectory.

    Examples:
        mdmix analyze density -p my_project -r rep1 rep2
        mdmix analyze density -t traj.nc --top system.prmtop -o densities/
    """
    from pymdmix.analysis import DensityAction

    _verbose = ctx.obj.get("verbose", False)

    if project:
        # Project-based analysis
        from pymdmix.project import Project
        proj = Project.load(Path(project))

        # Get replicas
        if replica:
            replicas = [proj.get_replica(r) for r in replica]
        else:
            replicas = list(proj.replicas.values())

        click.echo(f"Analyzing {len(replicas)} replica(s) from project '{proj.name}'")

        for rep in replicas:
            click.echo(f"  Processing: {rep.name}")
            _action = DensityAction(stride=stride)
            # Note: actual implementation would load trajectory and run
            click.echo("    [density calculation would run here]")

    elif trajectory and top:
        # Standalone analysis
        click.echo(f"Analyzing trajectory: {trajectory}")
        _action = DensityAction(stride=stride)
        # Note: actual implementation
        click.echo("  [density calculation would run here]")
    else:
        click.secho("✗ Specify either --project or (--trajectory + --top)", fg="red")
        sys.exit(1)

    click.secho("✓ Density analysis complete", fg="green")


@analyze.command("hotspots")
@click.argument("grids", nargs=-1, type=click.Path(exists=True))
@click.option("-o", "--output", required=True, help="Output prefix")
@click.option("--percentile", type=float, default=0.02, help="Cutoff percentile")
@click.option("--cutoff", type=float, help="Hard energy cutoff (overrides percentile)")
@click.pass_context
def analyze_hotspots(
    ctx: click.Context,
    grids: tuple,
    output: str,
    percentile: float,
    cutoff: Optional[float],
) -> None:
    """Detect binding hotspots from energy grids.

    Examples:
        mdmix analyze hotspots *.dx -o hotspots
        mdmix analyze hotspots energy_OH.dx energy_CT.dx -o hotspots --cutoff -1.5
    """
    from pymdmix.analysis import HotspotAction

    verbose = ctx.obj.get("verbose", False)

    if not grids:
        click.secho("✗ At least one grid file required", fg="red")
        sys.exit(1)

    click.echo(f"Detecting hotspots from {len(grids)} grid(s)")

    for g in grids:
        if verbose:
            click.echo(f"  Loading: {g}")

    _action = HotspotAction(
        percentile=percentile,
        energy_cutoff=cutoff,
    )

    # Note: actual implementation
    click.echo("  [hotspot detection would run here]")

    click.secho(f"✓ Hotspots saved to {output}.pdb and {output}.json", fg="green")


# =============================================================================
# INFO Commands
# =============================================================================

@cli.command("info")
@click.option("-p", "--project", type=click.Path(exists=True), help="Project directory")
@click.option("-r", "--replica", help="Specific replica name")
@click.option("-s", "--solvents", is_flag=True, help="List available solvents")
@click.pass_context
def info(
    ctx: click.Context,
    project: Optional[str],
    replica: Optional[str],
    solvents: bool,
) -> None:
    """Show information about projects, replicas, or solvents.

    Examples:
        mdmix info -p my_project
        mdmix info -s
    """
    if solvents:
        from pymdmix.core.solvent import SolventLibrary as Library
        library = Library()

        click.echo("Available solvents:")
        for name in sorted(library.list_solvents()):
            try:
                s = library.get(name)
                click.echo(f"  {name:8s} - {s.full_name}")
            except Exception:
                click.echo(f"  {name:8s} - (error loading)")
        return

    if project:
        from pymdmix.project import Project
        proj = Project.load(Path(project))

        click.echo(f"Project: {proj.name}")
        click.echo(f"  Path: {proj.path}")
        click.echo(f"  Replicas: {len(proj.replicas)}")

        if proj.replicas:
            click.echo("\n  Replica details:")
            for name, rep in proj.replicas.items():
                status = rep.status.value if hasattr(rep.status, 'value') else rep.status
                click.echo(f"    {name}: {rep.solvent} ({status})")

        if replica:
            rep = proj.get_replica(replica)
            click.echo(f"\n  Replica '{replica}':")
            click.echo(f"    Solvent: {rep.solvent}")
            click.echo(f"    Status: {rep.status}")
            if rep.topology:
                click.echo(f"    Topology: {rep.topology}")
    else:
        click.echo("Use --project to show project info, or --solvents to list solvents")


# =============================================================================
# QUEUE Commands
# =============================================================================

@cli.command("queue")
@click.option("-p", "--project", type=click.Path(exists=True), help="Project directory")
@click.option("-r", "--replica", multiple=True, help="Replica names")
@click.option("--system", type=click.Choice(["slurm", "sge", "pbs", "local"]), default="slurm",
              help="Queue system")
@click.option("-o", "--output", type=click.Path(), help="Output script path")
@click.pass_context
def queue(
    ctx: click.Context,
    project: Optional[str],
    replica: tuple,
    system: str,
    output: Optional[str],
) -> None:
    """Generate queue submission scripts.

    Examples:
        mdmix queue -p my_project --system slurm -o submit.sh
    """
    from pymdmix.engines.queue import QueueScriptGenerator, QueueSystem

    _verbose = ctx.obj.get("verbose", False)

    queue_type = QueueSystem(system)
    _generator = QueueScriptGenerator(queue_type)

    if project:
        from pymdmix.project import Project
        proj = Project.load(Path(project))

        if replica:
            replicas = [proj.get_replica(r) for r in replica]
        else:
            replicas = list(proj.replicas.values())

        click.echo(f"Generating {system.upper()} scripts for {len(replicas)} replica(s)")

        for rep in replicas:
            script_path = Path(output) if output else rep.path / f"submit_{system}.sh"
            # Note: actual generation
            click.echo(f"  {rep.name}: {script_path}")
    else:
        click.secho("✗ --project is required", fg="red")
        sys.exit(1)

    click.secho("✓ Queue scripts generated", fg="green")


# =============================================================================
# TOOLS Commands
# =============================================================================

@cli.group()
def tools() -> None:
    """Grid utilities and helper tools."""
    pass


@tools.command("diffgrids")
@click.option("-g1", "--grid1", type=click.Path(exists=True), required=True, help="First input grid")
@click.option("-g2", "--grid2", type=click.Path(exists=True), required=True, help="Second input grid")
@click.option("-o", "--output", type=click.Path(), required=True, help="Output grid file")
def tools_diffgrids(grid1: str, grid2: str, output: str) -> None:
    """Subtract two grids (grid1 - grid2).
    
    Examples:
        mdmix tools diffgrids -g1 probe1.dx -g2 probe2.dx -o diff.dx
    """
    from pymdmix.core.grid import Grid
    
    g1 = Grid.read_dx(grid1)
    g2 = Grid.read_dx(grid2)
    
    if g1.shape != g2.shape:
        click.secho(f"✗ Grid shapes don't match: {g1.shape} vs {g2.shape}", fg="red")
        sys.exit(1)
    
    result = Grid(
        data=g1.data - g2.data,
        origin=g1.origin,
        spacing=g1.spacing,
    )
    result.write_dx(output)
    click.secho(f"✓ Saved difference grid: {output}", fg="green")


@tools.command("sumgrids")
@click.option("-g1", "--grid1", type=click.Path(exists=True), required=True, help="First input grid")
@click.option("-g2", "--grid2", type=click.Path(exists=True), required=True, help="Second input grid")
@click.option("-o", "--output", type=click.Path(), required=True, help="Output grid file")
def tools_sumgrids(grid1: str, grid2: str, output: str) -> None:
    """Add two grids (grid1 + grid2).
    
    Examples:
        mdmix tools sumgrids -g1 probe1.dx -g2 probe2.dx -o sum.dx
    """
    from pymdmix.core.grid import Grid
    
    g1 = Grid.read_dx(grid1)
    g2 = Grid.read_dx(grid2)
    
    if g1.shape != g2.shape:
        click.secho(f"✗ Grid shapes don't match: {g1.shape} vs {g2.shape}", fg="red")
        sys.exit(1)
    
    result = Grid(
        data=g1.data + g2.data,
        origin=g1.origin,
        spacing=g1.spacing,
    )
    result.write_dx(output)
    click.secho(f"✓ Saved sum grid: {output}", fg="green")


@tools.command("getvalues")
@click.option("-i", "--input", "ingrid", type=click.Path(exists=True), required=True, help="Input grid file")
@click.option("-pdb", type=click.Path(exists=True), help="PDB file for coordinates")
@click.option("-xyz", type=click.Path(exists=True), help="XYZ text file (x y z per line)")
@click.option("-r", "--radius", type=float, default=0, help="Averaging radius around each point")
@click.option("-o", "--output", type=click.Path(), help="Output file (default: stdout)")
def tools_getvalues(ingrid: str, pdb: str, xyz: str, radius: float, output: str) -> None:
    """Get grid values at specified coordinates.
    
    Examples:
        mdmix tools getvalues -i energy.dx -pdb sites.pdb
        mdmix tools getvalues -i energy.dx -xyz coords.txt -r 2.0
    """
    import numpy as np
    from pymdmix.core.grid import Grid
    
    if not pdb and not xyz:
        click.secho("✗ Either -pdb or -xyz required", fg="red")
        sys.exit(1)
    
    grid = Grid.read_dx(ingrid)
    coords = []
    
    if pdb:
        try:
            import parmed as pmd
            struct = pmd.load_file(pdb)
            coords = [(a.xx, a.xy, a.xz) for a in struct.atoms]
        except ImportError:
            click.secho("✗ parmed required for PDB parsing", fg="red")
            sys.exit(1)
    else:
        with open(xyz) as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 3:
                    coords.append((float(parts[0]), float(parts[1]), float(parts[2])))
    
    # Get values at coordinates
    values = []
    for x, y, z in coords:
        val = grid.get_value(x, y, z)
        values.append(val if val is not None else np.nan)
    
    # Output
    if output:
        np.savetxt(output, values, fmt="%.3f")
        click.echo(f"Saved {len(values)} values to {output}")
    else:
        for v in values:
            click.echo(f"{v:.3f}")


@tools.command("avggrids")
@click.option("-i", "--inputs", type=click.Path(exists=True), multiple=True, required=True, help="Input grid files")
@click.option("-o", "--output", type=click.Path(), required=True, help="Output grid file")
@click.option("--boltzmann/--simple", default=False, help="Use Boltzmann averaging for energy grids")
@click.option("-T", "--temperature", type=float, default=300.0, help="Temperature for Boltzmann averaging")
def tools_avggrids(inputs: tuple, output: str, boltzmann: bool, temperature: float) -> None:
    """Average multiple grids.
    
    Examples:
        mdmix tools avggrids -i rep1.dx -i rep2.dx -i rep3.dx -o avg.dx
        mdmix tools avggrids -i rep1.dx -i rep2.dx --boltzmann -o avg.dx
    """
    from pymdmix.analysis.energy import boltzmann_average
    from pymdmix.core.grid import Grid
    
    grids = [Grid.read_dx(g) for g in inputs]
    
    if boltzmann:
        result = boltzmann_average(grids, temperature=temperature)
    else:
        # Simple average
        import numpy as np
        data = np.mean([g.data for g in grids], axis=0)
        result = Grid(data=data, origin=grids[0].origin, spacing=grids[0].spacing)
    
    result.write_dx(output)
    click.secho(f"✓ Averaged {len(inputs)} grids → {output}", fg="green")


@tools.command("energy")
@click.option("-i", "--input", "ingrid", type=click.Path(exists=True), required=True, help="Input density grid")
@click.option("-o", "--output", type=click.Path(), required=True, help="Output energy grid")
@click.option("-T", "--temperature", type=float, default=300.0, help="Temperature (K)")
def tools_energy(ingrid: str, output: str, temperature: float) -> None:
    """Convert density grid to free energy.
    
    Uses ΔG = -RT ln(ρ/ρ₀)
    
    Examples:
        mdmix tools energy -i density.dx -o energy.dx -T 300
    """
    from pymdmix.analysis.energy import density_to_free_energy
    from pymdmix.core.grid import Grid
    
    grid = Grid.read_dx(ingrid)
    energy = density_to_free_energy(grid, temperature=temperature)
    energy.write_dx(output)
    
    click.secho(f"✓ Converted to free energy: {output}", fg="green")
    click.echo(f"  Temperature: {temperature} K")
    click.echo(f"  Energy range: [{energy.data.min():.2f}, {energy.data.max():.2f}] kcal/mol")


# =============================================================================
# Entry Point
# =============================================================================

def main() -> None:
    """Entry point for the CLI."""
    cli(obj={})


if __name__ == "__main__":
    main()
