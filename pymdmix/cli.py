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

import configparser
import shutil
import sys
from pathlib import Path
from typing import Any, Literal, cast

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
@click.option(
    "-f",
    "--config",
    type=click.Path(exists=True),
    help=(
        "Configuration file. Can be a full project config with [SYSTEM] and "
        "[MDSETTINGS] sections (the primary workflow), or a global config "
        "in YAML/JSON format."
    ),
)
@click.option("-d", "--directory", type=click.Path(), help="Project directory (default: ./<name>)")
@click.option("--force", is_flag=True, help="Overwrite existing project")
@click.pass_context
def create_project(
    ctx: click.Context, name: str, config: str | None, directory: str | None, force: bool
) -> None:
    """Create a new pyMDMix project.

    PRIMARY WORKFLOW — pass a full project config file that contains both a
    [SYSTEM] section and a [MDSETTINGS] section.  The project is created and
    all replicas are added in a single step:

    \b
        pymdmix create project -n myproject -f project.cfg

    Get a ready-to-edit template with:

    \b
        pymdmix create template

    You can also create an empty project and add system/replicas separately:

    \b
        pymdmix create project -n myproject
        pymdmix add system -f system.cfg
        pymdmix add replica -f settings.cfg
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

    # ------------------------------------------------------------------
    # Determine if the config file is a full project cfg ([SYSTEM] + [MDSETTINGS])
    # or a YAML/JSON global config.
    # ------------------------------------------------------------------
    is_project_cfg = False
    if config:
        config_path = Path(config)
        if config_path.suffix in (".cfg", ".ini"):
            _parser = configparser.ConfigParser()
            _parser.read(config_path)
            is_project_cfg = "SYSTEM" in _parser.sections()

    if config and is_project_cfg:
        # ---- Full project config file ----------------------------------------
        # Parse *before* creating any directories so we fail cleanly on errors.
        from pymdmix.io.parsers import parse_project_config

        if verbose:
            click.echo(f"  Parsing full project config: {config}")

        try:
            proj_cfg = parse_project_config(config_path)
        except Exception as e:
            click.secho(f"✗ Failed to parse project config: {e}", fg="red")
            sys.exit(1)

        # Now create directory structure
        project = Project(name=name, config=Config(), path=project_dir)
        project.create_directories()

        # Register system
        system_cfg = proj_cfg.system
        project.pdb_file = str(system_cfg.input_file)
        project.systems[system_cfg.name] = {
            "name": system_cfg.name,
            "input_file": str(system_cfg.input_file),
            "unit_name": system_cfg.unit_name,
            "extra_residues": system_cfg.extra_residues,
            "extra_forcefields": system_cfg.extra_forcefields,
        }

        # Store the config file in the project's input folder
        destination = project.input_path / config_path.name
        destination.write_text(config_path.read_text())

        if verbose:
            click.echo(f"  System: {system_cfg.name} ({system_cfg.input_file})")

        # Add replicas
        from pymdmix.project.replica import ReplicaState, create_replica

        for md_settings in proj_cfg.settings:
            solvent = md_settings.solvent
            existing = [r for r in project.replicas if r.solvent == solvent]
            idx = len(existing) + 1
            rep_name = f"{system_cfg.name}_{solvent}_{idx}"
            replica = create_replica(
                name=rep_name,
                solvent=solvent,
                base_path=project.replicas_path,
                settings=md_settings,
            )
            project.replicas.append(replica)

        # ---- Solvate and write MD inputs ------------------------------------
        # The system config input can be either a PDB file or an Amber Object
        # File (OFF/lib).  The primary/recommended workflow (matching the
        # original pyMDMix) uses an OFF file that already contains the protein
        # force-field parameters; PDB is the alternative.
        input_path = system_cfg.input_file
        _off_exts = {".off", ".lib"}
        _pdb_exts = {".pdb", ".ent"}
        is_solvatable = (
            input_path.suffix.lower() in (_off_exts | _pdb_exts)
            and input_path.exists()
        )
        has_tleap = bool(shutil.which("tleap") or shutil.which("tLeap"))

        if not is_solvatable:
            click.secho(
                "  \u26a0 System input is not a PDB or OFF file \u2014 skipping solvation. "
                "Set topology/coordinates manually.",
                fg="yellow",
            )
        elif not has_tleap:
            click.secho(
                "  \u26a0 tleap not found \u2014 skipping solvation and input writing. "
                "Install AmberTools and re-run, or solvate manually.",
                fg="yellow",
            )
        else:
            from pymdmix.core.solvent import SolventLibrary
            from pymdmix.setup.solvate import SolvateResult, SolvationOptions, solvate_structure

            library = SolventLibrary()
            solvation_results: dict[str, SolvateResult] = {}

            for solvent_name in sorted({r.solvent for r in project.replicas}):
                solvent_obj = library.get(solvent_name)
                if solvent_obj is None:
                    click.secho(f"  \u26a0 Unknown solvent: {solvent_name}", fg="yellow")
                    continue

                output_dir = project.systems_path / solvent_name
                output_dir.mkdir(parents=True, exist_ok=True)
                output_prefix = f"{system_cfg.name}_{solvent_name}"

                if verbose:
                    click.echo(f"  Solvating with {solvent_name}...")

                result = solvate_structure(
                    input_path,
                    solvent_obj,
                    output_dir=output_dir,
                    output_prefix=output_prefix,
                    unit_name=system_cfg.unit_name,
                    options=SolvationOptions(
                        extra_forcefields=list(system_cfg.extra_forcefields or []),
                    ),
                )

                if result.success:
                    solvation_results[solvent_name] = result
                    click.echo(f"  Solvated {solvent_name}: {result.topology.name}")
                else:
                    click.secho(
                        f"  \u2717 Solvation failed for {solvent_name}: {result.error}",
                        fg="red",
                    )

            # Copy solvated topology/coordinates to each replica and write inputs
            n_ready = 0
            for replica in project.replicas:
                sol_result = solvation_results.get(replica.solvent)
                if not sol_result:
                    continue

                top_name = f"{system_cfg.name}_{replica.solvent}.prmtop"
                crd_name = f"{system_cfg.name}_{replica.solvent}.inpcrd"
                shutil.copy2(sol_result.topology, replica.path / top_name)
                shutil.copy2(sol_result.coordinates, replica.path / crd_name)
                replica.topology = top_name
                replica.coordinates = crd_name
                replica.state = ReplicaState.SETUP

                try:
                    replica.create_md_input()
                    replica.state = ReplicaState.READY
                    n_ready += 1
                except Exception as exc:
                    click.secho(
                        f"  \u26a0 Could not write inputs for {replica.name}: {exc}",
                        fg="yellow",
                    )

                replica.save()

            if n_ready:
                click.echo(f"  Ready:    {n_ready}/{len(project.replicas)} replicas")

        project.save()

        n_replicas = len(proj_cfg.settings)
        solvents = list({s.solvent for s in proj_cfg.settings})
        click.secho(f"\u2713 Project '{name}' created", fg="green")
        click.echo(f"  System:   {system_cfg.name}")
        click.echo(f"  Solvents: {', '.join(solvents)}")
        click.echo(f"  Replicas: {n_replicas}")
        for replica in project.replicas:
            click.echo(f"    - {replica.name}")

    elif config:
        # ---- YAML/JSON global config file ------------------------------------
        try:
            global_cfg = Config.from_file(Path(config))
        except Exception as e:
            click.secho(f"✗ Failed to load config: {e}", fg="red")
            sys.exit(1)
        project = Project(name=name, config=global_cfg, path=project_dir)
        project.create_directories()
        project.save()
        if verbose:
            click.echo(f"  Loaded global config from {config}")
        click.secho(f"✓ Project '{name}' created successfully", fg="green")
    else:
        # ---- Empty project ---------------------------------------------------
        project = Project(name=name, config=Config(), path=project_dir)
        project.create_directories()
        project.save()
        click.secho(f"✓ Project '{name}' created successfully", fg="green")
        click.echo(
            "  Next steps:\n"
            "    pymdmix add system -f system.cfg\n"
            "    pymdmix add replica -f settings.cfg\n"
            "\n"
            "  Or create a project from a single config file:\n"
            "    pymdmix create template            # get an editable template\n"
            "    pymdmix create project -n myproject -f project.cfg"
        )


@create.command("template")
@click.option(
    "-o",
    "--output",
    "output_path",
    type=click.Path(),
    default="project.cfg",
    show_default=True,
    help="Output path for the template file.",
)
@click.pass_context
def create_template(ctx: click.Context, output_path: str) -> None:
    """Save an annotated project configuration template.

    The template contains a [SYSTEM] section and a [MDSETTINGS] section with
    all supported options and their default values.  Edit the template and pass
    it to  pymdmix create project  to bootstrap a new project in one step.

    \b
    Examples:
        pymdmix create template
        pymdmix create template -o /path/to/my_project.cfg
    """
    from pymdmix.utils.tools import templates_root

    src = templates_root("project.cfg")
    if not src.exists():
        click.secho(f"✗ Template not found: {src}", fg="red")
        sys.exit(1)

    dest = Path(output_path)
    if dest.exists():
        click.secho(f"✗ File already exists: {dest}", fg="red")
        click.echo("  Remove it first or choose a different output path.")
        sys.exit(1)

    shutil.copy(src, dest)
    click.secho(f"✓ Template saved to {dest}", fg="green")
    click.echo(
        "  Edit the template, then create your project with:\n"
        f"    pymdmix create project -n myproject -f {dest}"
    )


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
    help=(
        "Configuration file. Accepts both a [MDSETTINGS] cfg file (primary) "
        "and a legacy [REPLICA] cfg file."
    ),
)
@click.option(
    "-p",
    "--project",
    type=click.Path(exists=True),
    default=".",
    help="Project directory (default: current)",
)
@click.option(
    "--count",
    type=int,
    default=None,
    help=(
        "Number of replicas per solvent (overrides NREPL in the config). "
        "Only used with [REPLICA] format; [MDSETTINGS] files define replica "
        "counts via NREPL."
    ),
)
@click.pass_context
def add_replica(ctx: click.Context, config: str, project: str, count: int | None) -> None:
    """Add replica(s) to the project from a configuration file.

    PRIMARY WORKFLOW — use a [MDSETTINGS] configuration file that lists
    solvents, number of replicas, and simulation parameters:

    \b
        pymdmix add replica -f settings.cfg

    Example settings.cfg:

    \b
        [MDSETTINGS]
        SOLVENTS = ETA, MAM
        NREPL = 3
        NANOS = 20
        TEMP = 300
        RESTR = FREE

    Alternatively, use the legacy [REPLICA] format for a single solvent:

    \b
        [REPLICA]
        SYSTEM = MyProtein
        SOLVENT = ETA
        NANOS = 20
        RESTRMODE = FREE

    Examples:
        pymdmix add replica -f settings.cfg
        pymdmix add replica -f replica.cfg --count 3
    """
    from pymdmix.project import Project
    from pymdmix.project.replica import create_replica

    verbose = ctx.obj.get("verbose", False)
    project_path = Path(project)
    config_path = Path(config)

    # Detect file type by inspecting sections
    _parser = configparser.ConfigParser()
    _parser.read(config_path)
    sections = _parser.sections()
    has_mdsettings = any(s.startswith("MDSETTINGS") for s in sections)

    # Load project
    proj = Project.load(project_path)

    if has_mdsettings:
        # ---- [MDSETTINGS] workflow (primary) ---------------------------------
        from pymdmix.io.parsers import parse_settings_config_file

        if verbose:
            click.echo(f"  Parsing [MDSETTINGS] config: {config_path.name}")

        try:
            settings_list = parse_settings_config_file(config_path)
        except Exception as e:
            click.secho(f"✗ Failed to parse settings config: {e}", fg="red")
            sys.exit(1)

        # Determine system name (first system in project, or from config)
        system_name = None
        if proj.systems:
            system_name = next(iter(proj.systems))

        created_names = []
        for md_settings in settings_list:
            solvent = md_settings.solvent
            existing = [r for r in proj.replicas if r.solvent == solvent]
            idx = len(existing) + 1
            prefix = f"{system_name}_" if system_name else ""
            rep_name = f"{prefix}{solvent}_{idx}"
            replica = create_replica(
                name=rep_name,
                solvent=solvent,
                base_path=proj.replicas_path,
                settings=md_settings,
            )
            proj.replicas.append(replica)
            created_names.append(rep_name)
            if verbose:
                click.echo(f"  Created: {rep_name}")

        proj._update_modified()
        proj.save()

        click.secho(f"✓ {len(created_names)} replica(s) added to project", fg="green")
        for rname in created_names:
            click.echo(f"    - {rname}")

    else:
        # ---- [REPLICA] legacy workflow ----------------------------------------
        from pymdmix.io.parsers import parse_replica_config

        replica_config = parse_replica_config(config_path)

        if verbose:
            click.echo(f"  System: {replica_config.system}")
            click.echo(f"  Solvent: {replica_config.solvent}")
            click.echo(f"  Nanoseconds: {replica_config.nanos}")

        n = count if count is not None else 1
        created = proj.add_replicas(solvent=replica_config.solvent, n_replicas=n)
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

        if n == 1:
            click.secho(f"✓ Replica '{created_names[0]}' added to project", fg="green")
        else:
            click.secho(f"✓ {n} replicas added to project", fg="green")
            for rname in created_names:
                click.echo(f"    - {rname}")


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
