"""
Trajectory Alignment
====================

Align trajectory frames to a reference structure using RMSD fitting.
Uses MDAnalysis for alignment or cpptraj for Amber trajectories.

Examples
--------
>>> from pymdmix.analysis.align import align_trajectory
>>> align_trajectory(
...     topology="system.prmtop",
...     trajectory="prod.nc",
...     reference="reference.pdb",
...     output="aligned.nc",
...     mask="@CA",
... )
"""

from __future__ import annotations

import logging
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

import numpy as np

log = logging.getLogger(__name__)


@dataclass
class AlignmentResult:
    """Result of trajectory alignment."""

    output_trajectory: Path
    n_frames: int
    rmsd_mean: float
    rmsd_std: float
    method: str  # "mdanalysis" or "cpptraj"
    mean_rmsd: float = 0.0  # alias kept for back-compat; equals rmsd_mean

    def __repr__(self) -> str:
        return (
            f"AlignmentResult(frames={self.n_frames}, "
            f"RMSD={self.rmsd_mean:.3f}±{self.rmsd_std:.3f} Å)"
        )


def align_trajectory(
    topology: str | Path,
    trajectory: str | Path,
    output: str | Path,
    reference: str | Path | None = None,
    mask: str = "@CA,C,N",
    method: str = "auto",
    start: int | None = None,
    stop: int | None = None,
    step: int = 1,
) -> AlignmentResult:
    """
    Align trajectory frames to a reference structure.

    Parameters
    ----------
    topology : Path
        Topology file (prmtop, psf, pdb)
    trajectory : Path
        Input trajectory file
    output : Path
        Output aligned trajectory file
    reference : Path, optional
        Reference structure for alignment.
        If None, uses first frame of trajectory.
    mask : str
        Atom selection for RMSD fitting.
        Amber mask syntax: "@CA,C,N" (backbone), "@CA" (alpha carbons)
        MDAnalysis syntax: "name CA C N" or "protein and name CA"
    method : str
        Alignment method: "auto", "mdanalysis", or "cpptraj"
    start : int, optional
        First frame to process
    stop : int, optional
        Last frame to process
    step : int
        Frame stride

    Returns
    -------
    AlignmentResult
        Alignment statistics
    """
    topology = Path(topology)
    trajectory = Path(trajectory)
    output = Path(output)

    if reference:
        reference = Path(reference)

    # Auto-select method
    if method == "auto":
        # Prefer cpptraj for Amber files (faster, native)
        if trajectory.suffix in (".nc", ".mdcrd", ".crd"):
            method = "cpptraj"
        else:
            method = "mdanalysis"

    if method == "cpptraj":
        return _align_cpptraj(topology, trajectory, output, reference, mask, start, stop, step)
    else:
        return _align_mdanalysis(topology, trajectory, output, reference, mask, start, stop, step)


def _align_mdanalysis(
    topology: Path,
    trajectory: Path,
    output: Path,
    reference: Path | None,
    mask: str,
    start: int | None,
    stop: int | None,
    step: int,
) -> AlignmentResult:
    """Align using MDAnalysis."""
    try:
        import MDAnalysis as mda
        from MDAnalysis.analysis import align
    except ImportError:
        raise ImportError("MDAnalysis required for alignment. Install with: pip install MDAnalysis")

    # Convert Amber mask to MDAnalysis selection if needed
    selection = _convert_mask_to_mda(mask)

    # Load universe
    u = mda.Universe(str(topology), str(trajectory))

    # Reference structure
    if reference:
        ref = mda.Universe(str(reference))
    else:
        ref = u.copy()
        ref.trajectory[0]  # First frame as reference

    # Perform alignment
    log.info(f"Aligning {len(u.trajectory)} frames using MDAnalysis...")

    aligner = align.AlignTraj(
        u,
        ref,
        select=selection,
        filename=str(output),
        start=start,
        stop=stop,
        step=step,
    )
    aligner.run()

    # Calculate RMSD statistics
    rmsd_values = aligner.results.rmsd[:, 2]  # RMSD column

    return AlignmentResult(
        output_trajectory=output,
        n_frames=len(rmsd_values),
        rmsd_mean=float(np.mean(rmsd_values)),
        rmsd_std=float(np.std(rmsd_values)),
        method="mdanalysis",
        mean_rmsd=float(np.mean(rmsd_values)),
    )


def _align_cpptraj(
    topology: Path,
    trajectory: Path,
    output: Path,
    reference: Path | None,
    mask: str,
    start: int | None,
    stop: int | None,
    step: int,
) -> AlignmentResult:
    """Align using cpptraj (Amber)."""

    # Build cpptraj script
    script_lines = [
        f"parm {topology}",
        f"trajin {trajectory}",
    ]

    # Frame selection
    if start is not None or stop is not None or step != 1:
        frame_args = []
        if start is not None:
            frame_args.append(f"start {start}")
        if stop is not None:
            frame_args.append(f"stop {stop}")
        if step != 1:
            frame_args.append(f"offset {step}")
        # Modify trajin line
        script_lines[-1] = f"trajin {trajectory} {' '.join(frame_args)}"

    # Reference
    if reference:
        script_lines.append(f"reference {reference}")
        ref_name = "reference"
    else:
        script_lines.append("reference [first]")
        ref_name = "[first]"

    # Alignment
    script_lines.extend(
        [
            f"rms ToRef {mask} ref {ref_name}",
            f"trajout {output}",
            "run",
        ]
    )

    script = "\n".join(script_lines)

    # Write script to temp file
    with tempfile.NamedTemporaryFile(mode="w", suffix=".in", delete=False) as f:
        f.write(script)
        script_path = f.name

    try:
        # Run cpptraj
        log.info("Running cpptraj alignment...")
        result = subprocess.run(
            ["cpptraj", "-i", script_path],
            capture_output=True,
            text=True,
        )

        if result.returncode != 0:
            raise RuntimeError(f"cpptraj failed: {result.stderr}")

        # Parse RMSD from output
        rmsd_values = _parse_cpptraj_rmsd(result.stdout)

        return AlignmentResult(
            output_trajectory=output,
            n_frames=len(rmsd_values),
            rmsd_mean=float(np.mean(rmsd_values)) if rmsd_values else 0.0,
            rmsd_std=float(np.std(rmsd_values)) if rmsd_values else 0.0,
            method="cpptraj",
            mean_rmsd=float(np.mean(rmsd_values)) if rmsd_values else 0.0,
        )

    finally:
        Path(script_path).unlink(missing_ok=True)


def _convert_mask_to_mda(mask: str) -> str:
    """Convert Amber mask syntax to MDAnalysis selection."""
    # If already looks like MDAnalysis syntax, return as-is
    if " and " in mask or " or " in mask or mask.startswith("name "):
        return mask

    # Simple Amber mask conversion
    # @CA,C,N -> name CA C N
    if mask.startswith("@"):
        atoms = mask[1:].split(",")
        return "name " + " ".join(atoms)

    # :WAT -> resname WAT
    if mask.startswith(":"):
        resnames = mask[1:].split(",")
        return "resname " + " ".join(resnames)

    # Default: assume protein backbone
    if mask in ("backbone", "bb"):
        return "protein and name CA C N O"

    return mask


def _parse_cpptraj_rmsd(output: str) -> list[float]:
    """Parse RMSD values from cpptraj output."""
    rmsd_values = []

    for line in output.split("\n"):
        # Look for RMSD output lines
        if "RMSD" in line and "=" in line:
            try:
                # Extract RMSD value
                parts = line.split("=")
                if len(parts) >= 2:
                    value = float(parts[-1].strip().split()[0])
                    rmsd_values.append(value)
            except (ValueError, IndexError):
                continue

    return rmsd_values


def align_replica(
    replica,
    reference: str | Path | None = None,
    mask: str = "@CA,C,N",
    output_suffix: str = "_aligned",
    **kwargs,
) -> AlignmentResult:
    """
    Align all trajectories in a replica.

    Parameters
    ----------
    replica : Replica
        Replica object with trajectory paths
    reference : Path, optional
        Reference structure. If None, uses replica's reference PDB.
    mask : str
        Atom selection for fitting
    output_suffix : str
        Suffix for output trajectory files

    Returns
    -------
    AlignmentResult
        Combined alignment statistics
    """
    from pymdmix.project.replica import Replica

    if not isinstance(replica, Replica):
        raise TypeError(f"Expected Replica, got {type(replica)}")

    if not replica.topology:
        raise ValueError("Replica has no topology file")

    if not replica.trajectory:
        raise ValueError("Replica has no trajectory file")

    # Use replica's reference if not provided
    if reference is None:
        reference = replica.reference

    # Build output path
    traj_path = replica.path / replica.trajectory
    output = traj_path.with_stem(traj_path.stem + output_suffix)

    return align_trajectory(
        topology=replica.path / replica.topology,
        trajectory=traj_path,
        output=output,
        reference=reference,
        mask=mask,
        **kwargs,
    )


# =============================================================================
# AlignAction wrapper class
# =============================================================================


class AlignAction:
    """
    Action wrapper for trajectory alignment.

    Parameters
    ----------
    mask : str, optional
        Atom selection for RMSD fitting (default: "@CA,C,N")
    reference : Path, optional
        Reference structure for alignment
    nprocs : int
        Number of processors (for API compatibility)
    """

    name = "align"
    description = "Align trajectory to reference structure"

    def __init__(
        self,
        mask: str | None = None,
        reference: Path | None = None,
        nprocs: int = 1,
    ):
        self.mask = mask or "@CA,C,N"
        self.reference = reference
        self.nprocs = nprocs

    def run(self, replica, step_range: tuple | None = None, **kwargs) -> AlignmentResult:
        """
        Run alignment on a replica.

        Parameters
        ----------
        replica : Replica
            Replica to align
        step_range : tuple, optional
            (start, end) nanoseconds to align

        Returns
        -------
        AlignmentResult
        """
        result = align_replica(
            replica,
            reference=self.reference,
            mask=self.mask,
            **kwargs,
        )
        # Ensure the mean_rmsd alias is in sync (align_replica may not set it)
        result.mean_rmsd = result.rmsd_mean
        return result
