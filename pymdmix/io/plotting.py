"""
Plotting utilities for pyMDMix.

Functions for visualizing:
- Energy grids (slices, projections)
- Convergence plots
- RMSD time series
- Hotspot distributions

Note:
    matplotlib is an optional dependency. Functions will raise ImportError
    if matplotlib is not available.
"""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    import matplotlib.pyplot as plt

    from pymdmix.core.grid import Grid

# Type aliases
ArrayLike = np.ndarray | Sequence


def _check_matplotlib():
    """Check if matplotlib is available."""
    try:
        import matplotlib.pyplot as plt

        return plt
    except ImportError:
        raise ImportError(
            "matplotlib is required for plotting. Install with: pip install matplotlib"
        )


# =============================================================================
# Energy Grid Visualization
# =============================================================================


def plot_energy_grid(
    grid: Grid,
    axis: int = 2,
    level: float | None = None,
    cmap: str = "RdBu_r",
    vmin: float | None = None,
    vmax: float | None = None,
    colorbar: bool = True,
    title: str | None = None,
    output: str | Path | None = None,
    figsize: tuple[float, float] = (8, 6),
) -> plt.Figure:
    """Plot a 2D slice of an energy grid.

    Args:
        grid: Grid to plot
        axis: Axis perpendicular to slice (0=X, 1=Y, 2=Z)
        level: Position along axis (Å). If None, uses middle.
        cmap: Matplotlib colormap name
        vmin: Minimum value for color scale
        vmax: Maximum value for color scale
        colorbar: Whether to show colorbar
        title: Plot title
        output: If given, save figure to this path
        figsize: Figure size in inches

    Returns:
        matplotlib Figure object

    Examples:
        >>> from pymdmix.core.grid import Grid
        >>> grid = Grid.from_dx("energy.dx")
        >>> fig = plot_energy_grid(grid, axis=2, level=0.0)
        >>> fig.savefig("slice.png")
    """
    plt = _check_matplotlib()

    from pymdmix.core.grid import Grid

    if not isinstance(grid, Grid):
        raise TypeError(f"Expected Grid, got {type(grid)}")

    # Determine slice position
    if level is None:
        # Middle of the grid
        idx = grid.shape[axis] // 2
    else:
        # Convert Angstrom to index
        idx = int((level - grid.origin[axis]) / grid.spacing[axis])
        idx = max(0, min(idx, grid.shape[axis] - 1))

    # Extract slice
    if axis == 0:
        slice_data = grid.data[idx, :, :].T
        xlabel, ylabel = "Y (Å)", "Z (Å)"
        extent = [
            grid.origin[1],
            grid.origin[1] + grid.shape[1] * grid.spacing[1],
            grid.origin[2],
            grid.origin[2] + grid.shape[2] * grid.spacing[2],
        ]
        axis_label = "X"
    elif axis == 1:
        slice_data = grid.data[:, idx, :].T
        xlabel, ylabel = "X (Å)", "Z (Å)"
        extent = [
            grid.origin[0],
            grid.origin[0] + grid.shape[0] * grid.spacing[0],
            grid.origin[2],
            grid.origin[2] + grid.shape[2] * grid.spacing[2],
        ]
        axis_label = "Y"
    else:  # axis == 2
        slice_data = grid.data[:, :, idx].T
        xlabel, ylabel = "X (Å)", "Y (Å)"
        extent = [
            grid.origin[0],
            grid.origin[0] + grid.shape[0] * grid.spacing[0],
            grid.origin[1],
            grid.origin[1] + grid.shape[1] * grid.spacing[1],
        ]
        axis_label = "Z"

    # Create figure
    fig, ax = plt.subplots(figsize=figsize)

    # Plot
    im = ax.imshow(
        slice_data,
        extent=extent,
        origin="lower",
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        aspect="equal",
    )

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    if colorbar:
        cbar = fig.colorbar(im, ax=ax)
        cbar.set_label("Energy (kcal/mol)")

    if title:
        ax.set_title(title)
    else:
        pos = grid.origin[axis] + idx * grid.spacing[axis]
        ax.set_title(f"{axis_label} = {pos:.1f} Å")

    fig.tight_layout()

    if output:
        fig.savefig(output, dpi=150, bbox_inches="tight")

    return fig


def plot_energy_projection(
    grid: Grid,
    axis: int = 2,
    method: str = "min",
    cmap: str = "RdBu_r",
    vmin: float | None = None,
    vmax: float | None = None,
    colorbar: bool = True,
    title: str | None = None,
    output: str | Path | None = None,
    figsize: tuple[float, float] = (8, 6),
) -> plt.Figure:
    """Plot a 2D projection of an energy grid.

    Args:
        grid: Grid to plot
        axis: Axis to project along (0=X, 1=Y, 2=Z)
        method: Projection method ("min", "max", "mean", "sum")
        cmap: Matplotlib colormap name
        vmin: Minimum value for color scale
        vmax: Maximum value for color scale
        colorbar: Whether to show colorbar
        title: Plot title
        output: If given, save figure to this path
        figsize: Figure size in inches

    Returns:
        matplotlib Figure object
    """
    plt = _check_matplotlib()

    # Project data
    if method == "min":
        proj = np.min(grid.data, axis=axis)
    elif method == "max":
        proj = np.max(grid.data, axis=axis)
    elif method == "mean":
        proj = np.mean(grid.data, axis=axis)
    elif method == "sum":
        proj = np.sum(grid.data, axis=axis)
    else:
        raise ValueError(f"Unknown projection method: {method}")

    # Setup axes
    if axis == 0:
        proj = proj.T
        xlabel, ylabel = "Y (Å)", "Z (Å)"
        extent = [
            grid.origin[1],
            grid.origin[1] + grid.shape[1] * grid.spacing[1],
            grid.origin[2],
            grid.origin[2] + grid.shape[2] * grid.spacing[2],
        ]
    elif axis == 1:
        proj = proj.T
        xlabel, ylabel = "X (Å)", "Z (Å)"
        extent = [
            grid.origin[0],
            grid.origin[0] + grid.shape[0] * grid.spacing[0],
            grid.origin[2],
            grid.origin[2] + grid.shape[2] * grid.spacing[2],
        ]
    else:
        proj = proj.T
        xlabel, ylabel = "X (Å)", "Y (Å)"
        extent = [
            grid.origin[0],
            grid.origin[0] + grid.shape[0] * grid.spacing[0],
            grid.origin[1],
            grid.origin[1] + grid.shape[1] * grid.spacing[1],
        ]

    fig, ax = plt.subplots(figsize=figsize)

    im = ax.imshow(
        proj,
        extent=extent,
        origin="lower",
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        aspect="equal",
    )

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    if colorbar:
        cbar = fig.colorbar(im, ax=ax)
        cbar.set_label("Energy (kcal/mol)")

    if title:
        ax.set_title(title)
    else:
        ax.set_title(f"{method.capitalize()} projection along {'XYZ'[axis]}")

    fig.tight_layout()

    if output:
        fig.savefig(output, dpi=150, bbox_inches="tight")

    return fig


# =============================================================================
# Convergence Analysis
# =============================================================================


def plot_convergence(
    values: ArrayLike,
    time: ArrayLike | None = None,
    window: int = 100,
    xlabel: str = "Frame",
    ylabel: str = "Value",
    title: str = "Convergence Analysis",
    output: str | Path | None = None,
    figsize: tuple[float, float] = (10, 4),
) -> plt.Figure:
    """Plot convergence of a property over time.

    Shows:
    - Raw values
    - Running average
    - Standard deviation bands

    Args:
        values: Array of property values
        time: Optional time axis values
        window: Rolling window size
        xlabel: X-axis label
        ylabel: Y-axis label
        title: Plot title
        output: If given, save figure to this path
        figsize: Figure size in inches

    Returns:
        matplotlib Figure object

    Examples:
        >>> energies = [...]  # Load from simulation
        >>> fig = plot_convergence(energies, ylabel="Energy (kcal/mol)")
    """
    plt = _check_matplotlib()

    values = np.asarray(values)
    n = len(values)

    if time is None:
        time = np.arange(n)
    else:
        time = np.asarray(time)

    # Calculate running statistics
    running_mean = np.zeros(n)
    running_std = np.zeros(n)

    for i in range(n):
        start = max(0, i - window + 1)
        window_data = values[start : i + 1]
        running_mean[i] = np.mean(window_data)
        running_std[i] = np.std(window_data)

    # Plot
    fig, ax = plt.subplots(figsize=figsize)

    # Raw values (light)
    ax.plot(time, values, alpha=0.3, color="gray", lw=0.5, label="Raw")

    # Running average
    ax.plot(time, running_mean, color="blue", lw=1.5, label=f"Running avg (n={window})")

    # Standard deviation bands
    ax.fill_between(
        time,
        running_mean - running_std,
        running_mean + running_std,
        alpha=0.2,
        color="blue",
        label="±1 σ",
    )

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend(loc="upper right")
    ax.grid(True, alpha=0.3)

    fig.tight_layout()

    if output:
        fig.savefig(output, dpi=150, bbox_inches="tight")

    return fig


def plot_rmsd(
    rmsd: ArrayLike,
    time: ArrayLike | None = None,
    labels: Sequence[str] | None = None,
    xlabel: str = "Time (ns)",
    ylabel: str = "RMSD (Å)",
    title: str = "RMSD",
    output: str | Path | None = None,
    figsize: tuple[float, float] = (10, 4),
) -> plt.Figure:
    """Plot RMSD time series.

    Args:
        rmsd: RMSD values (1D array or list of 1D arrays for multiple series)
        time: Time axis values
        labels: Labels for multiple series
        xlabel: X-axis label
        ylabel: Y-axis label
        title: Plot title
        output: If given, save figure to this path
        figsize: Figure size in inches

    Returns:
        matplotlib Figure object
    """
    plt = _check_matplotlib()

    # Handle single vs multiple series
    if isinstance(rmsd, np.ndarray) and rmsd.ndim == 1:
        rmsd = [rmsd]
    elif not isinstance(rmsd[0], (list, np.ndarray)):
        rmsd = [rmsd]

    fig, ax = plt.subplots(figsize=figsize)

    for i, r in enumerate(rmsd):
        r = np.asarray(r)
        t = time if time is not None else np.arange(len(r))
        label = labels[i] if labels else f"Series {i + 1}"
        ax.plot(t, r, label=label, lw=1)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.set_ylim(bottom=0)

    if len(rmsd) > 1:
        ax.legend()

    ax.grid(True, alpha=0.3)
    fig.tight_layout()

    if output:
        fig.savefig(output, dpi=150, bbox_inches="tight")

    return fig


# =============================================================================
# Hotspot Visualization
# =============================================================================


def plot_hotspot_energies(
    hotspots: Sequence[dict],
    xlabel: str = "Hotspot ID",
    ylabel: str = "Energy (kcal/mol)",
    title: str = "Hotspot Energies",
    output: str | Path | None = None,
    figsize: tuple[float, float] = (10, 4),
) -> plt.Figure:
    """Plot hotspot energies as a bar chart.

    Args:
        hotspots: List of hotspot dictionaries with 'id' and 'energy' keys
        xlabel: X-axis label
        ylabel: Y-axis label
        title: Plot title
        output: If given, save figure to this path
        figsize: Figure size in inches

    Returns:
        matplotlib Figure object
    """
    plt = _check_matplotlib()

    ids = [h.get("id", i) for i, h in enumerate(hotspots)]
    energies = [h.get("energy", 0.0) for h in hotspots]

    fig, ax = plt.subplots(figsize=figsize)

    colors = ["green" if e < 0 else "red" for e in energies]
    ax.bar(range(len(ids)), energies, color=colors, alpha=0.7)

    ax.set_xticks(range(len(ids)))
    ax.set_xticklabels(ids)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.axhline(0, color="black", lw=0.5)
    ax.grid(True, alpha=0.3, axis="y")

    fig.tight_layout()

    if output:
        fig.savefig(output, dpi=150, bbox_inches="tight")

    return fig


def plot_probe_distribution(
    probe_energies: dict,
    title: str = "Probe Energy Distribution",
    output: str | Path | None = None,
    figsize: tuple[float, float] = (10, 6),
) -> plt.Figure:
    """Plot distribution of energies for each probe type.

    Args:
        probe_energies: Dictionary mapping probe names to energy arrays
        title: Plot title
        output: If given, save figure to this path
        figsize: Figure size in inches

    Returns:
        matplotlib Figure object
    """
    plt = _check_matplotlib()

    n_probes = len(probe_energies)

    fig, axes = plt.subplots(1, n_probes, figsize=figsize, sharey=True)
    if n_probes == 1:
        axes = [axes]

    for ax, (name, energies) in zip(axes, probe_energies.items()):
        energies = np.asarray(energies)
        # Filter out masked values (e.g., 999)
        energies = energies[energies < 100]

        ax.hist(energies, bins=50, alpha=0.7, color="steelblue", edgecolor="white")
        ax.set_xlabel("Energy (kcal/mol)")
        ax.set_title(name)
        ax.axvline(0, color="red", ls="--", lw=1)

        # Statistics
        mean_e = np.mean(energies)
        _min_e = np.min(energies)
        ax.axvline(mean_e, color="orange", ls=":", lw=1, label=f"mean={mean_e:.2f}")

    axes[0].set_ylabel("Count")
    fig.suptitle(title)
    fig.tight_layout()

    if output:
        fig.savefig(output, dpi=150, bbox_inches="tight")

    return fig


def plot_replica_rmsd(
    replica,
    project_path: Path,
    output: str | Path | None = None,
    **kwargs,
) -> plt.Figure:
    """
    Plot RMSD for a replica from its alignment output files.

    Parameters
    ----------
    replica : Replica
        Replica object
    project_path : Path
        Path to project directory
    output : Path, optional
        Output file path

    Returns
    -------
    matplotlib.Figure
    """
    _check_matplotlib()

    replica_path = project_path / replica.name
    align_path = replica_path / "align"

    if not align_path.exists():
        raise FileNotFoundError(f"Alignment directory not found: {align_path}")

    # Find RMSD files
    rmsd_files = sorted(align_path.glob("*_bb_rmsd.dat")) + sorted(align_path.glob("*_bb_rmsd.out"))

    if not rmsd_files:
        raise FileNotFoundError(f"No RMSD files found in {align_path}")

    # Load RMSD data
    all_rmsd = []
    for f in rmsd_files:
        try:
            data = np.loadtxt(f, comments=["#", "@"])
            if data.ndim == 2:
                all_rmsd.append(data[:, 1])  # Second column is RMSD
            else:
                all_rmsd.append(data)
        except Exception:
            continue

    if not all_rmsd:
        raise ValueError(f"Could not load any RMSD data from {align_path}")

    # Concatenate
    rmsd = np.concatenate(all_rmsd)
    time = np.arange(len(rmsd)) * 0.001  # Assume 1 ps timestep, convert to ns

    return plot_rmsd(
        rmsd=rmsd,
        time=time,
        title=f"RMSD: {replica.name}",
        output=output,
        **kwargs,
    )
