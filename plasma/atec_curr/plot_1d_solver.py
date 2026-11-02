import sys
import h5py
import numpy as np
import matplotlib.pyplot as plt

VAR_LABELS = {"n": r"$n$ [m$^{-3}$]", "T": r"$T$ [eV]", "nn": r"$n_{neutral}$ [m$^{-3}$]"}


def plot_solution(
    filename,
    selected_times=None,
    n_slices=5,
    save=None,
):
    """
    Plot the n, T and neutral density profiles over time from an HDF5 solution file.

    Parameters
    ----------
    filename : str
        HDF5 solution file.

    selected_times : list[float] or None, optional
        Physical times [s] to plot. The closest stored
        timestep is used.

    n_slices : int, optional
        Number of equally spaced snapshots (in index
        space) to plot when selected_times is None.

    save : str or None, optional
        If given, the figure is saved to this file instead of shown.
    """

    with h5py.File(filename, "r") as f:
        xvals = f["xvals"][:]
        times = f["time"][:]
        profiles = {name: f[name][:] for name in VAR_LABELS if name in f}   # older files have no "nn"

    if selected_times is not None:
        indices = [np.argmin(np.abs(times - t)) for t in selected_times]
    else:
        indices = np.linspace(0, len(times) - 1, min(n_slices, len(times)), dtype=int)

    fig, axes = plt.subplots(1, len(profiles), figsize=(6 * len(profiles), 5), sharex=True)
    colors = plt.cm.viridis(np.linspace(0, 1, len(indices)))

    for ax, (name, values) in zip(axes, profiles.items()):
        for color, it in zip(colors, indices):
            ax.plot(xvals, values[it, :], color=color, label=f"t = {times[it]*1e3:.3g} ms")
        ax.set_xlabel("x [m]")
        ax.set_ylabel(VAR_LABELS[name])
        ax.grid(True)
    axes[0].legend()
    fig.suptitle("Solution evolution")
    fig.tight_layout()

    if save:
        fig.savefig(save, dpi=150)
    else:
        plt.show()


if __name__ == "__main__":
    plot_solution(sys.argv[1] if len(sys.argv) > 1 else "solution.h5", n_slices=6)
