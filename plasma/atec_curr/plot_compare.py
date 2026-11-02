import argparse
import h5py
import matplotlib.pyplot as plt

from plot_1d_solver import VAR_LABELS


def plot_compare(filenames, labels=None, log=False, save=None):
    """
    Plot the profiles at the last stored time of several HDF5 solution files on the same axes,
    one subplot per variable (n, T, nn). The particle balance of each file at that time is printed.

    Parameters
    ----------
    filenames : list[str]
        HDF5 solution files.

    labels : list[str] or None, optional
        Legend label of each file, by default the file name.

    log : bool, optional
        Logarithmic y axis for the densities.

    save : str or None, optional
        If given, the figure is saved to this file instead of shown.
    """

    if labels is None:
        labels = filenames
    if len(labels) != len(filenames):
        raise ValueError(f"{len(filenames)} files but {len(labels)} labels")

    fig, axes = plt.subplots(1, len(VAR_LABELS), figsize=(6 * len(VAR_LABELS), 5), sharex=True)

    for filename, label in zip(filenames, labels):
        with h5py.File(filename, "r") as f:
            xvals = f["xvals"][:]
            t_last = f["time"][-1]
            profiles = {name: f[name][-1, :] for name in VAR_LABELS if name in f}
            balance = {key: f["balance"][key][-1] for key in f["balance"]} if "balance" in f else None

        for ax, name in zip(axes, VAR_LABELS):
            if name in profiles:
                ax.plot(xvals, profiles[name], label=f"{label} (t = {t_last*1e3:.3g} ms)")

        print(f"{label}: t = {t_last*1e3:.4g} ms")
        if balance is not None:
            if "I_sat" in balance:     # older files have no I_sat
                print(f"  integral of J_sat = {balance['I_sat']:.4e} A/m")
            print(f"  inventory n = {balance['inventory_n']:.3e}, nn = {balance['inventory_nn']:.3e} [m^-2]")
            print(f"  in_left n = {balance['in_left_n']:+.3e}, nn = {balance['in_left_nn']:+.3e},"
                  f" ionization = {balance['ionization']:.3e}, sheath = {balance['sheath']:.3e} [m^-2 s^-1]")
            print(f"  net n = {balance['residual_n']:+.3e}, nn = {balance['residual_nn']:+.3e} [m^-2 s^-1]")

    for ax, name in zip(axes, VAR_LABELS):
        ax.set_xlabel("x [m]")
        ax.set_ylabel(VAR_LABELS[name])
        if log and name != "T":
            ax.set_yscale("log")
        ax.grid(True)
    axes[0].legend()
    fig.suptitle("Final profiles")
    fig.tight_layout()

    if save:
        fig.savefig(save, dpi=150)
    else:
        plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compare the last stored profiles of several solution files")
    parser.add_argument("files", nargs="+", help="HDF5 solution files")
    parser.add_argument("--labels", nargs="+", help="legend label of each file")
    parser.add_argument("--log", action="store_true", help="logarithmic y axis for the densities")
    parser.add_argument("--save", help="save the figure to this file instead of showing it")
    args = parser.parse_args()

    plot_compare(args.files, labels=args.labels, log=args.log, save=args.save)
