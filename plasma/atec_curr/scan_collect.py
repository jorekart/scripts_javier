"""
Collect the results of a scan prepared with scan_launch.py: I_sat (integral of J_sat) at the last stored time
as a function of the scanned parameter, and the final profiles of n, T and nn of every run.
The parameter value of each run is read from the input stored in its .h5 file, so any other input parameter
can be used as x axis with --param.

Example:
    python scan_collect.py /scratch/users/$USER/atec_scans/D_n --logx
"""
import argparse
from pathlib import Path

import h5py
import numpy as np
import yaml
import matplotlib.pyplot as plt

from plot_1d_solver import VAR_LABELS


def read_param(f, path):
    """Value of the input parameter section.name stored in the solution file"""
    section, name = path.rsplit(".", 1)
    return f["input/" + section.replace(".", "/")].attrs[name]


def read_run(filename, param):
    with h5py.File(filename, "r") as f:
        balance = f["balance"]
        sheath = balance["sheath"][-1]
        run = {"value": float(read_param(f, param)),
               "t": f["time"][-1],
               "xvals": f["xvals"][:],
               "profiles": {name: f[name][-1, :] for name in VAR_LABELS if name in f},
               "I_sat": balance["I_sat"][-1],
               # convergence: change of I_sat over the last stored interval, and particle balance residual
               "dI_sat": abs(balance["I_sat"][-1] - balance["I_sat"][-2]) / abs(balance["I_sat"][-1])
                         if len(balance["I_sat"]) > 1 else np.nan,
               "residual": max(abs(balance["residual_n"][-1]), abs(balance["residual_nn"][-1])) / max(sheath, 1.0)}
    return run


def main():
    parser = argparse.ArgumentParser(description="Collect and plot the results of a parameter scan")
    parser.add_argument("scan_dir", help="scan directory written by scan_launch.py (contains scan.yaml)")
    parser.add_argument("--param", help="parameter used as x axis, section.name (default: the scanned one)")
    parser.add_argument("--logx", action="store_true", help="logarithmic x axis for I_sat")
    parser.add_argument("--log", action="store_true", help="logarithmic y axis for the densities")
    parser.add_argument("--no-show", action="store_true", help="only save the figures in scan_dir, without showing them")
    args = parser.parse_args()

    scan_dir = Path(args.scan_dir)
    with open(scan_dir / "scan.yaml") as f:
        scan_info = yaml.safe_load(f)
    param = args.param or scan_info["param"]
    name = param.split(".")[-1]

    runs = []
    for filename in scan_info["files"]:
        if not (scan_dir / filename).exists():
            print(f"missing {filename}, skipped")
            continue
        runs.append(read_run(scan_dir / filename, param))
    if not runs:
        raise SystemExit("no solution files found")
    runs.sort(key=lambda run: run["value"])

    print(f"{name:>12} {'t [ms]':>8} {'I_sat [A/m]':>12} {'dI_sat/I_sat':>13} {'residual':>9}")
    for run in runs:
        print(f"{run['value']:12.4g} {run['t']*1e3:8.3g} {run['I_sat']:12.4e} {run['dI_sat']:13.2e} {run['residual']:9.2e}")
    print("dI_sat/I_sat: relative change over the last stored interval; residual: particle balance / sheath loss")

    values = np.array([run["value"] for run in runs])

    # I_sat as a function of the parameter
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.plot(values, [run["I_sat"] for run in runs], "o-")
    ax.set_xlabel(name)
    ax.set_ylabel(r"$\int J_{sat}\,dx$ [A/m]")
    if args.logx:
        ax.set_xscale("log")
    ax.grid(True)
    fig.tight_layout()
    fig.savefig(scan_dir / f"I_sat_vs_{name}.png", dpi=150)

    # final profiles of every run
    fig, axes = plt.subplots(1, len(VAR_LABELS), figsize=(6 * len(VAR_LABELS), 5), sharex=True)
    colors = plt.cm.viridis(np.linspace(0, 1, len(runs)))
    for color, run in zip(colors, runs):
        for ax, var in zip(axes, VAR_LABELS):
            if var in run["profiles"]:
                ax.plot(run["xvals"], run["profiles"][var], color=color, label=f"{name} = {run['value']:g}")
    for ax, var in zip(axes, VAR_LABELS):
        ax.set_xlabel("x [m]")
        ax.set_ylabel(VAR_LABELS[var])
        if args.log and var != "T":
            ax.set_yscale("log")
        ax.grid(True)
    axes[0].legend()
    fig.suptitle(f"Final profiles, scan of {name}")
    fig.tight_layout()
    fig.savefig(scan_dir / f"profiles_vs_{name}.png", dpi=150)

    print(f"figures saved in {scan_dir}: I_sat_vs_{name}.png, profiles_vs_{name}.png")
    if not args.no_show:
        plt.show()


if __name__ == "__main__":
    main()
