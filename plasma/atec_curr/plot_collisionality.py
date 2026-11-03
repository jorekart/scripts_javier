"""
Collisionality check of the final profiles: electron mean free path lambda_e compared with the parallel
connection length L_par(x). The Spitzer resistivity (and the local ohmic heating eta J_sat^2) assumes
lambda_e << L_par; where lambda_e / L_par approaches 1 the electrons cross the gap almost without collisions.
lambda_e = v_te * tau_ei uses the same collision time and Coulomb logarithm as eta_spitzer (plasma_functions.py),
L_par(x) the geometry stored in the input of each solution file.
It also checks the single temperature assumption T_i = T_e: the electron-ion energy equilibration time
tau_eq = (m_i / 2 m_e) tau_ei (Wesson) must be short compared with the ion parallel residence time
tau_par = L_par / (2 c_s), the particle confinement time of the sheath sink in the model.

Example:
    python plot_collisionality.py sol_a.h5 sol_b.h5 --labels a b --threshold 0.1
"""
import argparse

import h5py
import numpy as np
import matplotlib.pyplot as plt

from constants import *
from plasma_functions import L_parallel, e_collision_mean_free_path, collision_time_ei, eta_spitzer, c_sound

T_MIN = 0.5       # [eV]   floor of T for the mean free path (same as the floor of eta_spitzer)
N_MIN = 1.0e16    # [m^-3] floor of n for the mean free path (n can undershoot to ~0 at the cold end)


def collisionality(filename):
    with h5py.File(filename, "r") as f:
        if "input" not in f:
            raise SystemExit(f"{filename} has no stored input (older file), cannot evaluate L_par")
        x = f["xvals"][:]
        t = f["time"][-1]
        n = f["n"][-1, :]
        T = f["T"][-1, :]
        physics = dict(f["input/physics"].attrs)

    L_par = np.array([L_parallel(xi, physics["h_gap"], physics["theta_panel"], physics["d_gap_min"]) for xi in x])

    T_c = np.maximum(T, T_MIN)
    n_c = np.maximum(n, N_MIN)
    lambda_e = e_collision_mean_free_path(T_c, n_c, 1.0, n_c)

    # sound speed as in element_contributions.py (T floored at T_min_cs)
    T_cs = np.maximum(T, physics["T_min_cs"])
    cs = c_sound(T_cs, 0.0)[0]       # cold ions

    # ohmic power density with J = J_sat
    J_sat = E_CHARGE * n * cs
    eta = np.array([eta_spitzer(Ti, ni, 1)[0] for ni, Ti in zip(n, T)])
    P_ohm = eta * J_sat**2

    # electron-ion energy equilibration time vs ion parallel residence time
    tau_eq  = ION_MASS / (2.0 * E_MASS) * collision_time_ei(T_c, n_c, 1.0, n_c)
    tau_par = L_par / (2.0 * cs)

    return {"x": x, "t": t, "L_par": L_par, "lambda_e": lambda_e, "ratio": lambda_e / L_par, "P_ohm": P_ohm,
            "tau_eq": tau_eq, "tau_par": tau_par, "ratio_eq": tau_eq / tau_par}


def main():
    parser = argparse.ArgumentParser(description="Electron mean free path / L_par of the final profiles")
    parser.add_argument("files", nargs="+", help="HDF5 solution files")
    parser.add_argument("--labels", nargs="+", help="legend label of each file")
    parser.add_argument("--threshold", type=float, default=0.1,
                        help="lambda_e/L_par above which the plasma is considered not collisional (default 0.1)")
    parser.add_argument("--save", help="also save the figure to this file")
    args = parser.parse_args()

    labels = args.labels or args.files
    if len(labels) != len(args.files):
        raise ValueError(f"{len(args.files)} files but {len(labels)} labels")

    fig, axes = plt.subplots(2, 2, figsize=(13, 10), sharex=True)
    axes = axes.flatten()

    for filename, label in zip(args.files, labels):
        c = collisionality(filename)
        line, = axes[0].semilogy(c["x"], c["lambda_e"], label=f"{label}: $\\lambda_e$")
        axes[0].semilogy(c["x"], c["L_par"], "--", color=line.get_color(), label=f"{label}: $L_\\parallel$")
        axes[1].semilogy(c["x"], c["ratio"], color=line.get_color(), label=label)
        axes[2].plot(c["x"], c["P_ohm"], color=line.get_color(), label=label)
        axes[3].semilogy(c["x"], c["ratio_eq"], color=line.get_color(), label=label)

        # fraction of the ohmic power (J = J_sat) deposited where lambda_e/L_par > threshold
        P_tot = np.trapz(c["P_ohm"], c["x"])
        P_col = np.trapz(np.where(c["ratio"] > args.threshold, c["P_ohm"], 0.0), c["x"])
        print(f"{label}: t = {c['t']*1e3:.4g} ms, max lambda_e/L_par = {c['ratio'].max():.3g}"
              f" at x = {c['x'][np.argmax(c['ratio'])]*1e3:.3g} mm;"
              f" ohmic power where lambda_e/L_par > {args.threshold:g}: {P_col/P_tot*100 if P_tot > 0 else 0:.3g}%")
        i_max = np.argmax(c["ratio_eq"])
        print(f"   max tau_eq/tau_par = {c['ratio_eq'][i_max]:.3g} at x = {c['x'][i_max]*1e3:.3g} mm"
              f" (tau_eq = {c['tau_eq'][i_max]:.3g} s, tau_par = {c['tau_par'][i_max]:.3g} s);"
              f" ohmic power where tau_eq/tau_par > 1:"
              f" {np.trapz(np.where(c['ratio_eq'] > 1, c['P_ohm'], 0.0), c['x'])/P_tot*100 if P_tot > 0 else 0:.3g}%")

    axes[0].set_ylabel("length [m]")
    axes[0].set_title(r"electron mean free path and $L_\parallel$")
    axes[1].axhline(1.0, color="k", lw=1)
    axes[1].axhline(args.threshold, color="k", lw=1, ls=":")
    axes[1].set_ylabel(r"$\lambda_e / L_\parallel$")
    axes[1].set_title(f"collisionality (dotted: threshold {args.threshold:g})")
    axes[2].set_ylabel(r"$\eta J_{sat}^2$ [W/m$^3$]")
    axes[2].set_title(r"ohmic power density with $J = J_{sat}$")
    axes[3].axhline(1.0, color="k", lw=1)
    axes[3].set_ylabel(r"$\tau_{eq} / \tau_\parallel$")
    axes[3].set_title(r"$T_i = T_e$ check: $\tau_{eq}=(m_i/2m_e)\tau_{ei}$, $\tau_\parallel = L_\parallel/2c_s$")
    for ax in axes:
        ax.set_xlabel("x [m]")
        ax.grid(True, which="both", alpha=0.4)
        ax.legend(fontsize=8)
    fig.tight_layout()

    if args.save:
        fig.savefig(args.save, dpi=150)
    plt.show()


if __name__ == "__main__":
    main()
