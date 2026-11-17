"""Integrated ohmic power vs Jphi for the L_parallel scan of scan_Lpar_Jphi.py.

P_ohm(t) = int eta(Te(x,t)) * Jphi^2 dx   [W/m^2, per unit flux-tube cross-section]

with eta the same Spitzer resistivity used by the solver (Zeff = 1).
The plot uses the last stored time of each run; the time traces are also
plotted to check how close each case is to steady state.

Usage:
    python plot_scan_Lpar_Jphi.py [--scan-dir DIR] [--out FILE.png]
"""
import argparse
import os

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import yaml

from constants import S_GAUSS, W_GAUSS
from plasma_functions import eta_spitzer
from scan_Lpar_Jphi import L_PAR, JPHI, case_name


def ohmic_power(run_dir):
    with open(os.path.join(run_dir, "input.yaml")) as f:
        params = yaml.safe_load(f)
    ne = float(params["physics"]["ne"])
    J = float(params["physics"]["Jphi"])

    with h5py.File(os.path.join(run_dir, "solution.h5"), "r") as f:
        x = f["xvals"][:]
        t = f["time"][:]
        Te = f["varvals"][:]

    # Integrate with the solver's own 4-point Gauss quadrature. On each element
    # x(s) and Te(s) are cubic in s, so the 4 stored sub-points (s = 0, 1/3, 2/3, 1,
    # n_sub = 4 in hdf5_export) define them exactly. A trapezoidal rule on the stored
    # points is NOT accurate here: it puts weight on the wall point (Te = 1 eV), where
    # eta ~ Te^-1.5 is ~1000x the core value, and overestimates P by up to ~3x.
    vander = np.vander(np.linspace(0.0, 1.0, 4), 4, increasing=True)
    cx = np.linalg.solve(vander, x.reshape(-1, 4).T)                  # (4, n_el)
    P = np.zeros(len(t))
    for it in range(len(t)):
        cT = np.linalg.solve(vander, Te[it].reshape(-1, 4).T)
        for s, w in zip(S_GAUSS, W_GAUSS):
            T_g = cT[0] + cT[1] * s + cT[2] * s**2 + cT[3] * s**3
            x_s = cx[1] + 2.0 * cx[2] * s + 3.0 * cx[3] * s**2
            eta = np.array([eta_spitzer(T, ne, 1.0)[0] for T in T_g])
            P[it] += np.sum(w * eta * J**2 * x_s)
    return t, P, x, Te


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--scan-dir", default="/scratch/users/artolaj/1d_fem_Lpar_Jphi_scan")
    parser.add_argument("--out", default=None)
    args = parser.parse_args()
    out = args.out or os.path.join(args.scan_dir, "P_ohm_vs_Jphi.png")

    fig, axes = plt.subplots(1, 3, figsize=(16, 4.8))
    table = []
    for L in L_PAR:
        P_final = []
        for J in JPHI:
            t, P, x, Te = ohmic_power(os.path.join(args.scan_dir, case_name(L, J)))
            P_final.append(P[-1])
            table.append((L, J, P[-1], Te[-1].max(), (P[-1] - P[-2]) / P[-1]))
            axes[1].plot(t * 1e3, P / 1e6, color=f"C{L_PAR.index(L)}",
                         alpha=0.4 + 0.6 * JPHI.index(J) / (len(JPHI) - 1))
            axes[2].plot(x / L, Te[-1], color=f"C{L_PAR.index(L)}",
                         alpha=0.4 + 0.6 * JPHI.index(J) / (len(JPHI) - 1))
        axes[0].plot(np.array(JPHI) / 1e6, np.array(P_final) / 1e6, "o-",
                     color=f"C{L_PAR.index(L)}", label=rf"$L_\parallel$ = {L:g} m")

    axes[0].set_xlabel(r"$J_\phi$ [MA/m$^2$]")
    axes[0].set_ylabel(r"$\int \eta J_\phi^2\, dx$ [MW/m$^2$]")
    axes[0].set_title("Integrated ohmic power (final time)")
    axes[0].legend()
    axes[1].set_yscale("log")
    axes[1].set_xlabel("t [ms]")
    axes[1].set_ylabel(r"$P_{ohm}$ [MW/m$^2$]")
    axes[1].set_title("Time traces (darker = higher $J_\\phi$)")
    axes[2].set_xlabel(r"$x / L_\parallel$")
    axes[2].set_ylabel(r"$T_e$ [eV]")
    axes[2].set_title("Final Te profiles")
    for ax in axes:
        ax.grid(True)
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    print(f"Saved {out}")

    print(f"{'L_par[m]':>9} {'Jphi[A/m2]':>11} {'P_ohm[MW/m2]':>13} {'Te_max[eV]':>11} {'dP/P last out':>14}")
    for L, J, P, Tmax, dP in table:
        print(f"{L:9g} {J:11.2e} {P/1e6:13.4g} {Tmax:11.4g} {dP:14.2%}")


if __name__ == "__main__":
    main()
