#!/usr/bin/env python3
"""
Fit parameters for the model

    f(psi) = (rho_0 - rho_1)
             * (1 + rho_coef1*psi + rho_coef2*psi**2 + rho_coef3*psi**3)
             * (0.5 - 0.5 * tanh((psi - psi_barrier)/sig_n))
             + rho_1

from (psi, f) data points stored in a text/CSV file.

This script fits BOTH A := (rho_0 - rho_1) and rho_1 directly
(so you can also recover rho_0 = A + rho_1).

Usage examples
--------------
# basic usage (auto-detect delimiter, first two columns are psi and f)
python fit_psi_function.py data.txt

# specify columns (0-based), and save plot to custom path
python fit_psi_function.py data.csv --xcol 0 --ycol 2 --out fit.png

# file with a header row to skip
python fit_psi_function.py data.txt --skiprows 1

Dependencies: numpy, matplotlib, and scipy. (Install with: pip install numpy matplotlib scipy)
"""
from __future__ import annotations

import argparse
import math
import os
import sys
from dataclasses import dataclass
from typing import Tuple

import numpy as np
import matplotlib.pyplot as plt

try:
    from scipy.optimize import curve_fit
except ImportError as e:
    sys.stderr.write(
        "Error: scipy is required for nonlinear fitting. Install it with `pip install scipy`\n"
    )
    raise

# ----------------------------- Model --------------------------------------- #

def model(psi: np.ndarray,
          A: float,
          c1: float,
          c2: float,
          c3: float,
          barrier: float,
          sig: float,
          rho1: float) -> np.ndarray:
    """Model function.

    A = (rho_0 - rho_1)
    c1, c2, c3 are rho_coef(1:3)
    barrier = psi_barrier
    sig = sig_n (> 0)
    rho1 = rho_1 (offset / high-psi plateau)
    """
    poly = 1.0 + c1*psi + c2*psi**2 + c3*psi**3
    step = 0.5 - 0.5 * np.tanh((psi - barrier)/sig)
    return A * poly * step + rho1

# --------------------------- Utilities ------------------------------------- #

@dataclass
class FitResult:
    params: np.ndarray
    perr: np.ndarray
    cov: np.ndarray
    r2: float


def _sniff_delimiter(path: str) -> str | None:
    with open(path, 'r', encoding='utf-8', errors='ignore') as fh:
        for line in fh:
            s = line.strip()
            if not s or s.startswith('#'):
                continue
            # crude sniff
            if ';' in s:
                return ';'
            if ',' in s:
                return ','
            # default to whitespace
            return None
    return None


def read_xy(path: str,
            xcol: int = 0,
            ycol: int = 1,
            skiprows: int = 0,
            delimiter: str | None = None) -> Tuple[np.ndarray, np.ndarray]:
    """Read two columns (psi, f) from a text/CSV file.
    Falls back between whitespace and comma-separated formats.
    """
    if delimiter is None:
        delimiter = _sniff_delimiter(path)
    try:
        data = np.genfromtxt(
            path,
            delimiter=delimiter,
            comments='#',
            skip_header=skiprows,
            dtype=float,
        )
    except Exception as e:
        raise RuntimeError(f"Failed to read data from {path}: {e}")

    if data.ndim == 1:
        if data.size < max(xcol, ycol) + 1:
            raise ValueError("Data file appears to have too few columns.")
        x = np.asarray([data[xcol]], dtype=float)
        y = np.asarray([data[ycol]], dtype=float)
    else:
        if data.shape[1] <= max(xcol, ycol):
            raise ValueError(
                f"Requested columns xcol={xcol}, ycol={ycol} exceed file width {data.shape[1]}"
            )
        x = data[:, xcol].astype(float)
        y = data[:, ycol].astype(float)

    # Drop NaNs
    m = np.isfinite(x) & np.isfinite(y)
    x = x[m]
    y = y[m]
    if x.size < 7:
        raise ValueError("Need at least 7 valid data points to fit 7 parameters.")

    # sort by x for cleaner plotting/initialization
    order = np.argsort(x)
    return x[order], y[order]


def initial_guess(x: np.ndarray, y: np.ndarray) -> Tuple[np.ndarray, Tuple[np.ndarray, np.ndarray]]:
    """Heuristics to build a reasonable initial guess and bounds for curve_fit.

    Returns (p0, (lb, ub)). Now includes rho1.
    """
    n = x.size
    xr = x.max() - x.min()

    # Rough rho1 from high-psi tail (step ~ 0)
    k_tail = max(5, int(0.15 * n))
    tail_idx = np.argsort(x)[-k_tail:]
    rho1_0 = float(np.median(y[tail_idx]))

    # Use low-psi head where step ~ 1 to estimate A and polynomial coefficients
    head_idx = np.argsort(x)[:max(5, int(0.2 * n))]
    xh = x[head_idx]
    yh = y[head_idx] - rho1_0

    # Fit yh ≈ D0 + D1*x + D2*x^2 + D3*x^3 where Dk = A * ck, and c0=1 => D0 = A
    X = np.vstack([np.ones_like(xh), xh, xh**2, xh**3]).T
    try:
        D, *_ = np.linalg.lstsq(X, yh, rcond=None)
    except np.linalg.LinAlgError:
        D = np.array([np.median(yh), 0.0, 0.0, 0.0])

    A0 = float(D[0]) if abs(D[0]) > 1e-12 else (1.0 if np.allclose(y.mean(), 0) else y.mean())
    c10 = float(D[1]/A0) if abs(A0) > 1e-12 else 0.0
    c20 = float(D[2]/A0) if abs(A0) > 1e-12 else 0.0
    c30 = float(D[3]/A0) if abs(A0) > 1e-12 else 0.0

    # Two-plateau mid value for barrier guess
    y_low = float(np.median(y[head_idx]))
    y_high = float(np.median(y[tail_idx]))  # ~ rho1_0
    y_mid = 0.5 * (y_low + y_high)
    barrier0 = float(x[np.argmin((y - y_mid) ** 2)])

    # sig ~ fraction of x-range
    sig0 = max(1e-6, xr / 10.0)

    p0 = np.array([A0, c10, c20, c30, barrier0, sig0, rho1_0], dtype=float)

    # Bounds: sig > 0; barrier within a generous window; others free
    lb = np.array([-np.inf, -np.inf, -np.inf, -np.inf, x.min() - xr, 1e-12, -np.inf])
    ub = np.array([ np.inf,  np.inf,  np.inf,  np.inf, x.max() + xr,  np.inf,  np.inf])
    return p0, (lb, ub)


def fit_params(x: np.ndarray, y: np.ndarray) -> FitResult:
    p0, bounds = initial_guess(x, y)
    popt, pcov = curve_fit(
        model,
        x,
        y,
        p0=p0,
        bounds=bounds,
        maxfev=40000,
    )
    # Standard errors
    with np.errstate(invalid='ignore'):
        perr = np.sqrt(np.diag(pcov))

    # R^2
    yhat = model(x, *popt)
    ss_res = float(np.sum((y - yhat) ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan

    return FitResult(params=popt, perr=perr, cov=pcov, r2=r2)


def print_summary(res: FitResult) -> None:
    names = [
        "A = (rho_0 - rho_1)",
        "rho_coef1",
        "rho_coef2",
        "rho_coef3",
        "psi_barrier",
        "sig_n",
        "rho_1",
    ]
    print("\nFit results:")
    for name, val, err in zip(names, res.params, res.perr):
        if np.isfinite(err):
            print(f"  {name:>14s} = {val:.12g} ± {err:.3g}")
        else:
            print(f"  {name:>14s} = {val:.12g} (uncertainty unavailable)")

    A, *_, rho1 = res.params

    # Propagate uncertainty for rho_0 = A + rho_1
    try:
        varA = res.cov[0,0]
        varR1 = res.cov[-1,-1]
        covAR1 = res.cov[0,-1]
        var_rho0 = varA + varR1 + 2*covAR1
        err_rho0 = math.sqrt(var_rho0) if var_rho0 > 0 else float('nan')
    except Exception:
        err_rho0 = float('nan')

    rho0 = A + rho1
    if np.isfinite(err_rho0):
        print(f"\nDerived       rho_0 = A + rho_1 = {rho0:.12g} ± {err_rho0:.3g}")
    else:
        print(f"\nDerived       rho_0 = A + rho_1 = {rho0:.12g} (uncertainty unavailable)")

    print(f"\nGoodness of fit: R^2 = {res.r2:.6f}\n")


def plot_fit(x: np.ndarray,
             y: np.ndarray,
             res: FitResult,
             out_path: str | None = None,
             title: str | None = None) -> str:
    if out_path is None:
        out_path = "fit_psi.png"

    # Dense grid for smooth curve
    x_dense = np.linspace(x.min(), x.max(), 1000)
    y_dense = model(x_dense, *res.params)

    plt.figure(figsize=(7, 4.5))
    plt.scatter(x, y, s=16, label='data')
    plt.plot(x_dense, y_dense, linewidth=2, label='fit')
    plt.xlabel("psi_n")
    plt.ylabel("f(psi_n)")
    if title:
        plt.title(title)
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_path, dpi=160)
    try:
        plt.show()
    except Exception:
        pass
    finally:
        plt.close()
    return out_path


def main():
    ap = argparse.ArgumentParser(description="Fit f(psi_n) data to the provided model and plot the result.")
    ap.add_argument('infile', help='Path to input data file (two columns: psi, f).')
    ap.add_argument('--xcol', type=int, default=0, help='Zero-based column index for psi (default: 0).')
    ap.add_argument('--ycol', type=int, default=1, help='Zero-based column index for f(psi) (default: 1).')
    ap.add_argument('--skiprows', type=int, default=0, help='Number of header rows to skip (default: 0).')
    ap.add_argument('--delimiter', type=str, default=None, help='Field delimiter: ",", ";" or leave empty for auto/whitespace.')
    ap.add_argument('--out', type=str, default=None, help='Output plot path (default: <infile>_fit.png).')

    args = ap.parse_args()

    try:
        x, y = read_xy(args.infile, xcol=args.xcol, ycol=args.ycol, skiprows=args.skiprows, delimiter=args.delimiter)
    except Exception as e:
        sys.stderr.write(f"Failed to read input file: {e}\n")
        sys.exit(2)

    try:
        res = fit_params(x, y)
    except Exception as e:
        sys.stderr.write(f"Fitting failed: {e}\n")
        sys.exit(3)

    base, _ = os.path.splitext(os.path.basename(args.infile))
    out_path = args.out or f"{base}_fit.png"

    print_summary(res)

    try:
        saved = plot_fit(x, y, res, out_path=out_path, title=f"Fit for {base}")
        print(f"Plot saved to: {saved}")
    except Exception as e:
        sys.stderr.write(f"Plotting failed (continuing): {e}\n")


if __name__ == '__main__':
    main()
