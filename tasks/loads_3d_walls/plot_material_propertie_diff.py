from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt


# Functions for plotting in this notebook
def spec_heat_capacity_w(T):
    T = np.asarray(T, dtype=float)
    cp = 135.76 + 9.1159e-3*T + 2.3134e-9*T**3 - 6.5233e5 / T**2
    cp = np.maximum(cp, 129.0)
    cp = np.minimum(cp, 248.0)
    cp[T > T_melt_w] = np.nan
    return cp

def spec_heat_capacity_be(T):
    T = np.asarray(T, dtype=float)
    TC = T - 273.15
    cp = (
        -1.338948e-11*TC**4
        + 1.310600e-6*TC**3
        - 3.143247e-3*TC**2
        + 3.344647e0*TC
        + 1.741434e3
    )
    cp = np.where(TC > 1283.0, 3590.0, cp)
    cp = np.where(TC < 20.0, 1807.0, cp)
    cp[T > T_melt_be] = np.nan
    return cp

def heat_conductivity_w(T):
    T = np.asarray(T, dtype=float)
    k = 108.34 - 1.052e-2*T + 23419.9/T
    k = np.minimum(k, 173.0)
    k = np.maximum(k, 77.0)
    k[T > T_melt_w] = np.nan
    return k

def heat_conductivity_be(T):
    T = np.asarray(T, dtype=float)
    TC = T - 273.15
    k = (
        2.472945e-10*TC**4
        - 7.354535e-7*TC**3
        + 7.959136e-4*TC**2
        - 4.470410e-1*TC
        + 2.073906e2
    )
    k = np.minimum(k, 200.0)
    k = np.maximum(k, 60.0)
    k = np.where(TC > 1283.0, 84.7, k)
    k[T > T_melt_be] = np.nan
    return k

rho_w = 19254
rho_Be = 1750  # kg/m3

T0 = 273.15 + 25
T_melt_w = 3695.0
T_melt_be = 1560.0

T = np.linspace(T0, 4000, 1000)

# Plot 1: specific heat capacity
plt.figure(figsize=(6, 5))
plt.plot(T, heat_conductivity_w(T) / (spec_heat_capacity_w(T)*rho_w)*1e5, label="W", linewidth=2.5)
plt.plot(T, heat_conductivity_be(T) / (spec_heat_capacity_be(T)*rho_Be)*1e5, label="Be", linewidth=2.5)
plt.axvline(T_melt_w, linestyle="--", linewidth=2, alpha=0.8, label=r"W melt", color='lightblue')
plt.axvline(T_melt_be, linestyle="--", linewidth=2, alpha=0.8, label=r"Be melt", color='orange')
plt.xlabel("Temperature [K]", fontsize=16)
plt.ylabel(r"$\chi $ [$10^{-5}m^2/s$]$=K/(\rho c_p)$", fontsize=16)
plt.title("Heat diffusion coefficient vs T")
plt.grid(True, alpha=0.3)
plt.legend(frameon=True, fontsize=14)
#plt.yscale("log")
plt.tight_layout()
cp_path = "diffusion_W_Be.png"
plt.savefig(cp_path, dpi=200, bbox_inches="tight")
plt.show()
plt.close()

