import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

from cherab.core.atomic import argon, neon, hydrogen, argon, tungsten, krypton, xenon
from cherab.openadas import OpenADAS

matplotlib.use("TkAgg")

# -------------------------
# ADAS helpers
# -------------------------
def get_rates_recombination(element, adas):
    # z = 1..Z
    return {i: adas.recombination_rate(element, int(i))
            for i in range(1, element.atomic_number + 1)}

def get_rates_ionisation(element, adas):
    # z = 0..Z-1
    return {i: adas.ionisation_rate(element, int(i))
            for i in range(0, element.atomic_number)}

def get_LT_radiation(element, adas):
    # z = 0..Z-1
    return {i: adas.line_radiated_power_rate(element, int(i))
            for i in range(0, element.atomic_number)}

def get_BR_radiation(element, adas):
    # z = 1..Z
    return {i: adas.continuum_radiated_power_rate(element, int(i))
            for i in range(1, element.atomic_number + 1)}

def eval_rate_over_time(ratefunc, ne, Te):
    """
    Try vector call first; fallback to python loop if ADAS callable doesn't accept arrays.
    """
    try:
        out = ratefunc(ne, Te)
        out = np.asarray(out)
        if out.shape == ne.shape:
            return out
    except Exception:
        pass

    return np.fromiter((ratefunc(nei, Tei) for nei, Tei in zip(ne, Te)),
                       dtype=float, count=len(ne))

# -------------------------
# fast linear interpolation for vectors vs time (called inside RHS)
# -------------------------
def lerp_time_vec(t, tgrid, Y):
    """
    Linear interpolation of Y over tgrid at scalar time t.
    Y shape: (N, Nt) returns (N,)
    """
    if t <= tgrid[0]:
        return Y[:, 0]
    if t >= tgrid[-1]:
        return Y[:, -1]

    i = np.searchsorted(tgrid, t) - 1
    t0 = tgrid[i]
    t1 = tgrid[i + 1]
    a = (t - t0) / (t1 - t0)
    return (1.0 - a) * Y[:, i] + a * Y[:, i + 1]

def lerp_time_scalar(t, tgrid, y):
    """
    Linear interpolation of y over tgrid at scalar time t.
    y shape: (Nt,) returns scalar
    """
    if t <= tgrid[0]:
        return float(y[0])
    if t >= tgrid[-1]:
        return float(y[-1])

    i = np.searchsorted(tgrid, t) - 1
    t0 = tgrid[i]
    t1 = tgrid[i + 1]
    a = (t - t0) / (t1 - t0)
    return float((1.0 - a) * y[i] + a * y[i + 1])

def interp_matrix_to_time(t_out, tgrid, Mgrid):
    """
    Interpolate matrix Mgrid (N, Ntgrid) onto t_out (Ntout,).
    Returns Mout (N, Ntout).
    N is small (charge states), so a loop over N is fast.
    """
    N = Mgrid.shape[0]
    Mout = np.empty((N, len(t_out)))
    for k in range(N):
        Mout[k, :] = np.interp(t_out, tgrid, Mgrid[k, :])
    return Mout

# -------------------------
# Main
# -------------------------
def main():
    # ---- Atomic + ADAS
    adas = OpenADAS(permit_extrapolation=True)
    elem = xenon  # change if you want
    Z = elem.atomic_number
    n_charges = Z + 1

    rates_ion   = get_rates_ionisation(elem, adas)       # z=0..Z-1
    rates_recom = get_rates_recombination(elem, adas)    # z=1..Z
    rates_line  = get_LT_radiation(elem, adas)           # z=0..Z-1
    rates_cont  = get_BR_radiation(elem, adas)           # z=1..Z

    # ---- Time base (your original grid for output)
    n_times = int(1e4)
    t_final = 10e-3
    time = np.linspace(0.0, t_final, n_times)

    # ---- Time-varying Te(t), ne(t)
    # Replace these with your real profiles (must be arrays length n_times)
    # Example: Te decays, ne ramps slightly
    Te    = np.zeros(n_times);   ne    = np.zeros(n_times)
    Te[:] = 0.5e3;                  ne[:] = 3.5e19        # Replace with functions of time

    # ---- Rate cache grid (coarser than time -> big speedup)
    # Tune n_rate: 500..5000 typical. Increase if Te/ne change very sharply.
    n_rate = 2000
    t_rate = np.linspace(time[0], time[-1], n_rate)
    Te_rate = np.interp(t_rate, time, Te)
    ne_rate = np.interp(t_rate, time, ne)
    ne_rate_clip = np.minimum(ne_rate, 2e21)

    # ---- Precompute ionisation/recombination arrays on t_rate
    # Sion[z] valid z=0..Z-1, Sion[Z]=0
    # Srec[z] valid z=1..Z,   Srec[0]=0
    Sion_rate = np.zeros((n_charges, n_rate))
    Srec_rate = np.zeros((n_charges, n_rate))

    for z in range(n_charges - 1):
        Sion_rate[z, :] = eval_rate_over_time(rates_ion[z], ne_rate_clip, Te_rate)
    for z in range(1, n_charges):
        Srec_rate[z, :] = eval_rate_over_time(rates_recom[z], ne_rate_clip, Te_rate)

    # ---- Precompute radiation power-rate arrays on t_rate: Prad[z,t]
    # z=0: line only, z=Z: continuum only, middle: line+continuum
    Prad_rate = np.zeros((n_charges, n_rate))

    # z = 0
    Prad_rate[0, :] = eval_rate_over_time(rates_line[0], ne_rate_clip, Te_rate)

    # z = Z (fully stripped) -> continuum rate exists for z=Z
    Prad_rate[-1, :] = eval_rate_over_time(rates_cont[Z], ne_rate_clip, Te_rate)

    # middle charges
    for z in range(1, n_charges - 1):
        line = eval_rate_over_time(rates_line[z], ne_rate_clip, Te_rate)   # z=1..Z-1 ok
        cont = eval_rate_over_time(rates_cont[z], ne_rate_clip, Te_rate)   # z=1..Z-1 ok
        Prad_rate[z, :] = line + cont

    # -------------------------
    # Stiff ODE solve: dn/dt = ne(t) * A(t) * n
    # -------------------------
    def rhs(t, nvec):
        ne_t = lerp_time_scalar(t, t_rate, ne_rate)
        sion = lerp_time_vec(t, t_rate, Sion_rate)
        srec = lerp_time_vec(t, t_rate, Srec_rate)

        dn = -nvec * (sion + srec)
        dn[:-1] += nvec[1:]  * srec[1:]     # recomb from upper
        dn[1:]  += nvec[:-1] * sion[:-1]    # ion from lower
        return ne_t * dn

    def jac(t, nvec):
        ne_t = lerp_time_scalar(t, t_rate, ne_rate)
        sion = lerp_time_vec(t, t_rate, Sion_rate)
        srec = lerp_time_vec(t, t_rate, Srec_rate)

        J = np.zeros((n_charges, n_charges))
        np.fill_diagonal(J, -(sion + srec))
        for z in range(n_charges - 1):
            J[z, z + 1] = srec[z + 1]
        for z in range(1, n_charges):
            J[z, z - 1] = sion[z - 1]
        return ne_t * J

    # initial condition: all neutrals
    y0 = np.zeros(n_charges)
    y0[0] = 1.0

    sol = solve_ivp(
        rhs,
        (time[0], time[-1]),
        y0,
        method="BDF",     # stiff
        jac=jac,          # helps speed/stability
        t_eval=time,      # output on your original time grid
        rtol=1e-6,
        atol=1e-12,
    )

    ion_i_density = sol.y  # shape (n_charges, n_times)

    # optional cleanup: positivity + normalize (since we're tracking fractions)
    ion_i_density[ion_i_density < 0] = 0.0
    ion_i_density /= np.maximum(ion_i_density.sum(axis=0, keepdims=True), 1e-300)

    # -------------------------
    # Coronal equilibrium distribution at each time (from cached rates)
    # -------------------------
    Sion_time = interp_matrix_to_time(time, t_rate, Sion_rate)
    Srec_time = interp_matrix_to_time(time, t_rate, Srec_rate)

    ratios = Sion_time[:-1, :] / np.maximum(Srec_time[1:, :], 1e-300)  # shape (Z, n_times)
    ion_i_density_cor = np.zeros((n_charges, n_times))
    ion_i_density_cor[0, :] = 1.0
    ion_i_density_cor[1:, :] = np.cumprod(ratios, axis=0)
    ion_i_density_cor /= np.maximum(ion_i_density_cor.sum(axis=0, keepdims=True), 1e-300)

    # -------------------------
    # Radiated power: Lrad(t) = sum_z Prad(z,t) * n_z(t)
    # -------------------------
    Prad_time = interp_matrix_to_time(time, t_rate, Prad_rate)

    Lrad     = (Prad_time * ion_i_density).sum(axis=0)
    Lrad_cor = (Prad_time * ion_i_density_cor).sum(axis=0)

    ratio = Lrad / np.maximum(Lrad_cor, 1e-300)

    # -------------------------
    # Save + plot
    # -------------------------
    tplot = time * 1e3  # ms

    output_file = (
        f"{elem.name}_Te{Te[0]*1e-3:.2f}keV_"
        f"ne{ne[0]*1e-20:.2f}_"
        f"time{t_final*1e3:.2f}ms"
    )

    np.savetxt(output_file, np.column_stack([tplot, ratio]))

    plt.plot(tplot, ratio, "-", color="black")
    plt.xlabel("Time [ms]")
    plt.ylabel("Radiation / Coronal rad.")
    plt.yscale("log")
    plt.grid(True)
    #plt.title(f"{elem.name}: stiff solve (BDF) + cached ADAS rates")
    plt.show()


if __name__ == "__main__":
    main()
