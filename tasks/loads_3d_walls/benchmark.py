import numpy as np

def heat_conductivity_w(T):
    T = np.asarray(T, dtype=float)
    k = 108.34 - 1.052e-2*T + 23419.9/T
    k = np.minimum(k, 173.0)
    k = np.maximum(k, 77.0)
    return k
  
def spec_heat_capacity_w(T):
    T = np.asarray(T, dtype=float)
    cp = 135.76 + 9.1159e-3*T + 2.3134e-9*T**3 - 6.5233e5 / T**2
    cp = np.maximum(cp, 129.0)
    cp = np.minimum(cp, 248.0)
    return cp

# --- Material properties / parameters ---
nx = 120
L = 0.012
stab_param = 0.01
qheat = 1e7
tduration = 1e-3

# Example values: replace with your actual ones
heat_conduct = heat_conductivity_w(T=1000)     # W/(m K)
spec_heat_cap = spec_heat_capacity_w(T=1000)    # J/(kg K)
mat_density = 19254.0     # kg/m^3
alpha = heat_conduct/(mat_density*spec_heat_cap)         # m^2/s

# --- Derived quantities ---
dx = L / float(nx - 1)
print("heat diffusivity =", alpha, "m^2/s")
print("Radial grid width for T =", dx)

dTdx = qheat / heat_conduct
dt = stab_param * dx**2 / alpha
nt = int(tduration/ dt)

T_tri_hist = np.zeros(nx)

print("nt =", nt)
print("dt =", dt)

# --- Time evolution ---
for i_time in range(nt):
    T_new = T_tri_hist.copy()

    T_new[0] = T_tri_hist[0] + 2.0 * stab_param * (
        T_tri_hist[1] - T_tri_hist[0] + dTdx * dx
    )

    for i in range(1, nx - 1):
        T_new[i] = T_tri_hist[i] + stab_param * (
            T_tri_hist[i + 1] - 2.0 * T_tri_hist[i] + T_tri_hist[i - 1]
        )

    T_tri_hist = T_new

# --- Final normalization / comparison value ---
denom = (2.0 / np.sqrt(3.1415)) * qheat * np.sqrt(tduration) / np.sqrt(
    heat_conduct * spec_heat_cap * mat_density
)

print(f'The relative error in % = {np.abs((T_tri_hist[0]-denom)) / denom*100}')