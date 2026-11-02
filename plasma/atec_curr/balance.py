import numpy as np
from constants import *
from plasma_functions import c_sound
from atomic_rates import ionisation_rate, recombination_rate
from element_contributions import get_element_contributions, index_local

# Particle balance of the domain, per unit area perpendicular to x [m^-2 s^-1] (inventories in [m^-2]):
#   d/dt int(n)  = in_left_n  + in_right_n  + ionization - recombination - sheath
#   d/dt int(nn) = in_left_nn + in_right_nn - ionization + recombination + recycling
# "in_*" are the diffusive fluxes entering the domain through each end (positive = into the domain).
# In equilibrium both residuals vanish, and with R_recycling = 1: in_left_n + in_left_nn = 0 (closed cold end).
#
# The flux through a Dirichlet end is the consistent flux: from the weak form of the equation of the boundary
# node (the row replaced by the Dirichlet condition), int(v du/dt) = rhs_row + flux_in, taking du/dt = 0.
# It is exact in equilibrium, much more accurate than -D du/dx at the end, and makes the balance close.
# The flux through a reflective end is zero (natural condition).
# I_sat is the integral of J_sat = e n cs along the profile [A/m], the main result.
BALANCE_KEYS = ["I_sat", "inventory_n", "inventory_nn",
                "in_left_n", "in_right_n", "in_left_nn", "in_right_nn",
                "ionization", "recombination", "sheath", "recycling",
                "residual_n", "residual_nn"]


def is_dirichlet(input_params, name, side):
    """Same logic as construct_matrix.py: Dirichlet unless the variable has no entry or the side is reflective"""

    bnd_conditions = input_params["bnd_conditions"]
    return name in bnd_conditions and side not in bnd_conditions[name].get("reflective", [])


def boundary_fluxes(grid, input_params, H_gauss, H_s_gauss, H_ss_gauss, time):
    """Consistent fluxes [m^-2 s^-1] into the domain through the left and right ends, for n and nn"""

    n_ref = float(input_params["normalization"]["n_ref"])

    # element contributions per unit time (tstep = 1) at this time
    params = dict(input_params)
    params["time"] = dict(input_params["time"], tstep=1.0, t=time)

    fluxes = {}
    for side, element, i_v in (("left", grid.elements[0], 0), ("right", grid.elements[-1], NODES_PER_ELEM - 1)):
        nodes = [grid.nodes[i] for i in element.vertex_glob_index]
        a_loc, rhs_loc = get_element_contributions(element, nodes, H_gauss, H_s_gauss, H_ss_gauss, params)

        for i_var, name in ((I_N, "n"), (I_NN, "nn")):
            if is_dirichlet(input_params, name, side):
                fluxes[f"in_{side}_{name}"] = - rhs_loc[index_local(i_var, i_v, 0)] * n_ref
            else:
                fluxes[f"in_{side}_{name}"] = 0.0

    return fluxes


def particle_balance(grid, input_params, H_gauss, H_s_gauss, H_ss_gauss, time):

    L_par       = float(input_params["physics"]["L_par"])
    R_recycling = float(input_params["physics"]["R_recycling"])
    T_min_cs    = float(input_params["physics"]["T_min_cs"])
    n_ref = float(input_params["normalization"]["n_ref"])
    T_ref = float(input_params["normalization"]["T_ref"])

    balance = dict.fromkeys(BALANCE_KEYS, 0.0)

    # volume integrals, same rates as in element_contributions.py
    for element in grid.elements:
        nodes = [grid.nodes[i] for i in element.vertex_glob_index]
        sizes = element.sizes

        for igauss in range(N_GAUSS):
            H = H_gauss[igauss]
            H_s = H_s_gauss[igauss]

            x_jac = 0.0
            n0, T0, nn0 = 0.0, 0.0, 0.0
            for i_v in range(NODES_PER_ELEM):
                node = nodes[i_v]
                for i_dof in range(N_DOF_NODE):
                    x_jac += H_s[i_v,i_dof] * sizes[i_v,i_dof] * node.coord[i_dof]
                    n0    += H[i_v,i_dof]   * sizes[i_v,i_dof] * node.values[I_N,i_dof]
                    T0    += H[i_v,i_dof]   * sizes[i_v,i_dof] * node.values[I_T,i_dof]
                    nn0   += H[i_v,i_dof]   * sizes[i_v,i_dof] * node.values[I_NN,i_dof]

            ne_SI = n0 * n_ref
            Te_eV = T0 * T_ref
            n_neutral = nn0 * n_ref

            Te_cs = max(Te_eV, T_min_cs)
            cs0, cs0_T = c_sound(Te_cs, Te_cs)
            S_ion = ionisation_rate(ne_SI, Te_eV)[0]
            S_rec = recombination_rate(ne_SI, Te_eV)[0]

            dx = W_GAUSS[igauss] * x_jac
            balance["inventory_n"]   += ne_SI * dx
            balance["inventory_nn"]  += n_neutral * dx
            balance["ionization"]    += ne_SI * n_neutral * S_ion * dx
            balance["recombination"] += ne_SI**2 * S_rec * dx
            balance["sheath"]        += 2.0 * ne_SI * cs0 / L_par * dx
            balance["I_sat"]         += E_CHARGE * ne_SI * cs0 * dx       # same J_sat as in element_contributions.py

    balance["recycling"] = R_recycling * balance["sheath"]

    # fluxes through the ends, positive into the domain
    balance.update(boundary_fluxes(grid, input_params, H_gauss, H_s_gauss, H_ss_gauss, time))

    # net rate of change of the inventories (zero in equilibrium)
    balance["residual_n"]  = balance["in_left_n"] + balance["in_right_n"] \
                             + balance["ionization"] - balance["recombination"] - balance["sheath"]
    balance["residual_nn"] = balance["in_left_nn"] + balance["in_right_nn"] \
                             - balance["ionization"] + balance["recombination"] + balance["recycling"]

    return balance


def print_balance(balance):

    b = balance
    print(f"  Integral of J_sat: {b['I_sat']:.4e} A/m")
    print(f"  Particle balance [m^-2 s^-1]: inventory n = {b['inventory_n']:.3e}, nn = {b['inventory_nn']:.3e} [m^-2]")
    print(f"    n : in_left {b['in_left_n']:+.3e}  in_right {b['in_right_n']:+.3e}  ion {b['ionization']:+.3e}"
          f"  rec {-b['recombination']:+.3e}  sheath {-b['sheath']:+.3e}  -> net {b['residual_n']:+.3e}")
    print(f"    nn: in_left {b['in_left_nn']:+.3e}  in_right {b['in_right_nn']:+.3e}  ion {-b['ionization']:+.3e}"
          f"  rec {b['recombination']:+.3e}  recyc {b['recycling']:+.3e}  -> net {b['residual_nn']:+.3e}")
