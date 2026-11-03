import numpy as np
from constants import *
from plasma_functions import c_sound, L_parallel


def initial_conditions(grid, input_params):
    """
    Exponentially decaying profile from the left (hot end) to the right boundary value, for n and T:
    f(x) = right + (left - right) * (exp(-d/lambda) - exp(-L/lambda)) / (1 - exp(-L/lambda)),  d = x - x_min
    i.e. the boundary values at both ends and decay length lambda.
    The decay lengths are estimated from the balance of the perpendicular diffusion with the sheath losses,
    with L_par at the hot end and the sound speed at a fraction of the temperature of the hot end:
        n: D_perp   n'' = 2 n cs / L_par           -> lambda_n = sqrt(D_perp L_par / (2 cs))
        T: chi_perp T'' = 2 gamma_sheath T cs / L_par -> lambda_T = sqrt(chi_perp L_par / (2 gamma_sheath cs))
    Flat profile n_neutral_0 for the neutral density.
    """

    x_min = input_params["mesh"]["x_min"]
    x_max = input_params["mesh"]["x_max"]
    L0  = x_max - x_min

    D_perp   = float(input_params["physics"]["D_perp"])
    chi_perp = float(input_params["physics"]["chi_perp"])
    L_par    = L_parallel(x_min, float(input_params["physics"]["h_gap"]),
                          float(input_params["physics"]["theta_panel"]), float(input_params["physics"]["d_gap_min"]))
    gamma_sheath = float(input_params["physics"]["gamma_sheath"])

    # sound speed at half the temperature of the hot end
    T_est = 0.5 * float(input_params["bnd_conditions"]["T"]["left"])
    cs_est, cs_est_T = c_sound(T_est, 0.0)      # cold ions, as in the equations

    decay_length = {"n": np.sqrt(D_perp   * L_par / (2.0 * cs_est)),
                    "T": np.sqrt(chi_perp * L_par / (2.0 * gamma_sheath * cs_est))}
    print(f"Initial profiles: cs({T_est:.3g} eV) = {cs_est:.3e} m/s, decay lengths "
          f"lambda_n = {decay_length['n']*1e3:.3g} mm, lambda_T = {decay_length['T']*1e3:.3g} mm")

    for i_var, name in ((I_N, "n"), (I_T, "T")):
        left  = float(input_params["bnd_conditions"][name]["left"])
        right = float(input_params["bnd_conditions"][name]["right"])

        # normalized variables
        ref   = input_params["normalization"][REF_NAMES[name]]
        left  = left / ref
        right = right / ref

        lam = decay_length[name]
        e_L = np.exp(-L0 / lam)

        for node in grid.nodes:
            x = node.coord[0]
            x_s = node.coord[1]
            e_x = np.exp(-(x - x_min) / lam)
            node.values[i_var, 0] = right + (left - right) * (e_x - e_L) / (1.0 - e_L)
            node.values[i_var, 1] = - (left - right) * e_x / (lam * (1.0 - e_L)) * x_s

    # flat neutral density (normalized)
    nn_init = float(input_params["physics"]["n_neutral_0"]) / input_params["normalization"][REF_NAMES["nn"]]
    for node in grid.nodes:
        node.values[I_NN, 0] = nn_init
        node.values[I_NN, 1] = 0.0

    return grid
