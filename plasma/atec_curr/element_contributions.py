import numpy as np
from constants import *
from plasma_functions import c_sound, eta_spitzer
from atomic_rates import ionisation_rate, recombination_rate, line_radiation_coefficient, continuum_radiation_coefficient


def index_local(i_var, i_v, i_dof):

    return (i_var * NODES_PER_ELEM + i_v) * N_DOF_NODE + i_dof


def get_element_contributions(element, nodes,
                              H_gauss, H_s_gauss, H_ss_gauss, input_params):

    n_loc = NODES_PER_ELEM * N_DOF_NODE
    a_loc   = np.zeros((N_VAR * n_loc, N_VAR * n_loc))
    rhs_loc = np.zeros(N_VAR * n_loc)

    theta = input_params["time"]["theta"]
    zeta  = input_params["time"]["zeta"]
    tstep = input_params["time"]["tstep"]
    time  = input_params["time"]["t"]         # time at the start of the step

    # ohmic heating is switched on after t_ohmic_on
    t_ohmic_on = float(input_params["physics"]["t_ohmic_on"])
    ohmic_on = 1.0 if time >= t_ohmic_on else 0.0

    # physics: diffusion coefficients [m^2/s]
    D_perp   = float(input_params["physics"]["D_perp"])
    chi_perp = float(input_params["physics"]["chi_perp"])
    L_par    = float(input_params["physics"]["L_par"])
    gamma_sheath = float(input_params["physics"]["gamma_sheath"])
    n_tot        = float(input_params["physics"]["n_tot"])       # n + n_neutral [m^-3]

    # reference values: the solved variables are n/n_ref and T/T_ref
    n_ref = float(input_params["normalization"]["n_ref"])
    T_ref = float(input_params["normalization"]["T_ref"])

    # floor of n in the time derivative of the pressure equation (normalized)
    n_floor = float(input_params["physics"]["n_floor"]) / n_ref

    # minimum temperature used in the sound speed [eV]
    T_min_cs = float(input_params["physics"]["T_min_cs"])

    # no radiation losses (ADAS line + continuum) below this electron temperature [eV]
    T_min_rad = float(input_params["physics"]["T_min_rad"])

    # energy lost by the plasma per ionization [eV]
    E_ionization = float(input_params["physics"]["E_ionization"])

    # Get quantities at gaussian points
    x0_s      = np.zeros(N_GAUSS)
    values    = np.zeros((N_GAUSS, N_VAR))    # values[igauss, i_var]
    values_s  = np.zeros((N_GAUSS, N_VAR))    # derivatives with respect to the local coordinate

    sizes = element.sizes

    for igauss in range(N_GAUSS):
        H = H_gauss[igauss]
        H_s = H_s_gauss[igauss]

        for i_v in range(NODES_PER_ELEM):
            node = nodes[i_v]
            for i_dof in range(N_DOF_NODE):
                x0_s[igauss]             +=  H_s[i_v,i_dof] * sizes[i_v,i_dof] * node.coord[i_dof]
                values[igauss, I_N]      +=  H[i_v,i_dof]   * sizes[i_v,i_dof] * node.values[I_N,i_dof]
                values_s[igauss, I_N]    +=  H_s[i_v,i_dof] * sizes[i_v,i_dof] * node.values[I_N,i_dof]
                values[igauss, I_T]      +=  H[i_v,i_dof]   * sizes[i_v,i_dof] * node.values[I_T,i_dof]
                values_s[igauss, I_T]    +=  H_s[i_v,i_dof] * sizes[i_v,i_dof] * node.values[I_T,i_dof]


    for igauss in range(N_GAUSS):
        x_jac = x0_s[igauss]
        n0    = values[igauss, I_N]
        T0    = values[igauss, I_T]
        n0_x  = values_s[igauss, I_N] / x_jac
        T0_x  = values_s[igauss, I_T] / x_jac
        n0_corr = max(n0, n_floor)     # n multiplying dT/dt in the time derivative

        # physical quantities for the atomic / plasma functions
        ne_SI = n0 * n_ref           # electron density [m^-3]
        Te_eV = T0 * T_ref           # electron temperature [eV]

        # sound speed
        # (T used in cs is floored at T_min_cs, with zero derivative below it)
        Te_cs = max(Te_eV, T_min_cs)
        cs0, cs0_T = c_sound(Te_cs, Te_cs)
        if Te_eV < T_min_cs:
            cs0_T = 0.0
        cs0_T = cs0_T * T_ref                    # d cs / d(T/T_ref), unknown is the normalized T

        # Ionization and recombination rate coefficients [m^3/s] and their derivatives with respect
        # to the physical n_e and T_e, converted to derivatives with respect to the normalized n and T
        S_ion, S_ion_n, S_ion_T = ionisation_rate(ne_SI, Te_eV)
        S_rec, S_rec_n, S_rec_T = recombination_rate(ne_SI, Te_eV)
        S_ion_n, S_rec_n = S_ion_n * n_ref, S_rec_n * n_ref
        S_ion_T, S_rec_T = S_ion_T * T_ref, S_rec_T * T_ref

        # Neutral density from the constant total density, n_neutral = n_tot - n (never negative)
        n_neutral   = max(n_tot - ne_SI, 0.0)
        n_neutral_n = - n_ref if n_tot > ne_SI else 0.0     # d n_neutral / d(n/n_ref)

        # Particle source, normalized by n_ref: n_e*n_neutral*S_ion - n_e*n_i*alpha_rec  (n_i = n_e)
        src_n   = n0 * (n_neutral * S_ion - ne_SI * S_rec)
        src_n_n = (n_neutral * S_ion - 2.0 * ne_SI * S_rec) \
                  + n0 * (n_neutral_n * S_ion + n_neutral * S_ion_n - ne_SI * S_rec_n)
        src_n_T = n0 * (n_neutral * S_ion_T - ne_SI * S_rec_T)

        # Ionization energy loss in the pressure equation: E_ionization per ionization,
        # normalized by n_ref*T_ref like the pressure equation (ionization rate = n_e*n_neutral*S_ion)
        E_ion_norm = E_ionization / T_ref
        L_ion   = E_ion_norm * n0 * n_neutral * S_ion
        L_ion_n = E_ion_norm * (n_neutral * S_ion + n0 * (n_neutral_n * S_ion + n_neutral * S_ion_n))
        L_ion_T = E_ion_norm * n0 * n_neutral * S_ion_T

        # ADAS radiation, W*m^3 coefficients and their derivatives converted to the normalized n and T.
        # Powers are in W/m^3, normalized like the pressure equation (divided by e*n_ref*T_ref)
        c_power_norm = 1.0 / (E_CHARGE * n_ref * T_ref)

        PLT, PLT_n, PLT_T = line_radiation_coefficient(ne_SI, Te_eV)         # neutral line radiation
        PRB, PRB_n, PRB_T = continuum_radiation_coefficient(ne_SI, Te_eV)    # recombination + bremsstrahlung
        PLT_n, PRB_n = PLT_n * n_ref, PRB_n * n_ref
        PLT_T, PRB_T = PLT_T * T_ref, PRB_T * T_ref

        # The PRB radiation contains the ionization potential released in the recombination, which was
        # already removed from the electrons at the ionization (E_ionization): subtract it to avoid
        # double counting, P_rad = (L_PRB - S_rec*phi_ion) n_e n_i   (Eq. 18 of Phys. Plasmas 30, 042509)
        phi_ion = E_ionization * E_CHARGE      # [J]
        PRB   = PRB   - phi_ion * S_rec
        PRB_n = PRB_n - phi_ion * S_rec_n
        PRB_T = PRB_T - phi_ion * S_rec_T

        # line radiation: n_e * n_neutral * PLT
        P_line   = c_power_norm * ne_SI * n_neutral * PLT
        P_line_n = c_power_norm * (n_ref * n_neutral * PLT + ne_SI * n_neutral_n * PLT + ne_SI * n_neutral * PLT_n)
        P_line_T = c_power_norm * ne_SI * n_neutral * PLT_T

        # continuum radiation: n_e * n_i * PRB (n_i = n_e)
        P_cont   = c_power_norm * ne_SI**2 * PRB
        P_cont_n = c_power_norm * (2.0 * n_ref * ne_SI * PRB + ne_SI**2 * PRB_n)
        P_cont_T = c_power_norm * ne_SI**2 * PRB_T

        # radiation losses are switched off below T_min_rad
        if Te_eV < T_min_rad:
            P_line, P_line_n, P_line_T = 0.0, 0.0, 0.0
            P_cont, P_cont_n, P_cont_T = 0.0, 0.0, 0.0

        # Jsat and its derivatives with respect to the normalized n and T
        J_sat   = E_CHARGE * ne_SI * cs0
        J_sat_n = E_CHARGE * n_ref * cs0
        J_sat_T = E_CHARGE * ne_SI * cs0_T

        # Spitzer resistivity; eta_T converted to d eta / d(T/T_ref)
        eta, eta_T = eta_spitzer(Te_eV, ne_SI, Zeff=1)
        eta_T = eta_T * T_ref

        # Ohmic heating eta*J^2 [W/m^3], normalized like the pressure equation (divided by e*n_ref*T_ref)
        S_ohm   = ohmic_on * c_power_norm * eta * J_sat**2
        S_ohm_n = ohmic_on * c_power_norm * 2.0 * eta * J_sat * J_sat_n
        S_ohm_T = ohmic_on * c_power_norm * (eta_T * J_sat**2 + 2.0 * eta * J_sat * J_sat_T)

        wg = W_GAUSS[igauss]
        H = H_gauss[igauss]
        H_s = H_s_gauss[igauss]

        for i_v in range(NODES_PER_ELEM):
            for i_dof in range(N_DOF_NODE):
                index_n = index_local(I_N, i_v, i_dof)
                index_T = index_local(I_T, i_v, i_dof)

                v_test = H[i_v,i_dof] * sizes[i_v,i_dof]
                v_x    = H_s[i_v,i_dof] * sizes[i_v,i_dof] / x_jac

                # particle diffusion equation: dn/dt = d/dx( D_perp dn/dx )
                rhs_n   = - D_perp * v_x * n0_x
                rhs_n  += - v_test * n0 * 2.0 * cs0 / L_par
                rhs_n  += v_test * src_n

                rhs_loc[index_n] +=  rhs_n * wg * x_jac * tstep

                # pressure equation: 3/2 d(nT)/dt = d/dx( n chi_perp dT/dx ) - 2 gamma n T cs / L_par + eta J_sat^2
                rhs_T   = - chi_perp * n0 * v_x * T0_x
                rhs_T  += - v_test * 2.0 * gamma_sheath * n0 * T0 * cs0 / L_par
                rhs_T  += v_test * S_ohm
                rhs_T  += - v_test * L_ion
                rhs_T  += - v_test * (P_line + P_cont)

                rhs_loc[index_T] +=  rhs_T * wg * x_jac * tstep

                for j_v in range(NODES_PER_ELEM):
                    for j_dof in range(N_DOF_NODE):

                        index_n_kl = index_local(I_N, j_v, j_dof)
                        index_T_kl = index_local(I_T, j_v, j_dof)

                        bas = H[j_v,j_dof] * sizes[j_v,j_dof]
                        bas_x = H_s[j_v,j_dof] * sizes[j_v,j_dof] / x_jac

                        # particle diffusion equation
                        amat_n  = 0.0
                        amat_n += v_test * (1+zeta) * bas
                        amat_n += D_perp * v_x * bas_x * theta * tstep
                        amat_n += v_test * bas * 2.0 * cs0 / L_par * theta * tstep
                        amat_n += - v_test * bas * src_n_n * theta * tstep

                        a_loc[index_n,index_n_kl] +=  amat_n * wg * x_jac

                        # linearization of the sound speed in the particle sink: coupling to T
                        amat_nT  = 0.0
                        amat_nT += v_test * bas * 2.0 * n0 * cs0_T / L_par * theta * tstep
                        amat_nT += - v_test * bas * src_n_T * theta * tstep

                        a_loc[index_n,index_T_kl] +=  amat_nT * wg * x_jac

                        # pressure equation, linearized in (dn, dT): d(nT) = n dT + T dn
                        # coupling to T
                        amat_T  = 0.0
                        amat_T += v_test * 1.5 * (1+zeta) * n0_corr * bas
                        amat_T += chi_perp * n0 * v_x * bas_x * theta * tstep
                        amat_T += v_test * bas * 2.0 * gamma_sheath * n0 * (cs0 + T0 * cs0_T) / L_par * theta * tstep
                        amat_T += - v_test * bas * S_ohm_T * theta * tstep
                        amat_T += v_test * bas * L_ion_T * theta * tstep
                        amat_T += v_test * bas * (P_line_T + P_cont_T) * theta * tstep

                        a_loc[index_T,index_T_kl] +=  amat_T * wg * x_jac

                        # coupling to n
                        amat_Tn  = 0.0
                        amat_Tn += v_test * 1.5 * (1+zeta) * T0 * bas
                        amat_Tn += chi_perp * v_x * bas * T0_x * theta * tstep
                        amat_Tn += v_test * bas * 2.0 * gamma_sheath * T0 * cs0 / L_par * theta * tstep
                        amat_Tn += - v_test * bas * S_ohm_n * theta * tstep
                        amat_Tn += v_test * bas * L_ion_n * theta * tstep
                        amat_Tn += v_test * bas * (P_line_n + P_cont_n) * theta * tstep

                        a_loc[index_T,index_n_kl] +=  amat_Tn * wg * x_jac

    return a_loc, rhs_loc
