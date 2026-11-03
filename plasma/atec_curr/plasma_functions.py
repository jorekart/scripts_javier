import numpy as np
from constants import *


def v_thermal_e(Te):
    Te_SI    = Te * E_CHARGE / KB
    return (KB*Te_SI/E_MASS)**0.5

def v_perp_e(Te):
    return v_thermal_e(Te) * np.sqrt(2)

def v_thermal_i(Ti):
    Ti_SI    = Ti * E_CHARGE / KB
    return (KB*Ti_SI/ION_MASS)**0.5

def L_parallel(x, h_gap, theta_panel, d_gap_min):
    """
    Parallel connection length [m] between neighbouring roof-top panels at depth x [m] in the shadow region:
    L_par = max( (h_gap - x) / tan(theta_panel), d_gap_min ),  theta_panel in degrees
    """
    return max((h_gap - x) / np.tan(np.radians(theta_panel)), d_gap_min)

def c_sound(Te,Ti): # Temperatures in eV
    cs = ((Te+Ti)*E_CHARGE/ION_MASS)**0.5
    cs_T = cs / (2.0 * (Te+Ti))         # d cs / dTe [m/s/eV], at fixed Ti
    return cs, cs_T

def v_perp_i(Ti):
    return v_thermal_i(Ti) * np.sqrt(2)

def v_alfven(Bt, ni, ION_MASS):
    return Bt/np.sqrt(ION_MASS*ni*mu0)

def ion_gyro_freq(Bt):
    return Zion * E_CHARGE * Bt / ION_MASS

def e_gyro_freq(Bt):
    return E_CHARGE * Bt / E_MASS

def e_gyro_radius(Te, Bt):
    return v_perp_e(Te)/e_gyro_freq(Bt)

def ion_gyro_radius(Te, Bt):
    return v_perp_i(Te)/ion_gyro_freq(Bt)

COULOMB_LOG_NE_MIN = 5.0e13   # [m^-3] lower limit of the density used in the Coulomb logarithm
COULOMB_LOG_TE_MIN = 0.5      # [eV]   lower limit of the temperature used in the Coulomb logarithm
COULOMB_LOG_MIN = 1.0         # limits of the value of the Coulomb logarithm
COULOMB_LOG_MAX = 25.0

def coulomb_log(Te, ne):
    # limit the arguments (the density can undershoot below zero in the solver)
    Te = np.maximum(Te, COULOMB_LOG_TE_MIN)
    ne = np.maximum(ne, COULOMB_LOG_NE_MIN)
    if np.ndim(Te)>0 or np.ndim(ne)>0:
        Te, ne = np.broadcast_arrays(np.asarray(Te, dtype=float), np.asarray(ne, dtype=float))
        mask1 = Te < 10
        mask2 = ~mask1
        Coulomb_Log = np.copy(Te)
        Coulomb_Log[mask1] = 23.0000 - np.log((ne[mask1]*1e-6)**0.5*Te[mask1]**(-1.5))
        Coulomb_Log[mask2] = 24.1513 - np.log((ne[mask2]*1e-6)**0.5*Te[mask2]**(-1.0))
    else:
        if (Te<10):
            Coulomb_Log = 23.0000 - np.log((ne*1e-6)**0.5*Te**(-1.5))
        else:
            Coulomb_Log = 24.1513 - np.log((ne*1e-6)**0.5*Te**(-1.0))
    Coulomb_Log = np.clip(Coulomb_Log, COULOMB_LOG_MIN, COULOMB_LOG_MAX)
    #return Coulomb_Log
    return 16.6

def coulomb_log_ee(Te, ne):
    
    Coulomb_Log = 23.5 - np.log((ne*1e-6)**0.5*Te**(-5/4)) - (1e-5 + (np.log(Te)-2)**2/16)**0.5
    
    return Coulomb_Log

def eta_spitzer(Te_in, ne, Zeff):
    Te = np.max([Te_in, 0.5])
    coef_Zeff = Zeff*(1.+1.198*Zeff+0.222*Zeff**2)/(1.+2.966*Zeff+0.753*Zeff**2) / ((1.+1.198+0.222)/(1.+2.966+0.753))
    eta = 1.65e-9 * coulomb_log(Te, ne) * (Te*0.001)**(-1.5) * coef_Zeff # From Wesson
    eta_T = -1.5 * 1.65e-9 * coulomb_log(Te, ne) * (Te*0.001)**(-2.5) * coef_Zeff * 1e-3
    if (Te_in <= 0.5):
        eta_T = 0.0
    return eta, eta_T

def collision_time_ei(Te, ne, Zeff, ni):
    # Wesson 4th edition page 65
    return 1.09e16 * (Te*0.001)**1.5 / (ni * Zeff**2 * coulomb_log(Te, ne))

def collision_time_ee(Te, ne, Zeff, ni):
    # Wesson 4th edition page 65
    # + https://farside.ph.utexas.edu/teaching/plasma/Plasma/node41.html
    # seems different by a factor 2
    return collision_time_ei(Te, ne, Zeff, ni) * 2.0

def collision_time_ii(Ti, ne, Zeff, ni):
    # Wesson 4th edition page 65
    return 6.60e17 * (Te*0.001)**1.5 / (ni * Zeff**4 * coulomb_log(Te, ne)) * (ION_MASS/E_MASS)**0.5

def e_collision_mean_free_path(Te, ne, Zeff, ni):
    return v_thermal_e(Te) * collision_time_ei(Te, ne, Zeff, ni)

def heat_diff_spitzer(Te, ne, Te_min):
    TeKeV = np.max([Te, Te_min]) * 1e-3 
    chi = 3.6e29 * TeKeV**(2.5) / ne 
    chi_T = 9.0e29 * TeKeV**(1.5) / ne * 1e-3 
    if (Te <= Te_min):
        chi_T = 0.0
    return chi, chi_T