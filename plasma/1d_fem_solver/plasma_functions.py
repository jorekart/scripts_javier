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

def coulomb_log(Te, ne):
    if np.ndim(Te)>0:
        mask1 = Te < 10
        mask2 = Te > 10
        Coulomb_Log = np.copy(Te)
        Coulomb_Log[mask1] = 23.0000 - np.log((ne*1e-6)**0.5*Te[mask1]**(-1.5))
        Coulomb_Log[mask2] = 24.1513 - np.log((ne*1e-6)**0.5*Te[mask2]**(-1.0))
    else:
        if (Te<10):
            Coulomb_Log = 23.0000 - np.log((ne*1e-6)**0.5*Te**(-1.5))
        else:
            Coulomb_Log = 24.1513 - np.log((ne*1e-6)**0.5*Te**(-1.0))
    return Coulomb_Log

def coulomb_log_ee(Te, ne):
    
    Coulomb_Log = 23.5 - np.log((ne*1e-6)**0.5*Te**(-5/4)) - (1e-5 + (np.log(Te)-2)**2/16)**0.5
    
    return Coulomb_Log

def eta_spitzer(Te_in, ne, Zeff):
    Te = np.max([Te_in, 1.0])
    coef_Zeff = Zeff*(1.+1.198*Zeff+0.222*Zeff**2)/(1.+2.966*Zeff+0.753*Zeff**2) / ((1.+1.198+0.222)/(1.+2.966+0.753))
    eta = 1.65e-9 * coulomb_log(Te, ne) * (Te*0.001)**(-1.5) * coef_Zeff # From Wesson
    eta_T = -1.5 * 1.65e-9 * coulomb_log(Te, ne) * (Te*0.001)**(-2.5) * coef_Zeff * 1e-3
    if (Te_in <= 1.0):
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