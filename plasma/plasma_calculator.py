import numpy as np

# Input parameters
Te   = 10    # Electron temperature in eV
Ti   = 10     # Ion temperature in eV
ne   = 1e20   # Electron density in 1e20 m^-3 units
ni   = ne      # Ion density in 1e20 m^-3 units
Bt   = 5.3    # Toroidal field
Zion = 1.0    # Charge on main ions
Zeff = 2    # Effective charge
atomic_mass = 2.014  #Atomic mass the atom (2.014 for deuterium)
L_char      = 2  #Characteristic length

# Constants
mu0      = 4*np.pi*1e-7
e_charge = 1.60217662e-19
e_mass   = 9.1093837e-31
ion_mass = 1.67e-27 * atomic_mass
kB       = 1.380649 * 1e-23
varep_0  = 8.854e-12

# Functions
def v_thermal_e(Te):
    Te_SI    = Te * e_charge / kB
    return (kB*Te_SI/e_mass)**0.5

def v_perp_e(Te):
    return v_thermal_e(Te) * np.sqrt(2)

def v_thermal_i(Ti):
    Ti_SI    = Ti * e_charge / kB
    return (kB*Ti_SI/ion_mass)**0.5

def v_perp_i(Ti):
    return v_thermal_i(Ti) * np.sqrt(2)

def v_alfven(Bt, ni, ion_mass):
    return Bt/np.sqrt(ion_mass*ni*mu0)

def ion_gyro_freq(Bt):
    return Zion * e_charge * Bt / ion_mass

def e_gyro_freq(Bt):
    return e_charge * Bt / e_mass

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

def eta_spitzer(Te, ne, Zeff):
    coef_Zeff = Zeff*(1.+1.198*Zeff+0.222*Zeff**2)/(1.+2.966*Zeff+0.753*Zeff**2) / ((1.+1.198+0.222)/(1.+2.966+0.753))
    eta = 1.65e-9 * coulomb_log(Te, ne) * (Te*0.001)**(-1.5) * coef_Zeff # From Wesson
    return eta

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
    return 6.60e17 * (Te*0.001)**1.5 / (ni * Zeff**4 * coulomb_log(Te, ne)) * (ion_mass/e_mass)**0.5

def e_collision_mean_free_path(Te, ne, Zeff, ni):
    return v_thermal_e(Te) * collision_time_ei(Te, ne, Zeff, ni)
    
# Quantities
class quantity:
    def __init__(self, name, units, value):
        self.name  = name
        self.units = units
        self.value = value
quantities = []

quantities.append(  quantity('v_thermal_e        ','m/s   ', v_thermal_e(Te) )     )
quantities.append(  quantity('v_thermal_i        ','m/s   ', v_thermal_i(Ti) )     )
quantities.append(  quantity('v_alfven           ','m/s   ', v_alfven(Bt, ni, ion_mass) )   )
quantities.append(  quantity('alfven_time        ','s     ', L_char/v_alfven(Bt, ni, ion_mass) )   )
quantities.append(  quantity('ion_gyro_frequency ','Hz    ', ion_gyro_freq(Bt)    )   )
quantities.append(  quantity('e_gyro_frequency   ','Hz    ', e_gyro_freq(Bt)    )   )
quantities.append(  quantity('Ion-Ion col time   ','s     ', collision_time_ii(Ti, ne, Zeff, ni)    )   )
quantities.append(  quantity('Ion-ele col time   ','s     ', collision_time_ei(Ti, ne, Zeff, ni)    )   )
quantities.append(  quantity('Ele-ele col time   ','s     ', collision_time_ee(Ti, ne, Zeff, ni)    )   )
quantities.append(  quantity('ion_gyro_radius    ','m     ', ion_gyro_radius(Ti, Bt)    )   )
quantities.append(  quantity('e_gyro_radius      ','m     ', e_gyro_radius(Te, Bt)    )   )
quantities.append(  quantity('Ele mean free path ','m     ',  e_collision_mean_free_path(Te, ne, Zeff, ni)   )   )
quantities.append(  quantity('Coulomb_log        ','      ', coulomb_log(Te, ne)    )   )
quantities.append(  quantity('eta_Spitzer        ','Ohm m ', eta_spitzer(Te, ne, Zeff)    )   )

nquant = len(quantities)
for quant in quantities:
    print('%s = %.2e %s' % (quant.name, quant.value, quant.units))