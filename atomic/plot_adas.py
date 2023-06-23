import numpy as np
import matplotlib.pyplot as plt
from cherab.core.atomic import neon, hydrogen, argon, tungsten
from cherab.openadas import OpenADAS
from scipy.optimize import lsq_linear
from cherab.openadas.repository import populate

def get_rates_recombination(element):
    """
    load recombinatio rates for all ionic stages
    """
    coef_recom = {}
    for i in np.arange(1, elem.atomic_number + 1):
        coef_recom[i] = adas.recombination_rate(element, int(i))

    return coef_recom

def get_rates_ionisation(element):
    """
    load ionisation rates for all ionic stages
    :param element:
    :return:
    """
    coef_ionis = {}
    for i in np.arange(0, elem.atomic_number):
        coef_ionis[i] = adas.ionisation_rate(element, int(i))

    return coef_ionis


def get_LT_radiation(element):
    """
    load radiation recombinatio rates for all ionic stages
    """

    coef_rad = {}
    for i in np.arange(0, elem.atomic_number):
        coef_rad[i] = adas.line_radiated_power_rate(element, int(i))

    return coef_rad


def get_BR_radiation(element):
    """
    load radiation recombinatio rates for all ionic stages
    """

    coef_rad = {}
    for i in np.arange(1, elem.atomic_number+1):
        coef_rad[i] = adas.continuum_radiated_power_rate(element, int(i))

    return coef_rad


# Z distribution a la Di
def solve_ion_distribution(element, n_e, t_e, coef_ion, coef_recom):


    atomic_number = element.atomic_number
    prob = np.zeros(atomic_number+1)
    prob[0] = 1.0

    for i in range(1, atomic_number+1):
        ion_rate = coef_ion[i-1](n_e, t_e)
        rec_rate = coef_recom[i](n_e, t_e)

        prob[i]  = prob[i-1] * ion_rate / rec_rate 

    
    prob = prob / sum(prob)

    return prob


def coronal_output(elem, rates_ion, rates_recom, rates_line, rates_contin, nD, nimp, Te_array):

    n_states = elem.atomic_number + 1 

    # Iterate a few times to find consistent n_e and charge distribution
    ntries  = 5
    n_Te    = len(Te_array)
    ne0     = np.full((n_Te), nD+nimp)
    Z_av    = np.zeros(n_Te)
    Z_eff   = np.zeros(n_Te)
    Lrad    = np.zeros(n_Te)
    prob_Z  = np.zeros((elem.atomic_number + 1, len(Te_array)))
    
    for it in range(ntries):
      
        for i, te in enumerate(Te_array):
            prob_Z[:, i]  = solve_ion_distribution(elem, min(ne0[i],2e21),  te, rates_ion, rates_recom)
        
            # calculate average charge state as function of Te
            Z_av[i] = 0.0
            for j in range(0, n_states):
                Z_av[i] = Z_av[i] + float(j)*prob_Z[j,i]
        
            # Correct electron density
            ne0[i] = nD + nimp * Z_av[i]

            # Zeff
            numerator   = nD * 1.0**2 
            denominator = nD * 1.0 
            for j in range(0, n_states):
                numerator   = numerator   + nimp * float(j)**2 *prob_Z[j,i]
                denominator = denominator + nimp * float(j)    *prob_Z[j,i]

            Z_eff[i] = numerator / denominator 

            # Radiated power
            Lrad[i] = 0.0
            for j in range(0, n_states):
                if j==0:
                    Lrad[i] = Lrad[i] + (rates_line[j](ne0[i],te) ) * prob_Z[j,i]
                elif j==n_states-1:
                    Lrad[i] = Lrad[i] + (rates_contin[j](ne0[i],te)) * prob_Z[j,i]
                else:
                    Lrad[i] = Lrad[i] + (rates_line[j](ne0[i],te) + rates_contin[j](ne0[i],te)) * prob_Z[j,i]
                
    return prob_Z, Z_av, ne0, Z_eff, Lrad



# Te in eV and ne in SI
def eta_spitzer(Te_array, Z_eff, n_e):

    n_Te  = len(Te_array)
    eta   = np.zeros(n_Te)

    for i, Te in enumerate(Te_array):
        if Te < 10:
              Coulomb_Log = 23.0000 - np.log((n_e[i]*1e-6)**0.5*Te**(-1.5))  
        else:
              Coulomb_Log = 24.1513 - np.log((n_e[i]*1e-6)**0.5*Te**(-1.0))  

        coef_Zeff = Z_eff[i]*(1.+1.198*Z_eff[i]+0.222*Z_eff[i]**2)/(1.+2.966*Z_eff[i]+0.753*Z_eff[i]**2) / ((1.+1.198+0.222)/(1.+2.966+0.753))

        eta[i] = 1.65e-9 * Coulomb_Log * (Te*0.001)**(-1.5) * coef_Zeff

    return eta

# Gets equilibrium Te and others
def find_balance(elem, rates_ion, rates_recom, rates_line, rates_contin, nD, nimp, Te_array):

    coronal_data = coronal_output(elem, rates_ion, rates_recom, rates_line, rates_contin, nD, nimp, Te_array)
    
    prob_Z = coronal_data[0]
    Z_av   = coronal_data[1]
    n_e    = coronal_data[2]
    Z_eff  = coronal_data[3]
    Lrad   = coronal_data[4]
    eta    = eta_spitzer(Te_array, Z_eff, n_e) 
    
    ohmic         = eta * J_av**2
    radiation     = n_e * nimp * Lrad

    # Look minimum balance
    diff      = abs( ohmic - radiation )

    min_index = np.argmin(diff)
    Te_eq     = Te_array[min_index]
    eta_eq    = eta[min_index]
    ne_eq     = n_e[min_index]
    Zeff_eq   = Z_eff[min_index]

    return Te_eq, eta_eq, ne_eq, Zeff_eq

# initialise the atomic data provider
adas = OpenADAS(permit_extrapolation=True)

#elem = tungsten 
elem = neon 
temperature_steps = 100

Ip       = 15e6
Area_pol = 21
J_av     = Ip/Area_pol
J_av = 1e6

nD   = 5e19
nD2  = 1e20 
nimp = 2e18

n_states = elem.atomic_number + 1 

# Collect rate coefficients
rates_ion   = get_rates_ionisation(elem)
rates_recom = get_rates_recombination(elem)
rates_line  = get_LT_radiation(elem)
rates_contin= get_BR_radiation(elem)

Te_min =2 #rates_recom[1].raw_data["te"].min()
Te_max =1e4 # rates_recom[1].raw_data["te"].max()

# Initialize temperatures
Te_array = [10 ** x for x in np.linspace(np.log10(Te_min),
                                         np.log10(Te_max),
                                         num=temperature_steps)]

coronal_data = coronal_output(elem, rates_ion, rates_recom, rates_line, rates_contin, nD, nimp, Te_array)

prob_Z = coronal_data[0]
Z_av   = coronal_data[1]
n_e    = coronal_data[2]
Z_eff  = coronal_data[3]
Lrad   = coronal_data[4]
eta    = eta_spitzer(Te_array, Z_eff, n_e) 

ohmic         = eta * J_av**2
radiation     = n_e * nimp * Lrad

coronal_data = coronal_output(elem, rates_ion, rates_recom, rates_line, rates_contin, nD2, nimp, Te_array)

prob_Z2= coronal_data[0]
Z_av2  = coronal_data[1]
n_e2   = coronal_data[2]
Z_eff2 = coronal_data[3]
Lrad2  = coronal_data[4]
eta2   = eta_spitzer(Te_array, Z_eff2, n_e2) 

label1 = "n_D = "+ str(nD) + ", n_imp = "+ str(nimp) + " m^-3"
label2 = "n_D = "+ str(nD2)+ ", n_imp = "+ str(nimp) + " m^-3"

# Deuterium scan
nD_min = 5e19    #rates_recom[1].raw_data["te"].min()
nD_max = 1e22   #rates_recom[1].raw_data["te"].max()
nD_n   = 10
nD_array = [10 ** x for x in np.linspace(np.log10(nD_min), np.log10(nD_max), num=nD_n)]

eq_Te   = np.zeros(nD_n)
eq_eta  = np.zeros(nD_n)
eq_Zeff = np.zeros(nD_n)

for i, nD in enumerate(nD_array):

    balance = find_balance(elem, rates_ion, rates_recom, rates_line, rates_contin, nD, nimp, Te_array)

    eq_Te[i]   = balance[0]
    eq_eta[i]  = balance[1]
    eq_Zeff[i] = balance[3]


# Plot charge state distribution
plt.figure()
for i in range(elem.atomic_number + 1):
    pl = plt.semilogx(Te_array, prob_Z[i, :], label='{0} {1}+'.format(elem.symbol, i))
    plt.semilogx(Te_array, prob_Z2[i, :], '--',
                 color=pl[0].get_color(), lw=2)
plt.plot([], [], "k-",  label=label1 )
plt.plot([], [], "k--", label=label2 )
plt.xlabel("Electron Temperature (eV)")
plt.title('Fractional Abundance')
plt.legend()
plt.savefig(elem.name+'_charge_distribution.png')


# Plot average charge 
plt.figure()
plt.semilogx(Te_array,Z_av , label=label1)
plt.semilogx(Te_array,Z_av2, label=label2)
plt.grid(True)
plt.xlabel("Electron Temperature (eV)")
plt.title('average charge state of Neon')
plt.legend()
plt.savefig(elem.name+'_average_charge.png')


# Plot Zeff 
plt.figure()
plt.semilogx(Te_array,Z_eff,  label=label1)
plt.semilogx(Te_array,Z_eff2, label=label2)
plt.grid(True)
plt.xlabel("Electron Temperature (eV)")
plt.title('Z_eff')
plt.legend()
plt.savefig(elem.name+'_Zeff.png')


# Plot ne 
plt.figure()
plt.semilogx(Te_array,n_e,  label=label1)
plt.semilogx(Te_array,n_e2, label=label2)
plt.grid(True)
plt.xlabel("Electron Temperature (eV)")
plt.title('Electron density')
plt.legend()
plt.savefig(elem.name+'_ne.png')


# Plot Lrad 
plt.figure()
plt.plot(Te_array,Lrad,  label=label1)
plt.plot(Te_array,Lrad2, label=label2)
ax = plt.gca()
ax.set_xscale('log')
ax.set_yscale('log')
plt.grid(True)
plt.xlabel("Electron Temperature (eV)")
plt.title('Radiation power')
plt.legend()
plt.savefig(elem.name+'_Lrad.png')


# Plot eta_spitzer
plt.figure()
plt.plot(Te_array,eta ,  label=label1)
plt.plot(Te_array,eta2 , label=label2)
plt.xscale('log')
plt.yscale('log')
plt.grid(True)
plt.xlabel("Electron Temperature (eV)")
plt.title('Resistivity (Ohm.m)')
plt.legend()
plt.savefig(elem.name+'_eta_Spitzer.png')

# Neon balance 
plt.figure()
plt.plot(Te_array,ohmic ,  label='ohmic')
plt.plot(Te_array,radiation , label='radiation')
plt.xscale('log')
plt.yscale('log')
plt.grid(True)
plt.xlabel("Electron Temperature (eV)")
plt.title('Power competetion')
plt.legend()
plt.savefig(elem.name+'_balance.png')

# D balance 
plt.figure()
plt.plot(nD_array, eq_Te,  label='Te')
plt.plot(nD_array, eq_eta/eq_eta[0]*4,  label='eta/eta0 x4')
plt.plot(nD_array, eq_Zeff*4,  label='Zeff x4')
plt.xscale('log')
plt.grid(True)
plt.xlabel("n_D (m^-3)")
plt.title('Equilibrium parameters')
plt.legend()
plt.savefig(elem.name+'_nD_scan.png')


