import numpy as np
import matplotlib.pyplot as plt
from cherab.core.atomic import neon, hydrogen, argon, tungsten
from cherab.openadas import OpenADAS
from scipy.optimize import fsolve
from scipy.interpolate import interp1d
from scipy.optimize import lsq_linear
from cherab.openadas.repository import populate
import matplotlib
matplotlib.use('TKagg')


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
def find_balance(elem, rates_ion, rates_recom, rates_line, rates_contin, nD, nimp, Te_arr):

    coronal_data = coronal_output(elem, rates_ion, rates_recom, rates_line, rates_contin, nD, nimp, Te_arr)
    
    prob_Z = coronal_data[0]
    Z_av   = coronal_data[1]
    n_e    = coronal_data[2]
    Z_eff  = coronal_data[3]
    Lrad   = coronal_data[4]
    eta    = eta_spitzer(Te_arr, Z_eff, n_e) 
    
    ohmic         = eta * J_av**2
    radiation     = n_e * nimp * Lrad
#    gamma_sheath  = 8.0;   Rgeo = 6.2;   q_halo=2.0
    qpar   = 2*np.sqrt(2*gamma) * gamma_sheath * n_e * Te_arr**(1.5) / np.sqrt( m_i ) * e_ch**(1.5)
    par_loss      = qpar / Lpar 

    # Look minimum balance
    diff      = (ohmic - radiation - par_loss)**2

    diff_interp = interp1d(Te_arr, diff, kind='cubic', bounds_error=False)  # 'linear', 'quadratic', 'cubic'

    Te_arr2 = np.linspace(Te_arr[0], Te_arr[-1], num=1000)
    diff2   = diff_interp( Te_arr2 )
#    print( Te_eq)

    min_index = np.argmin(diff2)
    Te_eq     = Te_arr2[min_index]

    return Te_eq


# initialise the atomic data provider
adas = OpenADAS(permit_extrapolation=True)

# elem = tungsten 
elem = neon 

n_states = elem.atomic_number + 1 

# Collect rate coefficients
rates_ion   = get_rates_ionisation(elem)
rates_recom = get_rates_recombination(elem)
rates_line  = get_LT_radiation(elem)
rates_contin= get_BR_radiation(elem)

# Initialize temperatures
temperature_steps = 25
Te_min = 1  #rates_recom[1].raw_data["te"].min()
Te_max = 25 # rates_recom[1].raw_data["te"].max()
Te_array = np.linspace(Te_min, Te_max, num=temperature_steps)

# Neon scan
n_imp_min = 1e18    #rates_recom[1].raw_data["te"].min()
n_imp_max = 2e20   #rates_recom[1].raw_data["te"].max()
n_imp_n   = 6 
n_imp_array = [10 ** x for x in np.linspace(np.log10(n_imp_min), np.log10(n_imp_max), num=n_imp_n)]
n_imp_array[0] = 0

nD_min = 1e19    #rates_recom[1].raw_data["te"].min()
nD_max = 3e21   #rates_recom[1].raw_data["te"].max()
nD_n   = 300
nD_array = [10 ** x for x in np.linspace(np.log10(nD_min), np.log10(nD_max), num=nD_n)]

n_imp_array = np.array([0, 1e17, 1e18, 1e19, 1e20])

Ip       = 15e6
Area_pol = 21
#J_av     = Ip/Area_pol

e_ch   = 1.602e-19
mu0    = 4*np.pi*1e-7
m_i    = 1.67e-27 * 2
B_phi  = 5.3
R_geo  = 6.2
gamma_sheath = 8
gamma  = 5.0/3.0
#q_halo = 2.0 
ne_w   = 1e21
coupled = True

# D balance 
fig, axes = plt.subplots(2, 1, sharex=True, figsize=(5,7.5))
plt.subplots_adjust(hspace=.0)
r_colors = plt.cm.rainbow(np.linspace(0, 1, len(n_imp_array)))
q_array = np.array([1,2,3])
ltype   = ["-", "--", ":"]

for iq, q0 in enumerate(q_array):

    Lpar   = 2*np.pi * R_geo * q0
    J_av   = 2*B_phi / (q0 * R_geo * mu0)

    print(J_av)

    for i, n_imp in enumerate(n_imp_array):
    
        eff   = []
        Te_eq = []
    
        for j, nD in enumerate(nD_array):
    
            Te_bal = find_balance(elem, rates_ion, rates_recom, rates_line, rates_contin, nD, n_imp, Te_array)
        
            Te_bal = find_balance(elem, rates_ion, rates_recom, rates_line, rates_contin, nD, n_imp, np.linspace(0.5*Te_bal,2*Te_bal,num=20))

            coronal_data = coronal_output(elem, rates_ion, rates_recom, rates_line, rates_contin, nD, n_imp, np.array([Te_bal,Te_bal]))
        
            n_e    = coronal_data[2][0]
            Lrad   = coronal_data[4][0]
            Jsat   = e_ch * n_e * np.sqrt( 2 * gamma * e_ch * Te_bal / m_i )
        
            if coupled: 
                ne_heat = n_e 
            else:
                ne_heat = ne_w
        
            qpar   = 2*np.sqrt(2*gamma) * gamma_sheath * ne_heat * Te_bal**(1.5) / np.sqrt( m_i ) * e_ch**(1.5)
            rad    = n_e * n_imp * Lrad
            qloss  = qpar / Lpar
    
            if (J_av <= Jsat):
               eff.append( rad / ( rad + qloss ))
               Te_eq.append( Te_bal )
            else:
               eff.append(  np.nan )
               Te_eq.append( np.nan )

        if (iq==0): 
            axes[0].plot(nD_array, Te_eq, ltype[iq], color=r_colors[i], label=r'$n_{imp}$=%s m$^{-3}$'%("{:.0e}".format(n_imp) ))
            axes[1].plot(nD_array, eff,   ltype[iq], color=r_colors[i], label=r'$n_{imp}$=%s m$^{-3}$'%("{:.0e}".format(n_imp) ))
        else:
            axes[0].plot(nD_array, Te_eq, ltype[iq], color=r_colors[i], label='')
            axes[1].plot(nD_array, eff,   ltype[iq], color=r_colors[i], label='')

axes[-1].set_xlim(1e19, 3e21)
axes[-1].set_xscale('log')
#axes[1].set_yscale('log')
axes[0].grid(True, which="both", ls="-", linewidth=0.5)
axes[-1].grid(True, which="both", ls="-", linewidth=0.5)
axes[-1].set_xlabel(r"$n_D$ [m$^{-3}$ ]")
axes[0].set_ylabel(r"$T_e^{eq}$ [eV]")
axes[1].set_ylabel(r"radiated / Ohmic power")
axes[1].legend(loc='upper left', bbox_to_anchor=(0., 0.7), fontsize=9  )
#plt.show()
fig.suptitle(r"$q_0=1$: Solid, $q_0=2$: Dashed, $q_0=3$: Dotted", y=0.92)
fig.savefig("radiation_efficiency",dpi=300, bbox_inches='tight')

