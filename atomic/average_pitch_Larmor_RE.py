import numpy as np
import matplotlib.pyplot as plt
from cherab.core.atomic import neon, hydrogen, argon, tungsten
from cherab.openadas import OpenADAS
from scipy.optimize import lsq_linear
from cherab.openadas.repository import populate
from scipy import constants
from scipy.special import iv
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

nD   = 1e20
nimp = 2e19

n_states = elem.atomic_number + 1 

# Collect rate coefficients
rates_ion   = get_rates_ionisation(elem)
rates_recom = get_rates_recombination(elem)
rates_line  = get_LT_radiation(elem)
rates_contin= get_BR_radiation(elem)

Te_min =2 #rates_recom[1].raw_data["te"].min()
Te_max =1e1 # rates_recom[1].raw_data["te"].max()

# Initialize temperatures
Te_array = [10 ** x for x in np.linspace(np.log10(Te_min),
                                         np.log10(Te_max),
                                         num=temperature_steps)]




e_ch = constants.e
m_e  = constants.m_e
c    = constants.c
Bphi = 5.3

# Couloumb log, as in Hesslow 18
def coul_log(Te_eV, ne_SI):
    return 14.6 + 0.5 * np.log(Te_eV / (ne_SI*1e-20)) 

# Deflection frequency ( eq 5 in Hesslow 18), background imp is fully ionized
def nu_D(Te_eV, ne_SI, Z_eff, p_norm, nimp, prob_Z):

    n_Z   = len(prob_Z) - 1

    ln_aj = np.array( [4.7, 4.6, 4.5, 4.4, 4.3, 4.1, 4.0, 3.7, 3.2, 3.1] )

    ln_lamb_c  = coul_log(Te_eV, ne_SI)
    gamma      = np.sqrt( 1 + p_norm**2 )
    ln_lamb_ee = ln_lamb_c + 0.5 * np.log( gamma - 1 )
    ln_lamb_ei = ln_lamb_c +       np.log( p_norm * np.sqrt(2.0) )

    imp_screen = 0.0

    for j in range(n_Z):
        imp_screen += prob_Z[j] * ( (n_Z**2 - j**2) * (ln_aj[j] + np.log(p_norm)) - 2.0/3.0 * (n_Z-j)**2)

    imp_screen = imp_screen * nimp / ne_SI 

#    print( "nz =%s, lambc=%s, lambe=%s, lambi=%s, gamma=%s"%(n_Z, ln_lamb_c, ln_lamb_ee, ln_lamb_ei, gamma))

    nuD = ( ln_lamb_ee + Z_eff * ln_lamb_ei + imp_screen ) / ln_lamb_c    
    
    return nuD


def Ztot_star(Te_eV, ne_SI, Z_eff, p_norm, nimp, prob_Z):
    return nu_D(Te_eV, ne_SI, Z_eff, p_norm, nimp, prob_Z) - 1 

def E_crit_eff(Ztot_star, Bphi, Te_eV, ne_SI):
    tau_rad = 1.5*ne_SI * coul_log(Te_eV, ne_SI)*m_e/ (8.8e-12*Bphi**2)
    nuD     = Ztot_star +1 
    return 1 + ( nuD / np.sqrt(tau_rad) ) / (1/8+nuD**2/tau_rad)**(1/6)

# Ke is electron energy in eV, Ec is the electric field normalized to the critical electric field
def Af(Z_star, Ke, Ec):
    return (Ec + 1)/(Z_star+1) * (Ke*e_ch / (m_e*c**2) +1)

def average_pitch(Z_star, Ke, Ec):
    A = Af(Z_star, Ke, Ec)
    return 1/A - 1/np.tanh(A)

# Ke in eV, output in mm
def larmor_radius(Z_star, Ke, Ec):
    A = Af(Z_star, Ke, Ec)
    pnorm  = np.sqrt( (Ke*e_ch/(m_e*c**2))**2 - 1)
    p      = pnorm * m_e * c
    p_perp =  p * np.pi * iv(1, A) / (np.exp(A) - np.exp(-A))
    return p_perp / (e_ch*Bphi) * 1e3


nimp_array = [10 ** x for x in np.linspace(np.log10(1e18), np.log10(1e20), num=20)]
nD_array   = [10 ** x for x in np.linspace(np.log10(1e20), np.log10(1e21), num=2)]
T_array    = np.linspace(1, 5, num=2)
Ekin       = 0.5e7
pnorm      = np.sqrt( (Ekin*e_ch/(m_e*c**2))**2 - 1)
fig, axes = plt.subplots(4, 1, sharex=True, figsize=(5,10))
plt.subplots_adjust(hspace=.0)

# Coronal equilibrium data for different temperatures
for nD in nD_array:
    for Te in T_array:
        Ece = []
        pitch = []
        larmor = []
        Zstar = []
        for nimp in nimp_array:
  
            coronal_data = coronal_output(elem, rates_ion, rates_recom, rates_line, rates_contin, nD, nimp, np.array([Te,Te]))
            prob_Z = coronal_data[0]
            n_e    = coronal_data[2]
            Z_eff  = coronal_data[3]
  
            Zs  = Ztot_star(Te, n_e[0], Z_eff[0], pnorm, nimp, prob_Z[:,0])
            Ect = E_crit_eff(Zs, Bphi, Te, n_e[0])
            pt  = average_pitch(Zs, Ekin, Ect)
            lm  = larmor_radius(Zs, Ekin, Ect)

            Zstar.append(Zs)
            Ece.append(Ect)
            pitch.append(pt)
            larmor.append(lm)
  
        axes[0].plot(nimp_array, Zstar, label=r'$n_D=$%s, $T_e=$%s'%(nD,Te))  # Add label for legend
        axes[1].plot(nimp_array, Ece, label=r'$n_D=$%s, $T_e=$%s'%(nD,Te))  # Add label for legend
        axes[2].plot(nimp_array, pitch, label=r'$n_D=$%s, $T_e=$%s'%(nD,Te))  # Add label for legend
        axes[3].plot(nimp_array, larmor, label=r'$n_D=$%s, $T_e=$%s'%(nD,Te))  # Add label for legend

labels = [r'$Z_{*}$',r'$E / E_c$ ',  r'$<\xi>=<v_\parallel/v>$',  r'$<p_\perp>/(eB)$ [mm]']

for ax, lab in zip(axes, labels):
    ax.set_ylabel(lab,fontsize=14)
    ax.grid(True)  # Add grid lines
    ax.set_xscale('log')

axes[-1].set_xlabel(r'$n_{Neon}$ [m$^{-3}]$ ',fontsize=16)
axes[-1].legend(loc='center left', bbox_to_anchor=(1, 0.5))  # Place legend outside the plot
#plt.tight_layout()
#plt.show()        
fig.savefig("average_pitch_5MeV", bbox_inches='tight')
