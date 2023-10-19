# This python script calculates the RE electron pitch angle distribution
# and effective critical electric field, as given by the formulas 
# of the following article: 
#      L Hesslow et al 2018 Plasma Phys. Control. Fusion 60 074010
# Coronal equilibrium is assumed when using the impurity charge distribution

import numpy as np
import matplotlib.pyplot as plt
from cherab.core.atomic import neon, hydrogen, argon, tungsten
from cherab.openadas import OpenADAS
from cherab.openadas.repository import populate
from scipy import constants
from scipy.special import iv
import matplotlib
matplotlib.use('TKagg')

# Routines to get the recombination, ionisation & radiation rates from open ADAS
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

# Routine to calculate the fractional abundance of the different charge states
# of an impurity assuming coronal equilibrium
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

# Calculates different quantities (Zeff, ne, Lrad...) from the coronal equilibrium distribution
# for a set of temperatures
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


####################################################
# Formulas from  Hesslow et al 2018 PPCF 60 074010 #
####################################################
# Coulomb logarithm, Eq. 2
def coul_log(Te_eV, ne_SI):
    return 14.6 + 0.5 * np.log(Te_eV / (ne_SI*1e-20)) 

# Slowing-down frequency coefficients, Eqs. 11-12 
def nu_s_coef(Te_eV, ne_SI, Z_eff, nimp, prob_Z):

    n_Z       = len(prob_Z) - 1          # Charge number
    ln_lamb_c = coul_log(Te_eV, ne_SI)

    # Neon coefficients in Table A1
    ln_ovI    = np.array( [8.2, 8.0, 7.9, 7.7, 7.5, 7.3, 7.0, 6.6, 5.9, 5.8] )

    sum1 = 0.0; sum2 = 0.0
    for j in range(n_Z):
        Nej   = float(n_Z-j)
        nj    = nimp * prob_Z[j] 
        sum1 += nj / ne_SI * Nej * ( ln_ovI[j] - 1 )
        sum2 += nj / ne_SI * Nej * 3

    nu_s0 = 1 + sum1 / ln_lamb_c
    nu_s1 = 0.5 * (sum2 + 1) / ln_lamb_c
    
    return nu_s0, nu_s1

# Deflection frequency coefficients, Eqs. 7,8 
def nu_D_coeff(Te_eV, ne_SI, Z_eff, nimp, prob_Z, nD):

    n_Z       = len(prob_Z) - 1         # Charge number
    ln_lamb_c = coul_log(Te_eV, ne_SI)

    # Neon coefficients in Table A1
    ln_aj = np.array( [4.7, 4.6, 4.5, 4.4, 4.3, 4.1, 4.0, 3.7, 3.2, 3.1] )

    sum1 = 0.0
    for j in range(n_Z):
        nj    = nimp * prob_Z[j] 
        sum1 += nj / ne_SI * ( (n_Z**2 - j**2) * ln_aj[j]  - 2.0/3.0 * (n_Z-j)**2 ) 

    sum2  = ( nimp * n_Z**2 + nD  ) / ne_SI

    nu_D0 = 1.0 + Z_eff + sum1 / ln_lamb_c 
    nu_D1 = sum2 / ln_lamb_c
    
    return nu_D0, nu_D1

# Generalized normalized deflection frequency, Eq. 5
def nu_D(Te_eV, ne_SI, Z_eff, p_norm, nimp, prob_Z):

    n_Z       = len(prob_Z) - 1         # Charge number
    ln_lamb_c = coul_log(Te_eV, ne_SI)

    # Neon coefficients in Table A1
    ln_aj = np.array( [4.7, 4.6, 4.5, 4.4, 4.3, 4.1, 4.0, 3.7, 3.2, 3.1] )

    gamma      = np.sqrt( 1 + p_norm**2 )
    ln_lamb_ee = ln_lamb_c + 0.5 * np.log( gamma - 1 )
    ln_lamb_ei = ln_lamb_c +       np.log( p_norm * np.sqrt(2.0) )

    sum1 = 0.0
    for j in range(n_Z):
        nj    = nimp * prob_Z[j] 
        sum1 += nj / ne_SI * ( (n_Z**2 - j**2) * (ln_aj[j] + np.log(p_norm)) - 2.0/3.0 * (n_Z-j)**2)

    nuD = ( ln_lamb_ee + Z_eff * ln_lamb_ei + sum1 ) / ln_lamb_c    
    
    return nuD

# Normalized effective electric field, Eq. 23
def E_crit_eff(Te_eV, ne_SI, Z_eff, nimp, prob_Z, Bphi, nD):

    nu_s  = nu_s_coef(Te_eV, ne_SI, Z_eff, nimp, prob_Z)
    nu_s0 = nu_s[0]
    nu_s1 = nu_s[1]

    nu_D  = nu_D_coeff(Te_eV, ne_SI, Z_eff, nimp, prob_Z, nD)
    nu_D0 = nu_D[0]
    nu_D1 = nu_D[1]

    alpha = 1/137  # fine constant
    n_Z   = len(prob_Z) - 1
    
    phi_br0 = (alpha/coul_log(Te_eV, ne_SI)) * ( nimp * n_Z**2 + nD  ) / ne_SI * 0.35 # Eq. 18
    phi_br1 = (alpha/coul_log(Te_eV, ne_SI)) * ( nimp * n_Z**2 + nD  ) / ne_SI * 0.20 # Eq. 18

    tau_syn_inv = 1 / (15.44 * coul_log(Te_eV, ne_SI) ) * Bphi**2 / (ne_SI*1e-20)  # Eq 15.

    # Interate to converge to a solution
    E_ratio = 2.0
    for i in range(5):
        delta   = (nu_D0/nu_s1**2) * ( nu_D0 * tau_syn_inv / E_ratio + phi_br0 + phi_br1*np.log(nu_D0/(2*nu_s1)) ) # Eq. 24
        E_ratio = nu_s0 + nu_s1*( (1+nu_D1/nu_D0)*np.log(nu_D0/(2*nu_s1)) + np.sqrt(2*delta+1) ) 

    return E_ratio

# Ke is electron energy in eV, Ec is the electric field normalized to the critical electric field
def Ap(nuD, p_norm, Ec_ratio):
    return 2*Ec_ratio / nuD * p_norm**2 / np.sqrt( 1 + p_norm**2 )

# Average pitch angle
def average_pitch(nuD, p_norm, Ec_ratio):
    A = Ap(nuD, p_norm, Ec_ratio)
    return 1/A - 1/np.tanh(A)

# Output in mm
def average_larmor_radius(nuD, p_norm, Ec_ratio):
    A      = Ap(nuD, p_norm, Ec_ratio)
    p      = pnorm * m_e * c
    p_perp =  p * np.pi * iv(1, A) / (np.exp(A) - np.exp(-A))
    return p_perp / (e_ch*Bphi) * 1e3

# Main program
if __name__ == "__main__":

    # initialise the atomic data provider
    adas = OpenADAS(permit_extrapolation=True)

    #elem = tungsten 
    elem = neon 
    temperature_steps = 100

    n_states = elem.atomic_number + 1 

    # Collect rate coefficients
    rates_ion    = get_rates_ionisation(elem)
    rates_recom  = get_rates_recombination(elem)
    rates_line   = get_LT_radiation(elem)
    rates_contin = get_BR_radiation(elem)

    # Initialize temperatures in eV
    Te_min   = 1;    Te_max = 1e1 
    Te_array = [10 ** x for x in np.linspace(np.log10(Te_min),
                                            np.log10(Te_max),
                                            num=temperature_steps)]
    e_ch = constants.e
    m_e  = constants.m_e
    c    = constants.c
    Bphi = 5.3

    nimp_array = [10 ** x for x in np.linspace(np.log10(1e18), np.log10(1e20), num=60)]
    nD_array   = [10 ** x for x in np.linspace(np.log10(1e20), np.log10(1e21), num=2)]
    T_array    = np.linspace(1, 5, num=2)
    
    Ekin       = 1.0e7
    pnorm      = np.sqrt( (Ekin*e_ch/(m_e*c**2))**2 - 1)
    print("Energy = %s eV and p_norm = %s"%(Ekin, pnorm))
    
    fig, axes  = plt.subplots(4, 1, sharex=True, figsize=(5,10))
    plt.subplots_adjust(hspace=.0)

    # Coronal equilibrium data for different temperatures
    for nD in nD_array:
        for Te in T_array:
            Ece    = [];  pitch = [];  larmor = [];  nuD_norm = []
            for nimp in nimp_array:
    
                coronal_data = coronal_output(elem, rates_ion, rates_recom, rates_line, rates_contin, nD, nimp, np.array([Te,Te]))
                prob_Z = coronal_data[0][:,0]
                n_e    = coronal_data[2][0]
                Z_eff  = coronal_data[3][0]

                nuD = nu_D(Te, n_e, Z_eff, pnorm, nimp, prob_Z)
                Ect = E_crit_eff(Te, n_e, Z_eff, nimp, prob_Z, Bphi, nD)

                pt  = average_pitch(nuD, pnorm, Ect)
                lm  = average_larmor_radius(nuD, pnorm, Ect)

                nuD_norm.append(nuD)
                Ece.append(Ect)
                pitch.append(pt)
                larmor.append(lm)
    
            axes[0].plot(nimp_array, nuD_norm, label=r'$n_D=$%s, $T_e=$%s'%(nD,Te))  # Add label for legend
            axes[1].plot(nimp_array, Ece, label=r'$n_D=$%s, $T_e=$%s'%(nD,Te))  # Add label for legend
            axes[2].plot(nimp_array, np.abs(pitch), label=r'$n_D=$%s, $T_e=$%s'%(nD,Te))  # Add label for legend
            axes[3].plot(nimp_array, larmor, label=r'$n_D=$%s, $T_e=$%s'%(nD,Te))  # Add label for legend

    labels = [r'$\overline{\nu}_D$',r'$E / E_c$ ',  r'$<\xi>=<v_\parallel/v>$',  r'$<p_\perp>/(eB)$ [mm]']

    for ax, lab in zip(axes, labels):
        ax.set_ylabel(lab,fontsize=14)
        ax.grid(True)  # Add grid lines
        ax.set_xscale('log')

    axes[3].set_xlabel(r'$n_{Neon}$ [m$^{-3}]$ ',fontsize=16)
    axes[3].legend(loc='center left', bbox_to_anchor=(1, 0.5))  # Place legend outside the plot
    plt.tight_layout()
    plt.show()        
    #fig.savefig("average_pitch_Neon", bbox_inches='tight')
