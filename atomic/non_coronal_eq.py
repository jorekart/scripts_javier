import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
from cherab.core.atomic import neon, hydrogen, argon, tungsten
from cherab.openadas import OpenADAS
from scipy.optimize import lsq_linear
from cherab.openadas.repository import populate

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

n_states = elem.atomic_number + 1 

# Collect rate coefficients
rates_ion   = get_rates_ionisation(elem)
rates_recom = get_rates_recombination(elem)
rates_line  = get_LT_radiation(elem)
rates_contin= get_BR_radiation(elem)

def add_random_noise(array):
    noise = np.random.uniform(low=-0.001, high=0.001, size=array.shape)
    noisy_array = array + noise * array
    return noisy_array

# Time array
n_times = 90000
t_final = 30e-3 
time    = np.linspace(0,t_final,num=n_times)
tstep   = time[1]-time[0]   # Assuming a constant tstep

# Time evolution of Te[eV] and ne[m^-3]
Te    = np.zeros(n_times);   ne    = np.zeros(n_times)
Te[:] = 3;                  ne[:] = 1e20        # Replace with functions of time
Te =2e3* np.exp(-time/10e-3) + 10


# Information about charge states and initialization
n_charges          = n_states # Number of charge states including the neutral state (0)
ion_i_density      = np.zeros((n_charges, n_times))  # Higher charge index, more ionized (0=neutral)
ion_i_density_cor  = np.zeros((n_charges, n_times))  # Higher charge index, more ionized (0=neutral)
ion_i_density[:,0] = solve_ion_distribution(elem, min(ne[0],2e21), Te[0], rates_ion, rates_recom)
ion_i_density[:,0] = add_random_noise(ion_i_density[:,0]) 

# Ionization and recombination functions
def ionization_coefficient(z_ion, Te, ne):
    ioniz = 0.0
    if ((z_ion>= 0) and (z_ion<n_charges-1)):
        ioniz =rates_ion[z_ion](ne,Te) 
    return ioniz

def recombination_coefficient(z_ion, Te, ne):
    recomb = 0.0
    if ((z_ion>= 1) and (z_ion<=n_charges-1)):
        recomb = rates_recom[z_ion](ne,Te) 
    return recomb

# Time loop, evolve charge states
for i_time in range(n_times-1):

    ion_i_density_cor[:, i_time]  = solve_ion_distribution(elem, min(ne[i_time],2e21), Te[i_time], rates_ion, rates_recom)
    for z_ion in range(n_charges):

        # Get ionization and recombination coefficients from state i
        Sion_z = ionization_coefficient   (z_ion, Te[i_time], ne[i_time])
        Srec_z = recombination_coefficient(z_ion, Te[i_time], ne[i_time])
        
        # Get recombination coefficients from state i+1
        Srec_z_up = recombination_coefficient(z_ion+1, Te[i_time], ne[i_time])
        
        # Get ionization coefficients from state i-1
        Sion_z_dw = ionization_coefficient   (z_ion-1, Te[i_time], ne[i_time])
               
        ########## Evolve densities ##############################
        # Recombination and ionization losses
        RHS = -ion_i_density[z_ion, i_time] * (Sion_z + Srec_z) 
        
        # Gain from upper level recombination
        if ( z_ion+1 < n_charges ):
            RHS = RHS + ion_i_density[z_ion+1, i_time] * Srec_z_up

        # Gain from lower level ionization
        if ( z_ion-1 >= 0):
            RHS = RHS + ion_i_density[z_ion-1, i_time] * Sion_z_dw

        ion_i_density[z_ion, i_time+1] = ion_i_density[z_ion, i_time] + ne[i_time]*RHS*tstep 

# Get coronal equilibrium distribution
ion_i_density_cor[:, n_times-1]  = solve_ion_distribution(elem, min(ne[n_times-1],2e21), Te[n_times-1], rates_ion, rates_recom)

# Radiated power
Lrad = np.zeros(n_times)
for i_times in range(0, n_times):
    for j in range(0, n_states):
        if j==0:
            Lrad[i_times] = Lrad[i_times] + (rates_line[j](ne[i_times],Te[i_times]) ) * ion_i_density[j,i_times]
        elif j==n_states-1:
            Lrad[i_times] = Lrad[i_times] + (rates_contin[j](ne[i_times],Te[i_times])) * ion_i_density[j,i_times]
        else:
            Lrad[i_times] = Lrad[i_times] + (rates_line[j](ne[i_times],Te[i_times]) + rates_contin[j](ne[i_times],Te[i_times])) * ion_i_density[j,i_times]

Lrad_cor = np.zeros(n_times)
for i_times in range(0, n_times):
    for j in range(0, n_states):
        if j==0:
            Lrad_cor[i_times] = Lrad_cor[i_times] + (rates_line[j](ne[i_times],Te[i_times]) ) * ion_i_density_cor[j,i_times]
        elif j==n_states-1:
            Lrad_cor[i_times] = Lrad_cor[i_times] + (rates_contin[j](ne[i_times],Te[i_times])) * ion_i_density_cor[j,i_times]
        else:
            Lrad_cor[i_times] = Lrad_cor[i_times] + (rates_line[j](ne[i_times],Te[i_times]) + rates_contin[j](ne[i_times],Te[i_times])) * ion_i_density_cor[j,i_times]

 
# Plot densities as a function of time
cmap = get_cmap('hsv')
#cmap = get_cmap('rainbow')
colors = cmap(np.linspace(0, 1, n_charges))

tfact = 1e3
tplot = time*tfact
n_tot = np.zeros(n_times)
for i_time in range(n_times):
    n_tot[i_time] = sum(ion_i_density[:,i_time]) 

#for z_ion in range(n_charges):
#    plt.plot(tplot, ion_i_density[z_ion,:]/n_tot[:], color=colors[z_ion], label=f'$Z^{z_ion}$')  # Add label for legend
#    plt.plot(tplot, ion_i_density_cor[z_ion,:],'--', color=colors[z_ion])  # Add label for legend
#
#plt.plot(tplot, n_tot/n_tot[0], color='black', label=r'$n_{tot}$')  # Add label for legend
#plt.xlabel('Time [ms]')
#plt.ylabel('Fractional abundance')
#plt.yscale('log')
#plt.ylim(1e-5,2)
#plt.title('Charge state distribution')
#plt.grid(True)  # Add grid lines
#plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))  # Place legend outside the plot
#plt.show()



plt.plot(tplot, Lrad, color='black', label=r'Lrad')  # Add label for legend
plt.plot(tplot, Lrad_cor, '--', color='black', label=r'Lrad')  # Add label for legend
plt.xlabel('Time [ms]')
plt.ylabel('Radiation')
#plt.yscale('log')
#plt.ylim(1e-5,2)
plt.title('Charge state distribution')
plt.grid(True)  # Add grid lines
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))  # Place legend outside the plot
plt.show()
       
