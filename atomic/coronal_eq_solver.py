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


# Solve charge state distribution considering coronal equilibrium
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



# initialise the atomic data provider
adas = OpenADAS(permit_extrapolation=False)

#elem = tungsten 
elem = neon 

n_states = elem.atomic_number + 1 

# Collect rate coefficients
rates_ion   = get_rates_ionisation(elem)
rates_recom = get_rates_recombination(elem)
rates_line  = get_LT_radiation(elem)
rates_contin= get_BR_radiation(elem)

ne = 1e20  # electron density in m^-3
Te = 2   # electron temperature in eV

prob = solve_ion_distribution(elem, ne, Te, rates_ion, rates_recom)

print(prob)
