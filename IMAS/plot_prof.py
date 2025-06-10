import imas 
import matplotlib.pyplot as plt
import numpy as np
from imas import imasdef
import sys
sys.path.append('/home/ITER/artolaj/scripts_hub/scripts_javier/IMAS') 
from imas_custom_utils import parse_and_open_imas_entry
import numpy as np

entry_and_args = parse_and_open_imas_entry() 

input = entry_and_args["entry"]
args  = entry_and_args["args"]

qtty = args["quantity"]

# Available profiles
class profile:
    def __init__(self, name, ids, label, path):
        self.name  = name
        self.ids   = ids
        self.label = label
        self.path  = path
        self.data  = input.partial_get(ids, path, occurrence=0)

Te_name   = 'Te'
Jtor_name = 'Jphi'
ne_name   = 'ne'
psi_name  = 'psi'
transform_to_psiN = False
profile_ids = 'plasma_profiles'

if (int(args["version"]) == 3):
    profile_ids = 'core_profiles'

if (qtty==Te_name): 
    prof  = profile(Te_name,profile_ids,r'$T_e$ [eV]', 'profiles_1d(:)/electrons/temperature')
elif (qtty==Jtor_name): 
    prof  = profile(Jtor_name,profile_ids,r'$J_{tor}$ [A/m$^2$]', 'profiles_1d(:)/j_ohmic')
elif (qtty==ne_name): 
    prof  = profile(ne_name,profile_ids,r'$n_e$ [m$^{-3}$]', 'profiles_1d(:)/electrons/density')
elif (qtty==psi_name):
    prof = profile('psi',profile_ids,r'$\psi$ [Wb]', 'profiles_1d(:)/grid/psi')
    
n_prof  = len(prof.data[0,:])

if (prof.ids==profile_ids):

    # coord   = profile('psi','core_profiles',r'$\psi$ [Wb]', 'profiles_1d(:)/grid/psi')
    
    # if transform_to_psiN:
    #     psi_ax  = -input.partial_get('equilibrium', 'time_slice(:)/global_quantities/psi_axis', occurrence=0) 
    #     psi_bnd = -input.partial_get('equilibrium', 'time_slice(:)/boundary/psi', occurrence=0) 
    #     n_rad   = len(coord.data[:,0])
    #     for i in range(0,n_rad):
    #         coord.data[i,:] = (coord.data[i,:] - psi_ax) / (psi_bnd - psi_ax)  # Sign is wrong!! for DINA!!
    #     coord.name = 'PsiN'
    #     coord.label= r'$\psi_N$'
        
    coord = profile('rho_tor',profile_ids,r'$\rho_{norm}$', 'profiles_1d(:)/grid/rho_tor_norm')

time  = input.partial_get('equilibrium', 'time(:)', occurrence=0) 

input.close()

plt.xlabel(coord.label)
plt.ylabel(prof.label)

jump   = int(n_prof/10)

for ind in range(0,n_prof,jump):

    plt.plot(coord.data[:,ind],prof.data[:,ind], '-',label=f'{time[ind]:.2e}')

plt.legend(title='Time [s]')
plt.grid()
plt.show()
