import matplotlib
import os
import imas 
import matplotlib.pyplot as plt
import numpy as np
from imas import imasdef
import argparse
import getpass

print(" ")
print(" Example of usage: ")
print("    python read_prof.py -u public -d ITER -s 105033 -r 1")
print(" ")
print(" To see options and default values do: ")
print("    python read_prof.py -h ")
print(" ")

# Import shot
parser = argparse.ArgumentParser(description="Plot profiles from IMAS databases",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-s", "--shot", type=int, default=170, help="Shot number")
parser.add_argument("-r", "--run", type=int, default=1, help="Run number")
parser.add_argument("-u", "--user", type=str, default=getpass.getuser(),
                    help="Location of ~$USER/public/imasdb")
parser.add_argument("-d", "--database", type=str, default="test", help="Database name under public/imasdb/")
parser.add_argument("-o", "--occurrence", type=int, default=0, help="Occurrence number")
parser.add_argument("-p", "--profile", type=str, default="Te", help="Profile name, e.g. -p Te")
parser.add_argument("-f", "--backend", type=int, default=imasdef.MDSPLUS_BACKEND,
                    help="Database format: 12=MDSPLUS, 13=HDF5")
#parser.add_argument("-t", "--time", type=float, default=-1, help="The requested time in seconds")
args = parser.parse_args()

# Available profiles
class profile:
    def __init__(self, name, ids, label, path):
        self.name  = name
        self.ids   = ids
        self.label = label
        self.path  = path
        self.data  = input.partial_get(ids, path, occurrence=0)

# Open input datafile
shot,run,user,database =args.shot,args.run,args.user,args.database
input = imas.DBEntry(imas.imasdef.MDSPLUS_BACKEND,database,shot,run,user,data_version='3')
input.open()


Te_name   = 'Te'
Jtor_name = 'Jtor'
ne_name   = 'ne'
psi_name  = 'psi'
transform_to_psiN = True

if (args.profile==Te_name): 
    prof  = profile(Te_name,'core_profiles',r'$T_e$ [eV]', 'profiles_1d(:)/electrons/temperature')
elif (args.profile==Jtor_name): 
    prof  = profile(Jtor_name,'core_profiles',r'$J_{tor}$ [A/m$^2$]', 'profiles_1d(:)/j_ohmic')
elif (args.profile==ne_name): 
    prof  = profile(ne_name,'core_profiles',r'$n_e$ [m$^{-3}$]', 'profiles_1d(:)/electrons/density')
elif (args.profile==psi_name):
    prof = profile('psi','core_profiles',r'$\psi$ [Wb]', 'profiles_1d(:)/grid/psi')
    
n_prof  = len(prof.data[0,:])

if (prof.ids=='core_profiles'):

    # coord   = profile('psi','core_profiles',r'$\psi$ [Wb]', 'profiles_1d(:)/grid/psi')
    
    # if transform_to_psiN:
    #     psi_ax  = -input.partial_get('equilibrium', 'time_slice(:)/global_quantities/psi_axis', occurrence=0) 
    #     psi_bnd = -input.partial_get('equilibrium', 'time_slice(:)/boundary/psi', occurrence=0) 
    #     n_rad   = len(coord.data[:,0])
    #     for i in range(0,n_rad):
    #         coord.data[i,:] = (coord.data[i,:] - psi_ax) / (psi_bnd - psi_ax)  # Sign is wrong!! for DINA!!
    #     coord.name = 'PsiN'
    #     coord.label= r'$\psi_N$'
        
    coord = profile('rho_tor','core_profiles',r'$\rho_{norm}$', 'profiles_1d(:)/grid/rho_tor_norm')

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
