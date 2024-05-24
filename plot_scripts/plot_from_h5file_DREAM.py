import numpy as np
import h5py
import sys
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.widgets import Slider, Button
matplotlib.use("TkAgg")

if len(sys.argv) < 2:
    print("Usage: python plot_from_h5file_DREAM.py output_file_name.h5")
    sys.exit(1)

DREAM_file = sys.argv[1]

f = h5py.File(DREAM_file, 'r')

time  = f["/grid/t"][:]
rad   = f["/grid/r"][:]
Ip    = f["/eqsys/I_p"][:]
j_ohm = f["/eqsys/j_ohm"][:,:]   # ohmic current density in A/m²
j_tot = f["/eqsys/j_tot"][:,:]   # total current density in A/m²
Te    = f["/eqsys/T_cold"][:,:]  # electron temperature in eV
ne    = f["/eqsys/n_cold"][:,:]   # electron density in m^-3 

#  dataset    /eqsys/E_field
#  dataset    /eqsys/I_p
#  dataset    /eqsys/I_wall
#  dataset    /eqsys/T_cold
#  dataset    /eqsys/V_loop_w
#  dataset    /eqsys/W_cold
#  dataset    /eqsys/Y_p
#  dataset    /eqsys/j_hot
#  dataset    /eqsys/j_ohm
#  dataset    /eqsys/j_re
#  dataset    /eqsys/j_tot
#  dataset    /eqsys/n_cold
#  dataset    /eqsys/n_hot
#  dataset    /eqsys/n_i
#  dataset    /eqsys/n_re
#  dataset    /eqsys/n_tot
#  dataset    /eqsys/psi_edge
#  dataset    /eqsys/psi_p
#  dataset    /eqsys/psi_wall

prof0 = j_ohm
prof1 = j_tot
nt    = len(time)

# Plot profile with slider for time
fig, ax = plt.subplots()
plt.subplots_adjust(bottom=0.35)

l0, = plt.plot(rad, prof0[0,:],label='J_ohm') 
l1, = plt.plot(rad, prof1[0,:],label='J_tot') 
plt.legend()

allowed_steps = np.arange(nt)

axtime = plt.axes([0.25, 0.15, 0.65, 0.03])
profile = Slider(
    axtime, "Time [ms]", 0, nt,
    valinit=0, valstep=allowed_steps,
    color="green")
 

time_labels = [str(time[i]*1e3) for i in range(nt)]
profile.valtext.set_text(time_labels[0])  # Initial time label

def update(val):
    f = profile.val
    l0.set_ydata(prof0[f,:])
    l1.set_ydata(prof1[f,:])
    profile.valtext.set_text(time_labels[f])  # Update time label
 
profile.on_changed(update)



plt.show()
