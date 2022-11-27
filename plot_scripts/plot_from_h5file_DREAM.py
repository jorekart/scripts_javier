import numpy as np
import h5py
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.widgets import Slider, Button
matplotlib.use("TkAgg")

f = h5py.File('output_restart_CQ.h5', 'r')

time  = f["/grid/t"][:]
rad   = f["/grid/r"][:]
Ip    = f["/eqsys/I_p"][:]
j_ohm = f["/eqsys/j_ohm"][:,:]
j_tot = f["/eqsys/j_tot"][:,:]

#plt.plot(rad, j_ohm[0,:])
#plt.show()

prof0 = j_ohm
prof1 = j_tot
nt    = len(time)

# Plot profile with slider for time
fig, ax = plt.subplots()
plt.subplots_adjust(bottom=0.35)

l0, = plt.plot(rad, prof0[0,:]) 
l1, = plt.plot(rad, prof1[0,:]) 

allowed_steps = np.arange(nt)

axtime = plt.axes([0.25, 0.15, 0.65, 0.03])
freq = Slider(
    axtime, "Time", 0, nt,
    valinit=0, valstep=allowed_steps,
    color="green")
 
def update(val):
    f = freq.val
    l0.set_ydata(prof0[f,:])
    l1.set_ydata(prof1[f,:])
 
freq.on_changed(update)

plt.show()
