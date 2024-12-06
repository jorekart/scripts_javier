import imas
import matplotlib.pyplot as plt
import sys
sys.path.append('/home/ITER/artolaj/scripts_hub/scripts_javier/IMAS') 
from imas_custom_utils import parse_and_open_imas_entry
import numpy as np

entry_and_args = parse_and_open_imas_entry() 

imas_entry = entry_and_args["entry"]
args       = entry_and_args["args"]

qtty = args["quantity"]

summary = imas_entry.get('summary')

time = summary.time*1e3

if (qtty=="Ip"):
    value = summary.global_quantities.ip.value * 1e-6
    qtty_name = 'Ip [MA]'
elif (qtty=="Zaxis"):
    value = summary.local.magnetic_axis.position.z
    qtty_name = 'Zaxis [m]'
elif (qtty=="q95"):
    value =  summary.global_quantities.q_95.value 
    qtty_name = 'q95 [-]'
elif (qtty=="Wth"):
    value =  summary.global_quantities.energy_thermal.value 
    qtty_name = 'W_th [J]'
elif (qtty=="Te_av"):
    value =  summary.volume_average.t_e.value 
    qtty_name = 'Te_av [eV]'
elif (qtty=="ne_av"):
    value =  summary.volume_average.n_e.value 
    qtty_name = 'ne_av [eV]'

      
plt.plot(time,value,'*-')

plt.ylabel(qtty_name)
plt.xlabel('Time [ms]')

# Show the plot
plt.show()