import imas
import matplotlib.pyplot as plt
import sys
sys.path.append('/home/ITER/artolaj/scripts_hub/scripts_javier/IMAS') 
from imas_custom_utils import parse_and_open_imas_entry
import numpy as np
from scipy.interpolate import interp2d

vv = np.loadtxt('/home/ITER/artolaj/ITER_components/ITER_vessel_inner.txt')
fw = np.loadtxt('/home/ITER/artolaj/ITER_components/1st_wall.txt')

entry_and_args = parse_and_open_imas_entry() 

imas_entry = entry_and_args["entry"]
args       = entry_and_args["args"]

qtty = args["quantity"]
it   = args["time_index"]

idslist = {}

idslist['equilibrium'] = imas_entry.get('equilibrium')

r = idslist['equilibrium'].time_slice[it].ggd[0].r[0].values
z = idslist['equilibrium'].time_slice[it].ggd[0].z[0].values
Ic =idslist['equilibrium'].time_slice[it].ggd[0].j_phi[0].values


# psi_bnd = idslist['equilibrium'].time_slice[it].boundary.psi
# psi_sep = idslist['equilibrium'].time_slice[it].boundary_separatrix.psi
# psi_sep2 = idslist['equilibrium'].time_slice[it].boundary_secondary_separatrix.psi
# psi_axis = idslist['equilibrium'].time_slice[it].global_quantities.psi_axis

# psin = np.linspace(1.0, 1.1, num=40)
#contour_levels = sorted(psi_axis + psin * (psi_bnd - psi_axis ))

# Create a contour plot
#contour_plot = plt.contour(x, y, psi2d, levels=contour_levels, cmap='viridis')
#contour_plot = plt.contour(x, y, psi2d, cmap='viridis', levels=100)
plt.figure(figsize=(8, 6))
scatter = plt.scatter(r, z, c=Ic, cmap='viridis', s=10)  # s controls the marker size
plt.colorbar(scatter, label='Ic [units]')  # Add colorbar to indicate Ic values
plt.plot(vv[:,0]*1e-3,vv[:,1]*1e-3)
plt.plot(fw[:,0],fw[:,1])
print(idslist['equilibrium'].time[it])
plt.gca().set_aspect('equal', adjustable='box')
# Add labels and a colorbar

print(np.sum(Ic))

plt.xlabel('R [m]')
plt.ylabel('Z [m]')
plt.title('I_fila')
#plt.colorbar(contour_plot, label='Psi')

# Show the plot
plt.show()
