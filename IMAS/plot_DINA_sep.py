import imas
import matplotlib.pyplot as plt
import sys
import numpy as np
from imas import imasdef

uri1 = "imas:mdsplus?user=artolaj;pulse=100033;run=2;database=ITER_DISRUPTIONS;version=3"
uri2 = "imas:mdsplus?user=artolaj;pulse=100314;run=2;database=ITER_DISRUPTIONS;version=3"
uri3 = "imas:mdsplus?user=artolaj;pulse=100187;run=2;database=ITER_DISRUPTIONS;version=3"


uri_list = [uri1,uri2,uri3]
colors = ["blue", "black", "red"]

Zaxis_compare = [1,2,3,4]

for i, uri in enumerate(uri_list):
    imas_entry = imas.DBEntry(uri=uri, mode='r')
    imas_entry.open()
    
    summary = imas_entry.get('summary')
    
    Raxis = summary.boundary.magnetic_axis_r.value
    Zaxis = summary.boundary.magnetic_axis_z.value
    times = summary.time
    
    closest_indices = [np.argmin(np.abs(Zaxis - target)) for target in Zaxis_compare]
    
    for isep in closest_indices:
        eq = imas_entry.get_slice("equilibrium",times[isep],imasdef.PREVIOUS_INTERP)
        r_out = eq.time_slice[0].boundary.outline.r
        z_out = eq.time_slice[0].boundary.outline.z
        plt.plot(r_out, z_out, color=colors[i])
    
    plt.plot(Raxis, Zaxis, '--', color=colors[i])

    
    imas_entry.close()

uri_wall = "imas:mdsplus?user=public;pulse=116000;run=5;database=ITER_MD;version=3"
imas_entry = imas.DBEntry(uri=uri_wall, mode='r')
wall_ids = imas_entry.get('wall')

for unit in wall_ids.description_2d[0].limiter.unit:
    R_wall = unit.outline.r
    Z_wall = unit.outline.z
    #plt.plot(R_wall, Z_wall, label=unit.name)
    plt.plot(R_wall, Z_wall, label='', color='darkgreen')

plt.gca().set_aspect('equal', adjustable='box')
# Add labels and a colorbar

plt.legend()
plt.grid()
plt.xlabel('R [m]')
plt.ylabel('Z [m]')
#plt.title('PFC contour')

# Show the plot
plt.show()

