import imas
import matplotlib.pyplot as plt
import sys
sys.path.append('/home/ITER/artolaj/scripts_hub/scripts_javier/IMAS') 
from imas_custom_utils import parse_and_open_imas_entry
import numpy as np
from scipy.interpolate import interp2d

entry_and_args = parse_and_open_imas_entry() 

imas_entry = entry_and_args["entry"]
args       = entry_and_args["args"]


wall_ids = imas_entry.get('wall')

for unit in wall_ids.description_2d[0].limiter.unit:
    R_wall = unit.outline.r
    Z_wall = unit.outline.z
    plt.plot(R_wall, Z_wall, label=unit.name)

plt.gca().set_aspect('equal', adjustable='box')
# Add labels and a colorbar

plt.legend()
plt.grid()
plt.xlabel('R [m]')
plt.ylabel('Z [m]')
plt.title('PFC contour')

# Show the plot
plt.show()

