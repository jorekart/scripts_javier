import imas
import sys
sys.path.append('/home/ITER/artolaj/scripts_hub/scripts_javier/IMAS') 
from imas_custom_utils import parse_and_open_imas_entry
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
import numpy as np

vv = np.loadtxt('/home/ITER/artolaj/ITER_components/ITER_vessel_inner.txt')
fw = np.loadtxt('/home/ITER/artolaj/ITER_components/1st_wall.txt')

entry_and_args = parse_and_open_imas_entry() 

imas_entry = entry_and_args["entry"]
args       = entry_and_args["args"]

spi = imas_entry.get('spi')

it=0

# Initialize plot
fig, ax = plt.subplots()
plt.subplots_adjust(bottom=0.2)

ax.plot(vv[:, 0] * 1e-3, vv[:, 1] * 1e-3, label='Vessel')
ax.plot(fw[:, 0], fw[:, 1], label='First Wall')

# Plot dynamic SPI fragments
fragment_scatter = ax.scatter([], [], s=[], label='SPI Fragments')
ax.set_aspect('equal', adjustable='box')
ax.set_xlabel('R [m]')
ax.set_ylabel('Z [m]')
ax.legend()

def data_to_display_sizes(radius, ax):
    display_sizes = []
    for r in radius:
        size_in_data_coords = ax.transData.transform((r, 0)) - ax.transData.transform((0, 0))
        size_in_points = size_in_data_coords[0] *3e1
        display_sizes.append(size_in_points)
    return display_sizes

# Update function for the slider
def update_plot(val):
    global it
    it = int(slider.val)
    for injector in spi.injector:
        frag_R = []
        frag_Z = []
        frag_radius = []
        for frag in injector.fragment:
            frag_R.append(frag.position.r[it])
            frag_Z.append(frag.position.z[it])
            frag_rad = (frag.volume[it]*3/(4*np.pi))**(1/3.0)
            frag_radius.append(frag_rad) 
            
        display_sizes = data_to_display_sizes(frag_radius, ax)
        fragment_scatter.set_offsets(np.c_[frag_R, frag_Z])
        #fragment_scatter.set_sizes(np.array(frag_radius)*1e3)  # Size is proportional to the area
        fragment_scatter.set_sizes(display_sizes) 
    ax.set_title(f'SPI Fragments at Time={spi.time[it]*1e3} ms')
    fig.canvas.draw_idle()


# Add slider
ax_slider = plt.axes([0.2, 0.05, 0.6, 0.03], facecolor='lightgoldenrodyellow')
slider = Slider(ax_slider, 'Time Index', 0, len(spi.injector[0].fragment[0].position.r) - 1, valinit=it, valstep=1)
slider.on_changed(update_plot)

plt.show()