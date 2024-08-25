import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from matplotlib.colors import Normalize

dat   = np.loadtxt('fort.78')

R = dat[:,0]
Z = dat[:,1]
T_max = dat[:,2]
T_melt = 3695

threshold = 100  # Replace with your desired threshold value

# Create boolean masks
mask_above = T_max > threshold
mask_below = ~mask_above  # Equivalent to f <= threshold

mask_melt = T_max > 3695

# Define the grid for interpolation
# grid_R, grid_Z = np.mgrid[min(R):max(R):100j, min(Z):max(Z):100j]

# # Interpolate T_max values onto the grid
# grid_T_max = griddata((R, Z), T_max, (grid_R, grid_Z), method='linear')

lev = np.linspace(-0.015,0.02,30)

# Normalize the function values to [0, 1] for controlling alpha
norm = Normalize(vmin=np.min(T_max[mask_above]), vmax=np.max(T_max[mask_above]))
f_normalized = norm(T_max[mask_above])

# Invert the normalized values to set alpha (lower values -> lower opacity)
alphas = f_normalized  # or use 1 - f_normalized if you want to reverse


# Create a contour plot
plt.figure()
#cp = plt.contourf(grid_x, grid_y, grid_T_max, levels=lev, cmap='nipy_spectral')
# cp = plt.contourf(grid_R, grid_Z, grid_T_max, levels=300, cmap='nipy_spectral')
cp = plt.scatter(R[mask_above], Z[mask_above], c=T_max[mask_above], s=1.6, cmap='jet', alpha=alphas, vmin=0, vmax=T_melt)
plt.scatter(R[mask_melt], Z[mask_melt], s=10, color='black')
plt.colorbar(cp)
plt.title(r'$\Delta T_{max}$ [K]')
plt.xlabel(r'$\phi$', fontsize=16)
plt.ylabel(r'$\theta$',fontsize=16)
plt.show()
