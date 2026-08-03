import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from matplotlib.colors import Normalize

dat   = np.loadtxt('fort.78')
#fw    = np.loadtxt('wall.txt')

R = dat[:,0]
Z = dat[:,1]
T_max = dat[:,2] 
T_melt = 1556
#for i in range(len(T_max)):
#  if (T_max[i]>0):
#    T_max[i] = T_max[i] + 473.15
#Rw = fw[:,0]
#Zw = fw[:,1]

#thw = np.arctan2(-Zw,Rw-2.8)

threshold = 473.15 # Replace with your desired threshold value

# Create boolean masks
mask_above = T_max > threshold
mask_melt  = T_max > T_melt

Rf = R[mask_above]
Zf = Z[mask_above]
Tf = T_max[mask_above]

# Normalize the function values to [0, 1] for controlling alpha
norm = Normalize(vmin=np.min(Tf), vmax=np.max(Tf))
f_normalized = norm(Tf)

# Invert the normalized values to set alpha (lower values -> lower opacity)
alphas = f_normalized  # or use 1 - f_normalized if you want to reverse


# Create a contour plot
plt.figure()

cp = plt.scatter(Rf, Zf, c=Tf, s=1.6, cmap='jet', alpha=alphas, vmin=threshold, vmax=T_melt)

plt.scatter(R[mask_melt], Z[mask_melt], s=10, color='black')


plt.colorbar(cp)
plt.xlim(-np.pi*1.15,np.pi*1.15)
plt.ylim(-3,-0.5)
plt.title(r'$T_{max}$ [K]',fontsize=22)
plt.xlabel(r'$\phi$', fontsize=16)
plt.ylabel(r'$\theta$',fontsize=16)
plt.show()
