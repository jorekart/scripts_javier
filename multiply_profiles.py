import numpy as np
from scipy import interpolate

prof1 = np.loadtxt('prof1')
prof2 = np.loadtxt('prof2')

x1 = prof1[:,0]

p1_interp = interpolate.interp1d(prof1[:,0], prof1[:,1])
p2_interp = interpolate.interp1d(prof2[:,0], prof2[:,1])

mult = p1_interp(x1) * p2_interp(x1)

np.savetxt( 'multiplied',     np.transpose( [x1, mult] ) )
