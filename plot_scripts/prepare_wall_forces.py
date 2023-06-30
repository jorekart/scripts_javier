import numpy as np
import sys

input_file = sys.argv[1]

data = np.loadtxt(input_file)

Fx = data[:,3] 
Fy = data[:,4] 

Fh = np.sqrt( Fx**2 + Fy**2)
phase = np.arctan2(Fy,Fx)

nrow = np.shape(data)[0]
ncol = np.shape(data)[1]

data2 = np.zeros((nrow,ncol+2))

data2[:,0:ncol] = data
data2[:,ncol+0] = Fh
data2[:,ncol+1] = phase

np.savetxt(input_file, data2  )
