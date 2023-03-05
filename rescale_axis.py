import numpy as np
import sys

input_file = sys.argv[1]
index      = int(sys.argv[2])
factor     = float(sys.argv[3])

data = np.loadtxt(input_file,skiprows=1)

data[:, index] = data[:,index]*factor

np.savetxt(input_file, data  )
