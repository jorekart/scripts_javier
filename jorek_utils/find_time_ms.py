# python find_time_ms.py time_in_ms t_norm_factor

import numpy as np
import sys

fname = './times.dat'

try:
    data = np.loadtxt(fname, skiprows=1)
except:
    print('File times.dat not found')
    quit()


tsearch = float(sys.argv[1])
tnorm   = float(sys.argv[2])

times  = data[:,1]*tnorm*1e3

diff      = abs( times - tsearch )

min_value = min(diff)
min_index = np.argmin(diff)
print('')
print('    Restart ' + str(tsearch) + ' corresponds to restart ' + str(int(data[min_index,0])) + ' with t = ' + str(data[min_index,0])  )
print('')

