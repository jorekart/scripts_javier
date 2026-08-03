import numpy as np
import sys

fname = './times.dat'

try:
    data = np.loadtxt(fname, skiprows=1)
except:
    print('File times.dat not found')
    quit()


rsearch = float(sys.argv[1])

restarts  = data[:,1]

diff      = abs( restarts - rsearch )

min_value = min(diff)
min_index = np.argmin(diff)
print('')
print('    Restart ' + str(rsearch) + ' corresponds to restart ' + str(int(data[min_index,1])) + ' with t = ' + str(data[min_index,0])  )
print('')

