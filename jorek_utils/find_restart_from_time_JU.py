# This scripts takes a time in JOREK units given in the command line and 
# outputs the corresponding restart file
import numpy as np
import sys

fname = './times.dat'

try:
    data = np.loadtxt(fname, skiprows=1)
except:
    print('File times.dat not found')
    quit()


tsearch = float(sys.argv[1])

times     = data[:,1]

diff      = abs( times - tsearch )

min_value = min(diff)
min_index = np.argmin(diff)
print('')
print('    t = ' + str(tsearch) + ' corresponds to restart file ' + str(int(data[min_index,0])) + ' with t = ' + str(data[min_index,1])  )
print('')

