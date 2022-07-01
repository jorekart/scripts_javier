import numpy as np
import sys

input_file = sys.argv[1]
n_times    = int(sys.argv[2])

try:
    data = np.loadtxt(input_file, skiprows=1)
except:
    print('Input file not found')
    quit()

times_vec    = data[:,0]
restarts_vec = data[:,1]

t0  = times_vec[0]    
tf  = times_vec[-1]

times    = np.linspace(t0, tf, num=n_times)
restarts = np.zeros(n_times, dtype=int) 
times_f  = np.zeros(n_times) 

for i in range(0,len(times)):

    diff      = abs( times_vec - times[i] )
    min_value = min(diff)
    min_index = np.argmin(diff)

    restarts[i] = int(restarts_vec[min_index])
    times_f[i]   = times_vec[min_index]*1000 # In ms

 

np.savetxt( 'times_restarts.txt',  np.transpose( [times_f[:], restarts[:]] ), fmt='%.1f %05d' )
