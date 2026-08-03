import numpy as np

ffp = np.loadtxt('ffp_ref')
pp  = np.loadtxt('t_ref')

factor = 1.00 

np.savetxt( 'jorek_ffprime_scaled',     np.transpose( [ffp[:,0], ffp[:,1]] ) )
np.savetxt( 'jorek_temperature_scaled', np.transpose( [ pp[:,0],  pp[:,1]] ) )
