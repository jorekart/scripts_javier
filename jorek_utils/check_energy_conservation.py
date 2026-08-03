import numpy as np
import sys
import matplotlib.pyplot as plt

############# Documentation #########################################
##   
##  Call the script with two arguments 
##   
##     python check_energy_conservation.py 0D.dat dEdt.dat 
##   
##  where 
##   
##     dEdt.dat:   Is the file that appears after calling 
##                 plot_live_data.sh -q dEdt
##   
##     0D.dat:     Is the file that is produced with jorek2_postproc
##                 with the command zeroD_quantities. See below an 
##                 example of the jorek2_postproc input file
##
##                     namelist input
##                     si-units
##                     for step 0 to 99999 do
##                     zeroD_quantities
##                     done
##   
##  IMPORTANT: You also need to provide the correct time normalization factor
##             t_norm below!
##
#####################################################################

zD_file = sys.argv[1]
dE_file = sys.argv[2]

zD_data = np.loadtxt(zD_file)
dE_data = np.loadtxt(dE_file, skiprows=1)

t_norm  = 6.4836e-7

show_tot = True
show_mag = True 
show_thm = True 
show_kin = True

# Hardcoded parameters, columns
i_dWtot_dt    = 2  - 1 
i_dWmag_dt    = 3  - 1 
i_dWthm_dt    = 4  - 1 
i_dWkinpar_dt = 6  - 1 
i_ohmic       = 48 - 1
i_mag_source  = 47 - 1
i_JxB_dot_v   = 37 - 1
i_poynting    = 58 - 1
i_thm_work    = 38 - 1 
i_ext_heat    = 42 - 1 
i_pvn_flux    = 53 - 1 
i_qpar_flux   = 54 - 1 
i_qprp_flux   = 55 - 1 
i_kinpar_flux = 57 - 1 
i_viscopar_diss = 45 - 1 
i_viscopar_flux = 59 - 1 

# Time arrays
time_zD = zD_data[:,0] * 1e3            # In ms
time_dE = dE_data[:,0] * 1e3 * t_norm   # In ms

# Magnetic energy conservation data
dWmag_dt    =  dE_data[:,i_dWmag_dt   ] * -1e-6
ohmic       =  zD_data[:,i_ohmic      ] * -1e-6
mag_source  =  zD_data[:,i_mag_source ] *  1e-6
JxB_dot_v   =  zD_data[:,i_JxB_dot_v  ] *  1e-6
poynting    =  zD_data[:,i_poynting   ] *  1e-6

total_terms_mag = ohmic + mag_source + JxB_dot_v + poynting 

# Thermal energy conservation data
dWthm_dt    =  dE_data[:,i_dWthm_dt   ] * -1e-6
ohmic       =  zD_data[:,i_ohmic      ] *  1e-6
thm_work    =  zD_data[:,i_thm_work   ] *  1e-6
ext_heat    =  zD_data[:,i_ext_heat   ] *  1e-6
pvn_flux    =  zD_data[:,i_pvn_flux   ] * -1e-6
qpar_flux   =  zD_data[:,i_qpar_flux  ] * -1e-6
qprp_flux   =  zD_data[:,i_qprp_flux  ] * -1e-6

total_terms_thm = ohmic + thm_work + ext_heat + pvn_flux + qpar_flux + qprp_flux 

# Parallel kinetic energy conservation data
dWkinpar_dt   =  dE_data[:,i_dWkinpar_dt   ] * -1e-6
thm_work      =  zD_data[:,i_thm_work      ] * -1e-6
viscopar_diss =  zD_data[:,i_viscopar_diss ] * -1e-6 #* 0 
viscopar_flux =  zD_data[:,i_viscopar_flux ] *  1e-6
kinpar_flux   =  zD_data[:,i_kinpar_flux   ] * -1e-6

total_terms_kinpar = thm_work + viscopar_diss + viscopar_flux + kinpar_flux 

# total conservation
total_terms = mag_source + poynting + ext_heat + pvn_flux + qpar_flux + qprp_flux + viscopar_diss + viscopar_flux + kinpar_flux
dWtot_dt    =  dE_data[:,i_dWtot_dt   ] * -1e-6

# Plot total conservation 
if (show_tot):
  fig = plt.figure()
  ax = fig.add_axes([0.1, 0.1, 0.6, 0.8]) # main axes
  
  ax.set_title("Total energy conservation")
  ax.set_xlabel('Time [ms]', fontsize=16)
  ax.set_ylabel('Powers [MW]', fontsize=16)
  ax.grid()
  
  ax.plot(time_dE, dWtot_dt, '--',  label=r'$dW_{tot}/dt$',              color='black')
  ax.plot(time_zD, total_terms,     label='Sum of terms',                color='black')
  ax.legend(fontsize=12, loc='center left', bbox_to_anchor=(1.0, 0.5))
  
  plt.show()
  #fig.savefig('mag_energy_conservation.pdf', bbox_inches='tight')




# Plot mag energy conservation
if (show_mag):
  fig = plt.figure()
  ax = fig.add_axes([0.1, 0.1, 0.6, 0.8]) # main axes
  
  ax.set_title("Magnetic energy conservation")
  ax.set_xlabel('Time [ms]', fontsize=16)
  ax.set_ylabel('Powers [MW]', fontsize=16)
  ax.grid()
  
  ax.plot(time_dE, dWmag_dt, '--',  label=r'$dW_{mag}/dt$',              color='black')
  ax.plot(time_zD, total_terms_mag, label='Sum of terms',                color='black')
  ax.plot(time_zD,           ohmic, label=r'-$\eta J^2$',                 color='blue')
  ax.plot(time_zD,      mag_source, label=r' $\eta J\cdot J_{source}$',   color='pink')
  ax.plot(time_zD,       JxB_dot_v, label=r'-$J\times B\cdot v$',         color='green')
  ax.plot(time_zD,        poynting, label=r'-$E\times B\cdot n$',         color='red')
  ax.legend(fontsize=12, loc='center left', bbox_to_anchor=(1.0, 0.5))
  
  plt.show()
  #fig.savefig('mag_energy_conservation.pdf', bbox_inches='tight')


# Plot thermal energy conservation
if (show_thm):
  fig = plt.figure()
  ax = fig.add_axes([0.1, 0.1, 0.6, 0.8]) # main axes
  
  ax.set_title("Thermal energy conservation")
  ax.set_xlabel('Time [ms]', fontsize=16)
  ax.set_ylabel('Powers [MW]', fontsize=16)
  ax.grid()
  
  ax.plot(time_dE, dWthm_dt, '--',  label=r'$dW_{th}/dt$',               color='black')
  ax.plot(time_zD, total_terms_thm, label='Sum of terms',                color='black')
  ax.plot(time_zD,           ohmic, label=r' $\eta J^2$',                color='blue')
  ax.plot(time_zD,        ext_heat, label=r' $S_{ext}$',                 color='pink')
  ax.plot(time_zD,        thm_work, label=r' $v_\parallel \cdot \nabla p$',    color='orange')
  ax.plot(time_zD,       qpar_flux, label=r'-$q_\parallel \cdot n$',        color='red')
  ax.plot(time_zD,       qprp_flux, label=r'-$q_\perp\cdot n$',        color='green')
  ax.plot(time_zD,        pvn_flux, label=r'-$\gamma p v_\parallel\cdot n$',        color='brown')
  ax.legend(fontsize=12, loc='center left', bbox_to_anchor=(1.0, 0.5))
  
  plt.show()
  #fig.savefig('mag_energy_conservation.pdf', bbox_inches='tight')


# Plot thermal energy conservation
if (show_kin):
  fig = plt.figure()
  ax = fig.add_axes([0.1, 0.1, 0.6, 0.8]) # main axes
  
  ax.set_title("Kinetic energy conservation")
  ax.set_xlabel('Time [ms]', fontsize=16)
  ax.set_ylabel('Powers [MW]', fontsize=16)
  ax.grid()
  
  ax.plot(time_dE, dWkinpar_dt, '--',  label=r'$dW_{kin,\parallel}/dt$',               color='black')
  ax.plot(time_zD, total_terms_kinpar, label='Sum of terms',                color='black')
  ax.plot(time_zD,        thm_work, label=r'-$v_\parallel \cdot \nabla p$',    color='orange')
  ax.plot(time_zD,        kinpar_flux, label=r'$-\frac{1}{2}\rho v^2 v\cdot n$',        color='red')
  ax.plot(time_zD,       viscopar_flux, label=r'$\nu_\parallel\nabla v_\parallel\cdot n$',        color='blue')
  ax.plot(time_zD,       viscopar_diss, label=r'$-\nu_\parallel\nabla v_\parallel\cdot\nabla v_\parallel$', color='green')
  ax.legend(fontsize=12, loc='center left', bbox_to_anchor=(1.0, 0.5))
  
  plt.show()
  #fig.savefig('mag_energy_conservation.pdf', bbox_inches='tight')




