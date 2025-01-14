# -----------------------------------------
# CREATE DICTIONARIES FROM ASCII DATAFILES
# GENERATED BY DINA DISRUPTION SIMULATOR
# -----------------------------------------
import numpy as np
import string

# ------------------------------------------------------------------------------------------------

def dictionary_from_file(filename):
  ascii_file = open(filename,'r')
  data_from_file = ascii_file.read()
  ascii_file.close()
  # ALL THE LINES
  file_lines = data_from_file.split('\n')
  # NAME OF VARIABLES (FIRST LINE)
  file_list_of_names = file_lines[0]
  file_variables = file_list_of_names.split()
  # DATA ITSELF: EACH LINE / INDEX CORRESPONDS TO VALUES OF VARIABLES FOR A TIME SLICE
  file_data_lines = file_lines[1:]
  file_data = [None]*(len(file_data_lines)-1)
  for i in range(len(file_data_lines)-1):
    file_data[i] = file_data_lines[i].split()
  # CREATE DICTIONARY: KEYS = VARIBALES NAMES, VALUES = DATA ITSELF
  file_dict = {}
  for itime in range(len(file_data)):
    for jvar in range(len(file_variables)):
      if jvar==0:
        file_dict[itime] = dict.fromkeys([file_variables[jvar]],file_data[itime][jvar])
      else:
        file_dict[itime].update(dict.fromkeys([file_variables[jvar]],file_data[itime][jvar]))
  return file_dict

# ------------------------------------------------------------------------------------------------

def dictionary_from_dens(file1,file2,file3):

  ascii_file1 = open(file1,'r')
  data_from_file1 = ascii_file1.read()
  ascii_file1.close()

  ascii_file2 = open(file2,'r')
  data_from_file2 = ascii_file2.read()
  ascii_file2.close()

  ascii_file3 = open(file3,'r')
  data_from_file3 = ascii_file3.read()
  ascii_file3.close()

  file_lines1 = data_from_file1.split('\n')
  file_lines2 = data_from_file2.split('\n')
  file_lines3 = data_from_file3.split('\n')

  ntime = int(file_lines1[1])
  file_dict1={}
  file_dict1['z_ion'] = 1
  file_dict1['Time(s)'] = []
  file_dict1['Density(m-3)'] = []
  for itime in range(ntime):
    file_dict1['Time(s)'].append(float(file_lines1[3+itime].split()[0]))
    file_dict1['Density(m-3)'].append(float(file_lines1[3+itime].split()[1])*1e19)

  ntime = int(file_lines2[1].split()[0])
  file_dict2={}
  file_dict2['z_ion'] = int(file_lines2[1].split()[1])
  file_dict2['Time(s)'] = []
  file_dict2['Density(m-3)'] = []
  for itime in range(ntime):
    file_dict2['Time(s)'].append(float(file_lines2[3+itime].split()[0]))
    file_dict2['Density(m-3)'].append(float(file_lines2[3+itime].split()[1])*1e19)

  ntime = int(file_lines3[1].split()[0])
  file_dict3={}
  file_dict3['z_ion'] = int(file_lines3[1].split()[1])
  file_dict3['Time(s)'] = []
  file_dict3['Density(m-3)'] = []
  for itime in range(ntime):
    file_dict3['Time(s)'].append(float(file_lines3[3+itime].split()[0]))
    file_dict3['Density(m-3)'].append(float(file_lines3[3+itime].split()[1])*1e19)

  return file_dict1,file_dict2,file_dict3


# ------------------------------------------------------------------------------------------------

def dictionary_from_psi_data(filename):
  ascii_file = open(filename,'r')
  data_from_file = ascii_file.read()
  ascii_file.close()
  # ALL THE LINES
  file_lines = data_from_file.split('\n')
  nc = int(np.asfarray(file_lines[0].split())[0])
  data_all_time_slices = ''.join(file_lines).split(' '+str(nc)+' ')[1:]
  ntime = len(data_all_time_slices)
  nr = 33;
  nz = 65;
  file_dict = {}
  psi_fact  = -1e-5 * 2*np.pi;  size_fact = 1e-2;   time_fact = 1e-3; curr_fact = -1e7
  for itime in range(ntime):
    #print('itime = ',str(itime),'/',str(ntime))
    data_slice = np.asfarray(data_all_time_slices[itime].split(),float)
    ncontour = nc
    nrzc  = int(data_slice[0])
    nbndy = int(data_slice[2])
    ncurd = int(data_slice[3])
    data_slice = data_slice[4:]

    i_start = 0         
    i_end   = i_start + nrzc
    rc      = data_slice[i_start:i_end]
    
    i_start = i_end 
    i_end   = i_start + nrzc
    zc      = data_slice[i_start:i_end]
    
    i_start = i_end     
    i_end   = i_start + ncontour
    xu      = data_slice[i_start:i_end]*size_fact            # [M]
    
    i_start = i_end
    i_end   = i_start + ncontour
    yu      = data_slice[i_start:i_end]*size_fact            # [M]
    
    i_start  = i_end
    i_end    = i_start + 8
    scalar_data = data_slice[i_start:i_end]
    
    delta_R_for_psi_grid = scalar_data[0]*size_fact              # [M]
    delta_Z_for_psi_grid = scalar_data[1]*size_fact              # [M]
    poloidal_flux_at_magnetic_axis   = scalar_data[2]*psi_fact   # [Wb]
    poloidal_flux_at_plasma_boundary = scalar_data[3]*psi_fact   # [Wb]
    poloidal_flux_at_halo_boundary   = scalar_data[4]*psi_fact   # [Wb]
    R0_magnetic = scalar_data[5]*size_fact                       # [M]
    Z0_magnetic = scalar_data[6]*size_fact                       # [M]
    Time = scalar_data[7]*time_fact                              # [S]
    
    i_start  = i_end
    i_end    = i_start + nr*nz
    poloidal_flux_array = data_slice[i_start:i_end]*psi_fact     # [Wb]

    i_start  = i_end
    i_end    = i_start + nr
    R_of_the_flux = data_slice[i_start:i_end]*size_fact        # [M]

    i_start  = i_end
    i_end    = i_start + nz
    Z_of_the_flux = data_slice[i_start:i_end]*size_fact  # [M]

    i_start  = i_end
    i_end    = i_start + nbndy
    R_of_boundary_points = data_slice[i_start:i_end]*size_fact    # [M]
    
    i_start  = i_end
    i_end    = i_start + nbndy
    Z_of_boundary_points = data_slice[i_start:i_end]*size_fact    # [M]

    i_start  = i_end
    i_end    = i_start + ncurd
    R_of_cur_density = data_slice[i_start:i_end]*size_fact       # [M]

    i_start  = i_end
    i_end    = i_start + ncurd
    current_density = data_slice[i_start:i_end]*curr_fact      # [A/M2]

    file_dict[itime] = dict.fromkeys(['R_surrounding_limiter_contour[m]'],xu)
    file_dict[itime].update(dict.fromkeys(['nr'],nr))
    file_dict[itime].update(dict.fromkeys(['nz'],nz))
    file_dict[itime].update(dict.fromkeys(['Z_surrounding_limiter_contour[m]'],yu))
    file_dict[itime].update(dict.fromkeys(['delta_R_for_psi_grid[m]'],delta_R_for_psi_grid))
    file_dict[itime].update(dict.fromkeys(['delta_Z_for_psi_grid[m]'],delta_Z_for_psi_grid))
    file_dict[itime].update(dict.fromkeys(['poloidal_flux_at_magnetic_axis[Wb]'],poloidal_flux_at_magnetic_axis))
    file_dict[itime].update(dict.fromkeys(['poloidal_flux_at_plasma_boundary[Wb]'],poloidal_flux_at_plasma_boundary))
    file_dict[itime].update(dict.fromkeys(['poloidal_flux_at_halo_boundary[Wb]'],poloidal_flux_at_halo_boundary))
    file_dict[itime].update(dict.fromkeys(['R0_magnetic[m]'],R0_magnetic))
    file_dict[itime].update(dict.fromkeys(['Z0_magnetic[m]'],Z0_magnetic))
    file_dict[itime].update(dict.fromkeys(['Time[s]'],Time))
    file_dict[itime].update(dict.fromkeys(['poloidal_flux_array[Wb]'],poloidal_flux_array))
    file_dict[itime].update(dict.fromkeys(['R_of_the_flux[m]'],R_of_the_flux))
    file_dict[itime].update(dict.fromkeys(['Z_of_the_flux[m]'],Z_of_the_flux))
    file_dict[itime].update(dict.fromkeys(['R_of_boundary_points[m]'],R_of_boundary_points))
    file_dict[itime].update(dict.fromkeys(['Z_of_boundary_points[m]'],Z_of_boundary_points))
    file_dict[itime].update(dict.fromkeys(['R_of_cur_density[m]'],R_of_cur_density))
    file_dict[itime].update(dict.fromkeys(['current_density[A/m2]'],current_density))
  return file_dict

# ------------------------------------------------------------------------------------------------

def dictionary_from_profiles(filename):
  ascii_file = open(filename,'r')
  data_from_file = ascii_file.read()
  ascii_file.close()
  # PACK OF DATA FOR EACH TIME SLICE
  packdata = data_from_file.split('n,tt=')
  packdata = packdata[1:]
  ntime = len(packdata)
  file_dict = {}
  for itime in range(ntime):
    tslice_lines = packdata[itime].split('\n')
    firstline = tslice_lines[0].split()
    nradial = int(firstline[0])
    Time = float(firstline[1])*1.e-3
    values = np.asfarray(''.join(tslice_lines[1:]).split(),float)
    psi_normalized                   = values[0::10]       # 1
    electron_temperature             = values[1::10]       # EV
    electron_density                 = values[2::10]*1.e19 # M-3
    effective_charge                 = values[3::10]       # 1
    parallel_heat_flux               = values[4::10]*1.e9  # W/M2
    power_density_of_conductive_loss = values[5::10]*1.e12 # W/M3
    plasma_current_density           = values[6::10]*1.e7  # A/M2
    parallel_electric_field          = values[7::10]       # V/M
    radiated_power_density           = values[8::10]*1.e12 # W/M3
    runaway_current_density          = values[9::10]*1.e7  # A/M2
    file_dict[itime] = dict.fromkeys(['Time[s]'],Time)
    file_dict[itime].update(dict.fromkeys(['normalized_psi'],psi_normalized))
    file_dict[itime].update(dict.fromkeys(['Te[eV]'],electron_temperature))
    file_dict[itime].update(dict.fromkeys(['ne[m-3]'],electron_density))
    file_dict[itime].update(dict.fromkeys(['Zeff'],effective_charge))
    file_dict[itime].update(dict.fromkeys(['parallel_heat_flux[W/m2]'],parallel_heat_flux))
    file_dict[itime].update(dict.fromkeys(['power_density_of_conductive_loss[W/m3]'],power_density_of_conductive_loss))
    file_dict[itime].update(dict.fromkeys(['plasma_current_density[A/m2]'],plasma_current_density))
    file_dict[itime].update(dict.fromkeys(['parallel_electric_field[V/m]'],parallel_electric_field))
    file_dict[itime].update(dict.fromkeys(['radiated_power_density[W/m3]'],radiated_power_density))
    file_dict[itime].update(dict.fromkeys(['runaway_current_density[A/m2]'],runaway_current_density))
  return file_dict


# ------------------------------------------------------------------------------------------------

def dictionary_from_halo(filename):
  ascii_file = open(filename,'r')
  data_from_file = ascii_file.read()
  ascii_file.close()
  # SEPARATED DATA FOR EACH TIME SLICE (1ST = DESCRIPTION, REST = ALL TIME SLICES)
  packdata = data_from_file.split('Time [ms]')
  # DESCRIPTION = FIRST PACK
  description = packdata[0]
  version = str([line for line in description.split('\n') if "version" in line])
  # DATA CONTENT = ALL OTHER PACKS
  data = packdata[1:]
  # CONSTRUCTION OF THE DICTIONARY
  npoints   = 10
  nvar      = 5
  file_dict = {}
  for itime in range(len(data)):
    tslice_lines = data[itime].split('\n')
    Time    = float(tslice_lines[1])*1.e-3
    values  = tslice_lines[3].split()
    Rtch    = values[0]
    Ztch    = values[1]
    values  = (str.replace((','.join(np.asarray(tslice_lines[5:]))),',','')).split()
    Rcw   = [None]*npoints
    Zcw   = [None]*npoints
    Rccw  = [None]*npoints
    Zccw  = [None]*npoints
    Ihpol = [None]*npoints
    for ipoint in range(npoints):
      Rcw[ipoint]   = float(values[0+ipoint*nvar])
      Zcw[ipoint]   = float(values[1+ipoint*nvar])
      Rccw[ipoint]  = float(values[2+ipoint*nvar])
      Zccw[ipoint]  = float(values[3+ipoint*nvar])
      Ihpol[ipoint] = float(values[4+ipoint*nvar])
    file_dict[itime] = dict.fromkeys(['Time[s]'],Time)
    file_dict[itime].update(dict.fromkeys(['Rtch[cm]'],Rtch))
    file_dict[itime].update(dict.fromkeys(['Ztch[cm]'],Ztch))
    file_dict[itime].update(dict.fromkeys(['Rcw[cm]'],Rcw))
    file_dict[itime].update(dict.fromkeys(['Zcw[cm]'],Zcw))
    file_dict[itime].update(dict.fromkeys(['Rccw[cm]'],Rccw))
    file_dict[itime].update(dict.fromkeys(['Zccw[cm]'],Zccw))
    file_dict[itime].update(dict.fromkeys(['Ihpol[kA]'],Ihpol))    
  return file_dict,version

# ------------------------------------------------------------------------------------------------

def dictionary_from_for002(filename):
  ascii_file = open(filename,'r')
  data_from_file = ascii_file.read()
  ascii_file.close()

  data_lines = data_from_file.split('\n')

  file_dict = {}
  file_dict['nrad']       = int(data_lines[1].split()[0])    # number of transport grid points
  file_dict['mplasma']    = int(data_lines[1].split()[1])    # number of divisions in poloidal direction
  file_dict['next']       = int(data_lines[1].split()[2])    # number of initial time steps when simulations are carried out with fixed plasma current value
  file_dict['tt']         = float(data_lines[3].split()[0])  # initial simulation time moment, ms
  file_dict['tay']        = float(data_lines[3].split()[1])  # initial time step, ms .... (or the time step before thermal quench, ms ?)
  file_dict['t_end']      = float(data_lines[3].split()[2])  # final simulation time moment, ms
  file_dict['RS0']        = float(data_lines[3].split()[3])  # major radius where vacuum toroidal field bt0 (see below) is defined, cm
  file_dict['psend']      = float(data_lines[3].split()[4])  # special parameter

  file_dict['i_graph']    = int(data_lines[5].split()[0])    # key to switch on (=1) or off (=0) the plasma equilibrium interactive plotting during simulations
  file_dict['ALFA0']      = float(data_lines[7].split()[0])  # the values to be used to specify the plasma current profile j(rho) for initial equilibrium: initial normalized coefficient for initial plasma current profile
  file_dict['BETA']       = float(data_lines[7].split()[1])  # the values to be used to specify the plasma current profile j(rho) for initial equilibrium.
  file_dict['alfa1']      = float(data_lines[7].split()[2])  # 
  file_dict['omega']      = float(data_lines[7].split()[3])  # averaging in time plasma current profile coefficient

  file_dict['iread']      = int(data_lines[9].split()[0])    # ???
  file_dict['kzero']      = int(data_lines[9].split()[1])    # ???
  file_dict['IWRITE']     = int(data_lines[9].split()[2])    # ???
  file_dict['kEFIT']      = int(data_lines[9].split()[3])    # ???

#  file_dict['alfax1']     = float(data_lines[11].split()[0]) # not used
#  file_dict['alfax2']     = float(data_lines[11].split()[1]) # not used
#  file_dict['betax1']     = float(data_lines[11].split()[2]) # not used
#  file_dict['betax2']     = float(data_lines[11].split()[3]) # not used

#  file_dict['pw_1']       = float(data_lines[13].split()[0]) # not used
#  file_dict['pw_2']       = float(data_lines[13].split()[1]) # not used

  file_dict['te_a']       = float(data_lines[15].split()[0]) # initial electron and ion temperatures in plasma center, eV (if kcchp = 0  see below description)
#  file_dict['ti_a']       = float(data_lines[15].split()[1]) # not used
  file_dict['te_b']       = float(data_lines[15].split()[2]) # electron and ion temperatures in plasma boundary, eV
  file_dict['ti_b']       = float(data_lines[15].split()[3]) # ???
  file_dict['pw_e']       = float(data_lines[15].split()[4]) # power for initial temperature profile

  file_dict['pd0_a']      = float(data_lines[17].split()[0]) # density of D and T ions in plasma center, 10^13 cm-3
  file_dict['pt0_a']      = float(data_lines[17].split()[1]) # density of D and T ions in plasma center, 10^13 cm-3
  file_dict['pd0_b']      = float(data_lines[17].split()[2]) # density of D and T ions in plasma boundary, 10^13 cm-3
  file_dict['pt0_b']      = float(data_lines[17].split()[3]) # density of D and T ions in plasma boundary, 10^13 cm-3
  file_dict['pw_p']       = float(data_lines[17].split()[4]) # power for initial density profile

#  file_dict['zeff_a']     = float(data_lines[19].split()[0]) # not used
#  file_dict['zeff_b']     = float(data_lines[19].split()[1]) # not used

  file_dict['SIG0']       = float(data_lines[21].split()[0]) # unit coefficient for plasma resistance

  file_dict['zhib']       = float(data_lines[23].split()[0]) # normalized coefficient for thermo conductivity
  file_dict['tego']       = float(data_lines[23].split()[1]) # normalized coefficient for thermo conductivity
  file_dict['zalfa']      = float(data_lines[23].split()[2]) # Z value of alpha particle
#  file_dict['talfa']      = float(data_lines[23].split()[3]) # not used
#  file_dict['alp1']       = float(data_lines[23].split()[4]) # not used

  file_dict['ktp']        = int(data_lines[25].split()[0])   # normalized coefficients for particle pinch velocity
  file_dict['kpin']       = int(data_lines[25].split()[1])   # normalized coefficients for particle pinch velocity
  file_dict['ken']        = int(data_lines[25].split()[2])   # normalized coefficients for particle pinch velocity
  file_dict['ken1']       = int(data_lines[25].split()[3])   # normalized coefficients for particle pinch velocity
  file_dict['ken2']       = int(data_lines[25].split()[4])   # normalized coefficients for particle pinch velocity
  file_dict['kd2']        = int(data_lines[25].split()[5])   # normalized coefficients for particle pinch velocity
  file_dict['nal']        = int(data_lines[25].split()[6])   # normalized coefficient for alpha particle source

#  file_dict['alpy']       = float(data_lines[27].split()[0]) # not used
#  file_dict['ppp']        = float(data_lines[27].split()[1]) # not used
#  file_dict['eee']        = float(data_lines[27].split()[2]) # not used
#  file_dict['dd']         = float(data_lines[27].split()[3]) # not used
#  file_dict['dt']         = float(data_lines[27].split()[4]) # not used
#  file_dict['dh']         = float(data_lines[27].split()[5]) # not used
#  file_dict['df']         = float(data_lines[27].split()[6]) # not used

  file_dict['lt']         = int(data_lines[29].split()[0])   # indexes of boundary conditions for diffusion equations
  file_dict['ld']         = int(data_lines[29].split()[1])   # indexes of boundary conditions for diffusion equations
  file_dict['lh']         = int(data_lines[29].split()[2])   # indexes of boundary conditions for diffusion equations
  file_dict['ll']         = int(data_lines[29].split()[3])   # indexes of boundary conditions for diffusion equations
  file_dict['lm']         = int(data_lines[29].split()[4])   # indexes of boundary conditions for diffusion equations
  file_dict['it']         = int(data_lines[29].split()[5])   # indexes of boundary conditions for diffusion equations
  file_dict['id']         = int(data_lines[29].split()[6])   # indexes of boundary conditions for diffusion equations
  file_dict['ih']         = int(data_lines[29].split()[7])   # indexes of boundary conditions for diffusion equations

  file_dict['eps0']       = float(data_lines[31].split()[0]) # simulation accuracies
  file_dict['eps1']       = float(data_lines[31].split()[1]) # simulation accuracies
  file_dict['eps2']       = float(data_lines[31].split()[2]) # simulation accuracies

  file_dict['anoma_e']    = float(data_lines[33].split()[0]) # anomalous electron thermo conductivity coefficient
  file_dict['anoma_i']    = float(data_lines[33].split()[1]) # anomalous ion thermo conductivity coefficient
  file_dict['key_t11']    = float(data_lines[33].split()[2]) # key to switch on T-11 Ohmic scaling
  file_dict['kcchp']      = float(data_lines[33].split()[3]) # key to prescribe either plasma temperature or plasma density in the beginning (“0”  temperature, “1”  density, in this case an average plasma density value in 10^13 cm-3 is taken from the file dens.dat)

  file_dict['edope']      = float(data_lines[35].split()[0]) # normalized coefficient for additional heating
  file_dict['edopi']      = float(data_lines[35].split()[1]) # normalized coefficient for additional heating

  file_dict['udd']        = float(data_lines[37].split()[0]) # additional voltage in plasma boundary (for testing)

#  file_dict['k_ener']     = int(data_lines[39].split()[0])   # not used
#  file_dict['k_uv']       = int(data_lines[39].split()[1])   # not used

  file_dict['t_dop']      = float(data_lines[41].split()[0]) # time to switch on additional heating

  file_dict['r0']         = float(data_lines[43].split()[0]) # coordinates for plasma center for DINA initialization
  file_dict['z0']         = float(data_lines[43].split()[1]) # coordinates for plasma center for DINA initialization
#  file_dict['zref']       = float(data_lines[43].split()[2]) # not used

  file_dict['kzref']      = int(data_lines[45].split()[0])   # regularization keys for initial equilibrium
  file_dict['krref']      = int(data_lines[45].split()[1])   # regularization keys for initial equilibrium
  file_dict['key_b']      = int(data_lines[45].split()[2])   # key of boundary condition for Grad-Shafranov equation
#  file_dict['i_pf']       = int(data_lines[45].split()[3])   # not used

  file_dict['i_c']        = int(data_lines[47].split()[0])   # key to switch on/off moving grid solver

  file_dict['q_vde']      = float(data_lines[49].split()[0]) # the value of q95 when thermal quench during Hot VDE is happens (not used?)

  file_dict['tay_00']     = float(data_lines[51].split()[0]) # the time step after thermal quench, ms (not used?)
#  file_dict['tay_th']     = float(data_lines[51].split()[1]) # not used
  file_dict['t_disr']     = float(data_lines[51].split()[2]) # decay time of plasma temperature during thermal quench, ms (not used?)

#  file_dict['d_tpl']      = float(data_lines[53].split()[0]) # not used
#  file_dict['tpl_end']    = float(data_lines[53].split()[1]) # not used

#  file_dict['c_h']        = float(data_lines[55].split()[0]) # not used
#  file_dict['d_halo']     = float(data_lines[55].split()[1]) # not used

#  file_dict['kmaj']       = int(data_lines[57].split()[0])   # not used
#  file_dict['li_drop']    = int(data_lines[57].split()[1])   # not used
  file_dict['ndisrup']    = int(data_lines[57].split()[2])   # the number of steps, after which the thermal quench will stall if q_vde is not used
#  file_dict['n_dif']      = int(data_lines[57].split()[3])   # not used
  file_dict['nmix']       = int(data_lines[57].split()[4])   # the number of points along the radius for the plasma current profile flattening

#  file_dict['hpart']      = float(data_lines[59].split()[0]) # not used
#  file_dict['te_h']       = float(data_lines[59].split()[1]) # not used

#  file_dict['i_d3d']      = int(data_lines[61].split()[0])   # not used
  file_dict['i_iter']     = int(data_lines[61].split()[1])   # key to switch on the old L-mode scaling
  file_dict['i_smal']     = int(data_lines[61].split()[2])   # key to switch on Lacner-Gottardi scaling

#  file_dict['ngra']       = int(data_lines[63].split()[0])   # not used
  file_dict['i_ramp']     = int(data_lines[63].split()[1])   # special regularization keys
  file_dict['i_v']        = int(data_lines[63].split()[2])   # special regularization keys
  file_dict['i_con']      = int(data_lines[63].split()[3])   # special regularization keys

  file_dict['tpl']        = float(data_lines[65].split()[0]) # plasma current value, kA
  file_dict['bt0']        = float(data_lines[65].split()[1]) # vacuum toroidal magnetic field in rs0 (see above) major radius position, kGs
  file_dict['eu']         = float(data_lines[65].split()[2]) # plasma minor radius for DINA initialization, cm
  file_dict['elong']      = float(data_lines[65].split()[3]) # plasma elongation for DINA initialization

  file_dict['e_sep']      = float(data_lines[67].split()[0]) # special accuracy value for Grad-Shafranov equation solver

#  file_dict['i_beta']     = int(data_lines[69].split()[0])   # not used
  file_dict['i_gap5']     = int(data_lines[69].split()[1])   # key to switch on/off thermal source in electron and ion energy

#  file_dict['i_br']       = int(data_lines[71].split()[0])   # not used

  file_dict['ind_r1']     = int(data_lines[73].split()[0])   # keys to control initial equilibrium in R direction
  file_dict['ind_r2']     = int(data_lines[73].split()[1])   # keys to control initial equilibrium in R direction
  file_dict['ind_z1']     = int(data_lines[73].split()[2])   # keys to control initial equilibrium in Z direction
  file_dict['ind_z2']     = int(data_lines[73].split()[3])   # keys to control initial equilibrium in Z direction

  file_dict['i_feed']     = int(data_lines[75].split()[0])   # keys to control initial equilibrium
  file_dict['i_ecoil']    = int(data_lines[75].split()[1])   # keys to control initial equilibrium

#  file_dict['te_a0']      = float(data_lines[77].split()[0]) # not used

  file_dict['alf_dis']    = float(data_lines[79].split()[0]) # key to control an expansion of halo area

  file_dict['alf_pas1']   = int(data_lines[81].split()[0])   # ???
  file_dict['alf_pas2']   = int(data_lines[81].split()[1])   # ???
  file_dict['alf_pas3']   = int(data_lines[81].split()[2])   # ???

  file_dict['i_qpar']     = int(data_lines[83].split()[0])   # ???

  file_dict['key_run']    = int(data_lines[85].split()[0])   # ???

  return file_dict

# ------------------------------------------------------------------------------------------------

def dictionary_from_toroidal_currents(filename):
  ascii_file = open(filename,'r')
  data_from_file = ascii_file.read()
  ascii_file.close()
  # SEPARATED DATA FOR EACH TIME SLICE (1ST = NOT TIME-DEPENDENT, REST = TIME EVOLVING VARIABLES)
  packdata = data_from_file.split('Time [ms]')
  # FIRST_DATA_SET = FIRST PACK (VARIABLES NOT TIME-DEPENDENT)
  first_data_set = packdata[0]
  first_data_set_lines = first_data_set.split('\n')
  R_coordinate_strings  = first_data_set_lines[14].split() + first_data_set_lines[15].split()
  Z_coordinate_strings  = first_data_set_lines[18].split() + first_data_set_lines[19].split()
  radial_size_strings   = first_data_set_lines[22].split() + first_data_set_lines[23].split()
  vertical_size_strings = first_data_set_lines[26].split() + first_data_set_lines[27].split()
  name_of_coils = ['CS3U','CS2U','CS1U','CS1L','CS2L','CS3L','PF1','PF2','PF3','PF4','PF5','PF6']
  for icoil in range(len(name_of_coils)):
    if icoil==0:
      dict_R_coordinate  = dict([(name_of_coils[icoil],float(R_coordinate_strings[icoil]))])
      dict_Z_coordinate  = dict([(name_of_coils[icoil],float(Z_coordinate_strings[icoil]))])
      dict_radial_size   = dict([(name_of_coils[icoil],float(radial_size_strings[icoil]))])
      dict_vertical_size = dict([(name_of_coils[icoil],float(vertical_size_strings[icoil]))])
    else:
      dict_R_coordinate.update(dict([(name_of_coils[icoil],float(R_coordinate_strings[icoil]))]))
      dict_Z_coordinate.update(dict([(name_of_coils[icoil],float(Z_coordinate_strings[icoil]))]))
      dict_radial_size.update(dict([(name_of_coils[icoil],float(radial_size_strings[icoil]))]))
      dict_vertical_size.update(dict([(name_of_coils[icoil],float(vertical_size_strings[icoil]))]))
  nfilaments = int(first_data_set_lines[32])
  filament_data = first_data_set_lines[36:]
  nlines = len(filament_data)
  ifil = 0
  r1 = []
  z1 = []
  r2 = []
  z2 = []
  turn = []
  nturn = []
  for iline in range(nlines):
    if filament_data[iline].strip()=='1':
      fil_data = filament_data[iline+1].split()
      r1.append(float(fil_data[0]))
      z1.append(float(fil_data[1]))
      r2.append(float(fil_data[2]))
      z2.append(float(fil_data[3]))
      turn.append(float(fil_data[4]))
      nturn.append(1)
    if filament_data[iline].strip()=='2':
      fil_data = filament_data[iline+1].split()+filament_data[iline+2].split()
      r1.append(float(fil_data[0]))
      z1.append(float(fil_data[1]))
      r2.append(float(fil_data[2]))
      z2.append(float(fil_data[3]))
      turn.append(float(fil_data[4]))
      nturn.append(2)
      r1.append(float(fil_data[5]))
      z1.append(float(fil_data[6]))
      r2.append(float(fil_data[7]))
      z2.append(float(fil_data[8]))
      turn.append(float(fil_data[9]))
      nturn.append(-2) # Trick to differentiate the second turn of \
                       # the same element when I read it in dina_data_2_imas
  # SECOND_DATA_SET = SECOND PACK (TIME EVOLVING VARIABLES)
  second_data_set = packdata[1:]
  ntime=len(second_data_set)
  dict_time = {}
  for itime in range(ntime):
    tslice_lines = second_data_set[itime].split('\n')
    Time = float(tslice_lines[1])*1.e-3
    values = ''.join(tslice_lines[3:-4]).split('Number of plasma current filaments')
    values_passive = values[0].split()
    current_filament_passive = [None]*len(values_passive)
    for icur in range(len(values_passive)):
      current_filament_passive[icur] = float(values_passive[icur])
    values_plasma  = values[1].split()
    n_plasma_filaments = int(values_plasma[0])
    values_filaments = values_plasma[9:]
    plasma_current_filaments_r = [None]*n_plasma_filaments
    plasma_current_filaments_z = [None]*n_plasma_filaments
    plasma_current_filaments_i = [None]*n_plasma_filaments
    for ifil in range(n_plasma_filaments):
      plasma_current_filaments_r[ifil] = values_filaments[ifil*3]
      plasma_current_filaments_z[ifil] = values_filaments[ifil*3+1]
      plasma_current_filaments_i[ifil] = values_filaments[ifil*3+2]
    full_coil_currents = ''.join(tslice_lines[-4:]).split()
    values_coil_currents = full_coil_currents
    if itime == ntime-1:
      values_coil_currents = full_coil_currents[5:]
    for icoil in range(len(name_of_coils)):
      if icoil==0:
        dict_coils_currents = dict([(name_of_coils[icoil],float(values_coil_currents[icoil]))])
      else:
        dict_coils_currents.update(dict([(name_of_coils[icoil],float(values_coil_currents[icoil]))]))
    if itime==0:
      dict_time.update([('Time',[Time])])
      dict_time.update([('current_filament_passive[kA]',[current_filament_passive])])
      dict_time.update([('plasma_current_filaments_r[cm]',[plasma_current_filaments_r])])
      dict_time.update([('plasma_current_filaments_z[cm]',[plasma_current_filaments_z])])
      dict_time.update([('plasma_current_filaments_i[kA]',[plasma_current_filaments_i])])
      dict_time.update([('coils_currents[kA*turns]',[dict_coils_currents])])
    else:
      dict_time['Time'].append(Time)
      dict_time['current_filament_passive[kA]'].append(current_filament_passive)
      dict_time['plasma_current_filaments_r[cm]'].append(plasma_current_filaments_r)
      dict_time['plasma_current_filaments_z[cm]'].append(plasma_current_filaments_z)
      dict_time['plasma_current_filaments_i[kA]'].append(plasma_current_filaments_i)
      dict_time['coils_currents[kA*turns]'].append(dict_coils_currents)
  # FIRST DATA PACK; NOT TIME-DEPENDENT
  file_dict = dict([('coils_R_coordinate[cm]',dict_R_coordinate)])
  file_dict.update(dict([('coils_Z_coordinate[cm]',dict_Z_coordinate)]))
  file_dict.update(dict([('coils_radial_size[cm]',dict_radial_size)]))
  file_dict.update(dict([('coils_vertical_size[cm]',dict_vertical_size)]))
  file_dict.update([('r1[cm]',r1)])
  file_dict.update([('r2[cm]',r2)])
  file_dict.update([('z1[cm]',z1)])
  file_dict.update([('z2[cm]',z2)])
  file_dict.update([('turn',turn)])
  file_dict.update([('nturn',nturn)])
  # SECOND DATA PACK: TIME-DEPENDENT
  file_dict.update([('dict_time',dict_time)])
  return file_dict

# ----
# To list the keys in dict_time sub-dictionary:
# file_dict.get('dict_time').keys()
# 
# To access coils_vertical_size:
# file_dict.get('coils_vertical_size[cm]')
#
# To access current_filament_passive values:
# file_dict.get('dict_time').get('current_filament_passive[kA]')
#
# To access coils_currents values for 1st time slice:
# file_dict.get('dict_time').get('coils_currents[kA*turns]')[0].get('CS3U')
# ----

