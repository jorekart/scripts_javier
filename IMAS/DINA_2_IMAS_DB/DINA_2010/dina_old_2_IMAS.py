# -----------------------------------
# PURPOSE:
# TO READ DINA INPUT DATA FROM FILES
# AND CONVERT THEM TO IMAS
# -----------------------------------
def dina_ds_2_imas(database,user_or_path,shot,run,reference_name):

  # OUTPUT DATABASE AND FILE
  #import os
  #database,user_or_path = 'dina_disruption',os.getenv('USER')
  #shot,run = 8877,1

  # NECESSARY MODULES
  import imas, os, datetime
  import numpy as np
  from read_dina_ds_files import dictionary_from_file, dictionary_from_dens, \
    dictionary_from_psi_data, dictionary_from_profiles, dictionary_from_halo, \
    dictionary_from_for002, dictionary_from_toroidal_currents

  # LOCAL DATABASE ENVIRONMENT (FOR IMAS STORAGE)
  output = imas.DBEntry(imas.imasdef.HDF5_BACKEND,database,shot,run,user_or_path)
  output.create()

  # WHAT IS THE CURRENT DATE AND TIME?
  date_and_time = datetime.datetime.now()
  now = date_and_time.strftime("%Y-%m-%d %H:%M")

  # -------------------------------------------------
  # CREATE DICTIONARIES FROM DINA OUTPUT ASCII FILES
  # -------------------------------------------------
  print('Read DINA ascii files.')
  print('--> Read plasma.dat')
  plasma_dict            = dictionary_from_file('plasma.dat')
  # print('--> Read powers.dat')
  # powers_dict            = dictionary_from_file('powers.dat')
  print('--> Read halo.dat')
  halo_dict,version      = dictionary_from_halo('halo.dat')
  print('--> Read toroidal_currents.dat')
  toroidal_currents_dict = dictionary_from_toroidal_currents('toroidal_currents.dat')
  # print('--> Read ksuprof')
  # profiles_dict          = dictionary_from_profiles('ksuprof')
  print('--> Read psi_data')
  psi_dict               = dictionary_from_psi_data('psi_data')
  # print('--> Read for002')
  # for002_dict            = dictionary_from_for002('for002')
  # print('--> Read dens.dat, dens_imp,dat, dens_imp2.dat')
  # dens1,dens2,dens3      = dictionary_from_dens('dens.dat','dens_imp.dat','dens_imp2.dat')
  # print('--> Read beforeTQ.eqdsk ')
  # try:
  #     eqdsk_dict         = read_eqdsk('beforeTQ.eqdsk')
  #     eqdsk_found        = True
  # except: 
  #     print('--> The DINA file beforeTQ.eqdsk was not found ')
  #     eqdsk_found        = False
  # print('Done reading DINA ascii files.')

  ntime_plasma            = len(plasma_dict)
  #ntime_powers            = len(powers_dict)
  #ntime_profiles          = len(profiles_dict)
  ntime_toroidal_currents = len(toroidal_currents_dict.get('dict_time').get('Time'))
  ntime_psi_data          = len(psi_dict)
  ntime_halo              = len(halo_dict)
  fact_Ip = -1.0  # Reverse sign of all currents?

  # TIME ARRAYS
  time_plasma  = [None]*ntime_plasma
  for itime in range(ntime_plasma):
    time_plasma[itime] = float(plasma_dict[itime].get('t[ms]'))*1.e-3

  # time_powers  = [None]*ntime_powers
  # for itime in range(ntime_powers):
  #   time_powers[itime] = float(powers_dict[itime].get('t'))*1.e-3

  time_psi_data  = [None]*ntime_psi_data
  for itime in range(ntime_psi_data):
    time_psi_data[itime] = float(psi_dict[itime].get('Time[s]'))

  time_halo  = [None]*ntime_halo
  for itime in range(ntime_halo):
    time_halo[itime] = float(halo_dict[itime].get('Time[s]'))

  # time_profiles  = [None]*ntime_profiles
  # for itime in range(ntime_profiles):
  #   time_profiles[itime] = float(profiles_dict[itime].get('Time[s]'))

  # ------------------------------------------------------------------------------------------------

  def find_nearest(a, a0):
      'Element in nd array `a` closest to the scalar value `a0`'
      idx = np.abs(a - a0).argmin()
      return a.flat[idx],idx

  # ------------------------------------------------------------------------------------------------

  # --------------------
  # FILL IDSS WITH DATA
  # --------------------
  summary        = imas.summary()
  # core_profiles  = imas.core_profiles()
  # core_transport = imas.core_transport()
  pf_active      = imas.pf_active()
  pf_passive     = imas.pf_passive()
  equilibrium    = imas.equilibrium()
  wall           = imas.wall()
  disruption     = imas.disruption()

  # SUMMARY IDS
  print('------> Fill summary IDS')
  summary.code.name                         = 'DINA ITER Disruption Simulator'
  summary.ids_properties.comment            = '2010 simulation'
  summary.ids_properties.homogeneous_time   = 1
  summary.ids_properties.source             = reference_name
  summary.ids_properties.provider           = 'Generated from DINA output raw data'
  summary.ids_properties.creation_date      = now

  summary.time.resize(ntime_plasma)
  summary.global_quantities.ip.value.resize(ntime_plasma)
  summary.global_quantities.beta_pol.value.resize(ntime_plasma)
  summary.global_quantities.volume.value.resize(ntime_plasma)
  summary.runaways.current.value.resize(ntime_plasma)
  summary.volume_average.n_e.value.resize(ntime_plasma)
  summary.volume_average.t_e.value.resize(ntime_plasma)
  summary.volume_average.zeff.value.resize(ntime_plasma)
  summary.boundary.minor_radius.value.resize(ntime_plasma)
  summary.local.magnetic_axis.q.value.resize(ntime_plasma)
  # summary.boundary.magnetic_axis_r.value.resize(ntime_plasma)
  # summary.boundary.magnetic_axis_z.value.resize(ntime_plasma)
  summary.boundary.elongation.value.resize(ntime_plasma)
  # summary.local.separatrix.q.value.resize(ntime_plasma)
  summary.global_quantities.q_95.value.resize(ntime_plasma)
  summary.global_quantities.li_3.value.resize(ntime_plasma)
  # summary.global_quantities.power_ohm.value.resize(ntime_plasma)
  # summary.line_average.n_i_total.value.resize(ntime_plasma)
  # summary.local.magnetic_axis.e_field_parallel.value.resize(ntime_plasma)
  # summary.local.separatrix.e_field_parallel.value.resize(ntime_plasma)

  for itime in range(ntime_plasma):
    summary.time[itime]                                       = float(plasma_dict[itime].get('t[ms]'))*1.e-3
    summary.global_quantities.ip.value[itime]                 = float(plasma_dict[itime].get('Ip[kA]'))*1.e3*fact_Ip
    summary.global_quantities.beta_pol.value[itime]           = float(plasma_dict[itime].get('BETAp'))
    summary.global_quantities.volume.value[itime]             = float(plasma_dict[itime].get('V_sep[cm3]')) * 1e-6
    summary.runaways.current.value[itime]                     = float(plasma_dict[itime].get('Ire[kA]'))*1.e3*fact_Ip
    summary.volume_average.n_e.value[itime]                   = float(plasma_dict[itime].get('Ne[13,cm-3]'))*1.e19
    summary.volume_average.t_e.value[itime]                   = float(plasma_dict[itime].get('Te[eV]'))
    summary.volume_average.zeff.value[itime]                  = float(plasma_dict[itime].get('Zeff'))
    summary.boundary.minor_radius.value[itime]                = float(plasma_dict[itime].get('a[cm]'))*1e-2
    #summary.boundary.magnetic_axis_r.value[itime]             = float(plasma_dict[itime].get('rmag')) # replace with centroid
    #summary.boundary.magnetic_axis_z.value[itime]             = float(plasma_dict[itime].get('zmag')) # replace with centroid
    summary.boundary.elongation.value[itime]                  = float(plasma_dict[itime].get('k'))
    summary.local.magnetic_axis.q.value[itime]                = float(plasma_dict[itime].get('qa'))  # This is the value at magnetic axis
    summary.global_quantities.q_95.value[itime]               = float(plasma_dict[itime].get('q95'))
    summary.global_quantities.li_3.value[itime]                = float(plasma_dict[itime].get('Li'))
    #summary.line_average.n_i_total.value[itime]               = float(plasma_dict[itime].get('Nimp'))*1.e19
    #summary.local.magnetic_axis.e_field_parallel.value[itime] = float(plasma_dict[itime].get('E_ax'))
    #summary.local.separatrix.e_field_parallel.value[itime]    = float(plasma_dict[itime].get('E-bnd'))

  # for itime in range(ntime_powers):
  #   [tc,it] = find_nearest(np.asarray(time_plasma),time_powers[itime])
  #   summary.global_quantities.power_ohm.value[it] = float(powers_dict[itime].get('P_jh[MW]'))

  # CORE_PROFILES IDS
  # print('------> Fill core_profiles IDS')
  # core_profiles.code = summary.code
  # core_profiles.ids_properties = summary.ids_properties

  # core_profiles.time.resize(ntime_profiles)
  # core_profiles.profiles_1d.resize(ntime_profiles)

  # for itime in range(ntime_profiles):
  #   core_profiles.time[itime] = float(profiles_dict[itime].get('Time[s]'))
  #   psi_norm         = profiles_dict[itime].get('normalized_psi')
  #   # denormalize psi profile
  #   psi_axis         = psi_dict[itime].get('poloidal_flux_at_magnetic_axis[Wb]')
  #   psi_bnd_halo     = psi_dict[itime].get('poloidal_flux_at_halo_boundary[Wb]')
  #   psi_LCFS         = psi_dict[itime].get('poloidal_flux_at_plasma_boundary[Wb]')
  #   if (psi_bnd_halo == 0):   # Normalized to halo boundary, when halo is active
  #     psi_bnd = psi_LCFS
  #   else:
  #     psi_bnd = psi_bnd_halo
  #   psi = psi_norm*(psi_bnd-psi_axis) + psi_axis     
  #   core_profiles.profiles_1d[itime].grid.psi = psi
  #   core_profiles.profiles_1d[itime].electrons.temperature = profiles_dict[itime].get('Te[eV]')
  #   core_profiles.profiles_1d[itime].electrons.density = profiles_dict[itime].get('ne[m-3]')
  #   core_profiles.profiles_1d[itime].e_field.parallel = profiles_dict[itime].get('parallel_electric_field[V/m]')
  #   core_profiles.profiles_1d[itime].j_total = profiles_dict[itime].get('plasma_current_density[A/m2]')
  #   core_profiles.profiles_1d[itime].zeff = profiles_dict[itime].get('Zeff')

  # # CORE_TRANSPORT IDS
  # print('------> Fill core_transport IDS')
  # core_transport.code = summary.code
  # core_transport.ids_properties = summary.ids_properties

  # core_transport.time.resize(ntime_profiles)
  # core_transport.model.resize(1)
  # core_transport.model[0].profiles_1d.resize(ntime_profiles)

  # for itime in range(ntime_profiles):
  #   core_transport.time[itime] = float(profiles_dict[itime].get('Time[s]'))
  #   core_transport.model[0].profiles_1d[itime].electrons.energy.flux = profiles_dict[itime].get('parallel_heat_flux[W/m2]')

  # PF_ACTIVE IDS
  print('------> Fill pf_active IDS')
  pf_active.code = summary.code
  pf_active.ids_properties = summary.ids_properties
  pf_active.time.resize(ntime_toroidal_currents)

  ncoils = len(toroidal_currents_dict.get('coils_R_coordinate[cm]').keys())
  pf_active.coil.resize(ncoils)

  for icoil in range(ncoils):
    pf_active.coil[icoil].element.resize(1)

  for icoil in range(ncoils):
    pf_active.coil[icoil].name = list(toroidal_currents_dict.get('coils_R_coordinate[cm]'))[icoil]
    pf_active.coil[icoil].element[0].turns_with_sign = 1.0
    pf_active.coil[icoil].element[0].geometry.rectangle.r = list(toroidal_currents_dict.get('coils_R_coordinate[cm]').values())[icoil]*1.e-2
    pf_active.coil[icoil].element[0].geometry.rectangle.z = list(toroidal_currents_dict.get('coils_Z_coordinate[cm]').values())[icoil]*1.e-2
    pf_active.coil[icoil].element[0].geometry.rectangle.width = list(toroidal_currents_dict.get('coils_radial_size[cm]').values())[icoil]*1.e-2
    pf_active.coil[icoil].element[0].geometry.rectangle.height = list(toroidal_currents_dict.get('coils_vertical_size[cm]').values())[icoil]*1.e-2

  pf_active.time.resize(ntime_toroidal_currents)

  for icoil in range(ncoils):
    pf_active.coil[icoil].current.data.resize(ntime_toroidal_currents)

  for itime in range(ntime_toroidal_currents):
    pf_active.time[itime] = toroidal_currents_dict.get('dict_time').get('Time')[itime]
    for icoil in range(ncoils):
      pf_active.coil[icoil].current.data[itime] = \
      toroidal_currents_dict.get('dict_time').get('coils_currents[kA*turns]')[itime].get(pf_active.coil[icoil].name)*1.e3*fact_Ip

  # PF_PASSIVE IDS
  print('------> Fill pf_passive IDS')
  pf_passive.code = summary.code
  pf_passive.ids_properties = summary.ids_properties

  nloop=0
  for ifilament in range(len(toroidal_currents_dict.get('r1[cm]'))):
    if toroidal_currents_dict.get('nturn')[ifilament]==1:
      nloop = nloop + 1
    else:
      if toroidal_currents_dict.get('nturn')[ifilament]==2:
        nloop = nloop + 1
  pf_passive.loop.resize(nloop)

  # nloop = total number of elements (one vv elements has 1 turn, one bm element has 2 turns)
  # (trick: i have put nturn=-2 for the second turn of bms to evaluate nloop more easily and to better distinguish 1st/2nd turns
  nloop = 0
  for ifilament in range(len(toroidal_currents_dict.get('r1[cm]'))):
    if toroidal_currents_dict.get('nturn')[ifilament]==1:
      nloop = nloop + 1
      pf_passive.loop[nloop-1].element.resize(1)
      pf_passive.loop[nloop-1].element.resize(1)
      pf_passive.loop[nloop-1].element[0].geometry.outline.r.resize(2)
      pf_passive.loop[nloop-1].element[0].geometry.outline.z.resize(2)
    else:
      if toroidal_currents_dict.get('nturn')[ifilament]==2:
        nloop = nloop + 1
        pf_passive.loop[nloop-1].element.resize(2)
        pf_passive.loop[nloop-1].element[0].geometry.outline.r.resize(2)
        pf_passive.loop[nloop-1].element[0].geometry.outline.z.resize(2)
        pf_passive.loop[nloop-1].element[1].geometry.outline.r.resize(2)
        pf_passive.loop[nloop-1].element[1].geometry.outline.z.resize(2)

  iloop = 0
  for ifilament in range(len(toroidal_currents_dict.get('r1[cm]'))):
    if toroidal_currents_dict.get('nturn')[ifilament]==1:
      iloop = iloop + 1
      pf_passive.loop[iloop-1].name = 'Vacuum Vessel'
      pf_passive.loop[iloop-1].element[0].area = 0. # TO MAKE IT CLEAR THAT IT IS A SEGMENT AND NOT AN ACTUAL OBJECT
      pf_passive.loop[iloop-1].element[0].geometry.outline.r = \
      np.array([toroidal_currents_dict.get('r1[cm]')[ifilament]*1.e-2,toroidal_currents_dict.get('r2[cm]')[ifilament]*1.e-2])
      pf_passive.loop[iloop-1].element[0].geometry.outline.z = \
      np.array([toroidal_currents_dict.get('z1[cm]')[ifilament]*1.e-2,toroidal_currents_dict.get('z2[cm]')[ifilament]*1.e-2])
      pf_passive.loop[iloop-1].element[0].turns_with_sign = toroidal_currents_dict.get('turn')[ifilament]
    else:
      if toroidal_currents_dict.get('nturn')[ifilament]==2:
        iloop = iloop + 1
        pf_passive.loop[iloop-1].name = 'Blanket Module'
        pf_passive.loop[iloop-1].element[0].area = 0. # TO MAKE IT CLEAR THAT IT IS A SEGMENT AND NOT AN ACTUAL OBJECT
        pf_passive.loop[iloop-1].element[0].geometry.outline.r = \
        np.array([toroidal_currents_dict.get('r1[cm]')[ifilament]*1.e-2,toroidal_currents_dict.get('r2[cm]')[ifilament]*1.e-2])
        pf_passive.loop[iloop-1].element[0].geometry.outline.z = \
        np.array([toroidal_currents_dict.get('z1[cm]')[ifilament]*1.e-2,toroidal_currents_dict.get('z2[cm]')[ifilament]*1.e-2])
        pf_passive.loop[iloop-1].element[0].turns_with_sign = toroidal_currents_dict.get('turn')[ifilament]
        pf_passive.loop[iloop-1].element[1].area = 0. # TO MAKE IT CLEAR THAT IT IS A SEGMENT AND NOT AN ACTUAL OBJECT
        pf_passive.loop[iloop-1].element[1].geometry.outline.r = \
        np.array([toroidal_currents_dict.get('r1[cm]')[ifilament+1]*1.e-2,toroidal_currents_dict.get('r2[cm]')[ifilament+1]*1.e-2])
        pf_passive.loop[iloop-1].element[1].geometry.outline.z = \
        np.array([toroidal_currents_dict.get('z1[cm]')[ifilament+1]*1.e-2,toroidal_currents_dict.get('z2[cm]')[ifilament+1]*1.e-2])
        pf_passive.loop[iloop-1].element[1].turns_with_sign = toroidal_currents_dict.get('turn')[ifilament+1]

  pf_passive.time.resize(ntime_toroidal_currents)

  for iloop in range(nloop):
    pf_passive.loop[iloop].current.resize(ntime_toroidal_currents)

  for itime in range(ntime_toroidal_currents):
    pf_passive.time[itime] = toroidal_currents_dict.get('dict_time').get('Time')[itime]
    iloop = 0
    n_fila = len(toroidal_currents_dict.get('r1[cm]'))
    for ifilament in range(n_fila):
      if toroidal_currents_dict.get('nturn')[ifilament]==1 or toroidal_currents_dict.get('nturn')[ifilament]==2:
        iloop = iloop + 1
        pf_passive.loop[iloop-1].current[itime] = toroidal_currents_dict.get('dict_time').get('current_filament_passive[kA]')[itime][iloop-1]*1.e3*fact_Ip

  # EQUILIBRIUM IDS
  print('------> Fill equilibrium IDS')
  equilibrium.code = summary.code
  equilibrium.ids_properties = summary.ids_properties
  

  # if eqdsk_found: 
  #     equilibrium.vacuum_toroidal_field.r0 = float(eqdsk_dict['rcentr'])
  #     equilibrium.vacuum_toroidal_field.b0.resize(ntime_toroidal_currents)

  equilibrium.grids_ggd.resize(ntime_toroidal_currents)

  equilibrium.time.resize(ntime_toroidal_currents)
  equilibrium.vacuum_toroidal_field.b0.resize(ntime_toroidal_currents)
  equilibrium.time_slice.resize(ntime_toroidal_currents)
  
  equilibrium.vacuum_toroidal_field.r0 = 6.2
  equilibrium.vacuum_toroidal_field.b0[:] = -5.3
  
  for itime in range(ntime_toroidal_currents):

    # if eqdsk_found: 
    #     equilibrium.vacuum_toroidal_field.b0[itime] = float(eqdsk_dict['bcentr'])
    equilibrium.grids_ggd[itime].grid.resize(1)
    equilibrium.grids_ggd[itime].grid[0].identifier.index = -1
    equilibrium.grids_ggd[itime].grid[0].identifier.description = \
      '2D current distribution as toroidal filaments, each filament has the ggd.j_phi current in A'

    equilibrium.time[itime] = toroidal_currents_dict.get('dict_time').get('Time')[itime]
    equilibrium.grids_ggd[itime].time = equilibrium.time[itime]
    nplasma_filaments = len(toroidal_currents_dict.get('dict_time').get('plasma_current_filaments_r[cm]')[itime])
    equilibrium.time_slice[itime].ggd.resize(1)
    equilibrium.time_slice[itime].ggd[0].r.resize(1)
    equilibrium.time_slice[itime].ggd[0].z.resize(1)
    equilibrium.time_slice[itime].ggd[0].j_phi.resize(1)
    equilibrium.time_slice[itime].ggd[0].r[0].values.resize(nplasma_filaments)
    equilibrium.time_slice[itime].ggd[0].z[0].values.resize(nplasma_filaments)
    equilibrium.time_slice[itime].ggd[0].j_phi[0].values.resize(nplasma_filaments)


  for itime in range(ntime_toroidal_currents):
    nplasma_filaments = len(toroidal_currents_dict.get('dict_time').get('plasma_current_filaments_r[cm]')[itime])
    for iplasma_filament in range(nplasma_filaments):
      equilibrium.time_slice[itime].ggd[0].r[0].values[iplasma_filament] = \
      float(toroidal_currents_dict.get('dict_time').get('plasma_current_filaments_r[cm]')[itime][iplasma_filament])*1.e-2
      equilibrium.time_slice[itime].ggd[0].z[0].values[iplasma_filament] = \
      float(toroidal_currents_dict.get('dict_time').get('plasma_current_filaments_z[cm]')[itime][iplasma_filament])*1.e-2
      equilibrium.time_slice[itime].ggd[0].j_phi[0].values[iplasma_filament] = \
      float(toroidal_currents_dict.get('dict_time').get('plasma_current_filaments_i[kA]')[itime][iplasma_filament])*1.e3*fact_Ip

  for itime in range(ntime_psi_data):
    if itime<len(equilibrium.time):
      if psi_dict[itime].get('Time[s]') == equilibrium.time[itime]:
        equilibrium.time_slice[itime].global_quantities.psi_axis     = psi_dict[itime].get('poloidal_flux_at_magnetic_axis[Wb]')
        equilibrium.time_slice[itime].global_quantities.psi_boundary = psi_dict[itime].get('poloidal_flux_at_plasma_boundary[Wb]')
        equilibrium.time_slice[itime].global_quantities.magnetic_axis.r = psi_dict[itime].get('R0_magnetic[m]')
        equilibrium.time_slice[itime].global_quantities.magnetic_axis.z = psi_dict[itime].get('Z0_magnetic[m]')
        equilibrium.time_slice[itime].profiles_2d.resize(1)
        equilibrium.time_slice[itime].profiles_2d[0].grid_type.name = 'RZ grid'
        equilibrium.time_slice[itime].profiles_2d[0].grid_type.index = 1
        equilibrium.time_slice[itime].profiles_2d[0].grid_type.description = 'Rectangular RZ grid'
        psi_map = np.reshape(psi_dict[itime].get('poloidal_flux_array[Wb]'),[psi_dict[itime].get('nz'),psi_dict[itime].get('nr')])
        equilibrium.time_slice[itime].profiles_2d[0].psi = np.transpose(psi_map)  # To be compatible with equiplot and other IMAS tools
        equilibrium.time_slice[itime].profiles_2d[0].grid.dim1 = psi_dict[itime].get('R_of_the_flux[m]')
        equilibrium.time_slice[itime].profiles_2d[0].grid.dim2 = psi_dict[itime].get('Z_of_the_flux[m]')
        equilibrium.time_slice[itime].boundary.psi = psi_dict[itime].get('poloidal_flux_at_plasma_boundary[Wb]')
        equilibrium.time_slice[itime].boundary.outline.r = psi_dict[itime].get('R_of_boundary_points[m]')
        equilibrium.time_slice[itime].boundary.outline.z = psi_dict[itime].get('Z_of_boundary_points[m]')
        # equilibrium.time_slice[itime].profiles_1d.r_outboard = psi_dict[itime].get('R_of_cur_density[m]')
        # equilibrium.time_slice[itime].profiles_1d.j_phi = psi_dict[itime].get('current_density[A/m2]')

  # WALL IDS
  print('------> Fill wall IDS')
  wall.code = summary.code
  wall.ids_properties = summary.ids_properties
  wall.description_2d.resize(1)
  wall.time.resize(1)
  #wall.global_quantities.current_tor.resize(ntime_psi_data)
  itime = 0

  wall.time[itime] = psi_dict[itime].get('Time[s]')
  wall.description_2d[itime].limiter.unit.resize(1)
  wall.description_2d[itime].limiter.unit[0].outline.r = psi_dict[itime].get('R_surrounding_limiter_contour[m]')
  wall.description_2d[itime].limiter.unit[0].outline.z = psi_dict[itime].get('Z_surrounding_limiter_contour[m]')
  # for itime in range(ntime_plasma):
  #   [tc,it] = find_nearest(np.asarray(time_psi_data),time_plasma[itime])
  #   wall.global_quantities.current_tor[it] = float(plasma_dict[itime].get('Ivv')) * 1e6  # From MA to A

  # DISRUPTION IDS
  print('------> Fill disruption IDS')
  disruption.code           = summary.code
  disruption.ids_properties = summary.ids_properties
  disruption.time.resize(ntime_plasma)
  disruption.global_quantities.current_halo_pol.resize(ntime_plasma)
  disruption.global_quantities.current_halo_phi.resize(ntime_plasma)
  # disruption.global_quantities.power_ohm.resize(ntime_plasma)
  # disruption.global_quantities.power_ohm_halo.resize(ntime_plasma)
  # disruption.global_quantities.power_radiated_electrons_impurities.resize(ntime_plasma)
  # disruption.global_quantities.power_radiated_electrons_impurities_halo.resize(ntime_plasma)
  # disruption.global_quantities.energy_ohm.resize(ntime_plasma)
  # disruption.global_quantities.energy_ohm_halo.resize(ntime_plasma)
  # disruption.global_quantities.energy_radiated_electrons_impurities.resize(ntime_plasma)
  # disruption.global_quantities.energy_radiated_electrons_impurities_halo.resize(ntime_plasma)
  disruption.global_quantities.psi_halo_boundary.resize(ntime_plasma)
  disruption.halo_currents.resize(ntime_plasma)

  npoints = len(halo_dict[0]['Rcw[cm]'])
  for itime in range(ntime_plasma):
    disruption.halo_currents[itime].area.resize(npoints)
    disruption.time[itime]                                                         = float(plasma_dict[itime].get('t[ms]'))*1.e-3
    disruption.global_quantities.current_halo_pol[itime]                           = float(plasma_dict[itime].get('Ihpol[kA]'))*1.e3
    disruption.global_quantities.current_halo_phi[itime]                           = float(plasma_dict[itime].get('Ihtor[kA]'))*1.e3*fact_Ip

  #for itime in range(ntime_powers):
    #[tc,it] = find_nearest(np.asarray(time_plasma),time_powers[itime])
    # disruption.global_quantities.power_ohm[it]                                  = float(powers_dict[itime].get('P_jh[MW]'))*1.e6
    # disruption.global_quantities.power_ohm_halo[it]                             = float(powers_dict[itime].get('P_jh_halo[MW]'))*1.e6
    # disruption.global_quantities.power_radiated_electrons_impurities[it]        = float(powers_dict[itime].get('P_imp[MW]'))*1.e6
    # disruption.global_quantities.power_radiated_electrons_impurities_halo[it]   = float(powers_dict[itime].get('P_imp_halo[MW]'))*1.e6
    # disruption.global_quantities.energy_ohm[it]                                 = float(powers_dict[itime].get('E_jh[MJ]'))*1.e6
    # disruption.global_quantities.energy_ohm_halo[it]                            = float(powers_dict[itime].get('E_jh_halo[MJ]'))*1.e6
    # disruption.global_quantities.energy_radiated_electrons_impurities[it]       = float(powers_dict[itime].get('E_imp[MJ]'))*1.e6
    # disruption.global_quantities.energy_radiated_electrons_impurities_halo[it]  = float(powers_dict[itime].get('E_imp_halo[MJ]'))*1.e6
    #disruption.global_quantities.psi_halo_boundary[it]                          = float(psi_dict[itime].get('poloidal_flux_at_halo_boundary[Wb]'))
  
  for itime in range(ntime_plasma):
    disruption.global_quantities.psi_halo_boundary[itime]  = float(psi_dict[itime].get('poloidal_flux_at_halo_boundary[Wb]')) 

  for itime in range(ntime_plasma):
    [tc,it] = find_nearest(np.asarray(time_halo),time_plasma[itime])
    #print(it)
    disruption.halo_currents[itime].active_wall_point.r = float(halo_dict[it]['Rtch[cm]'])*1.e-2
    disruption.halo_currents[itime].active_wall_point.z = float(halo_dict[it]['Ztch[cm]'])*1.e-2
    for ipoint in range(npoints):
      disruption.halo_currents[itime].area[ipoint].start_point.r    = float(halo_dict[it]['Rcw[cm]'][ipoint])*1.e-2
      disruption.halo_currents[itime].area[ipoint].start_point.z    = float(halo_dict[it]['Zcw[cm]'][ipoint])*1.e-2
      disruption.halo_currents[itime].area[ipoint].end_point.r      = float(halo_dict[it]['Rccw[cm]'][ipoint])*1.e-2
      disruption.halo_currents[itime].area[ipoint].end_point.z      = float(halo_dict[it]['Zccw[cm]'][ipoint])*1.e-2
      disruption.halo_currents[itime].area[ipoint].current_halo_pol = float(halo_dict[it]['Ihpol[kA]'][ipoint])*1e3

  # nradial = len(profiles_dict[0]['normalized_psi'])
  # disruption.profiles_1d.resize(ntime_plasma)
  # for itime in range(ntime_profiles):
  #   [tc,it] = find_nearest(np.asarray(time_plasma),time_profiles[itime])
  #   disruption.profiles_1d[it].grid.psi.resize(nradial)
  #   disruption.profiles_1d[it].j_runaways.resize(nradial)
  #   disruption.profiles_1d[it].power_density_conductive_losses.resize(nradial)
  #   disruption.profiles_1d[it].power_density_radiative_losses.resize(nradial)
  #   psi_norm         = profiles_dict[itime].get('normalized_psi')
  #   # denormalize psi profile
  #   psi_axis         = psi_dict[itime].get('poloidal_flux_at_magnetic_axis[Wb]')
  #   psi_bnd_halo     = psi_dict[itime].get('poloidal_flux_at_halo_boundary[Wb]')
  #   psi_LCFS         = psi_dict[itime].get('poloidal_flux_at_plasma_boundary[Wb]')
  #   if (psi_bnd_halo == 0):   # Normalized to halo boundary, when halo is active
  #     psi_bnd = psi_LCFS
  #   else:
  #     psi_bnd = psi_bnd_halo
  #   psi = psi_norm*(psi_bnd-psi_axis) + psi_axis     
  #   disruption.profiles_1d[it].grid.psi   = psi
  #   disruption.profiles_1d[it].j_runaways = profiles_dict[itime]['runaway_current_density[A/m2]']
  #   disruption.profiles_1d[it].power_density_conductive_losses = profiles_dict[itime]['power_density_of_conductive_loss[W/m3]']
  #   disruption.profiles_1d[it].power_density_radiative_losses = profiles_dict[itime]['power_density_of_conductive_loss[W/m3]']

  print('Write output IDSs')
  print('--> Write summary IDS')
  output.put(summary)
  # print('--> Write core_profiles IDS')
  # output.put(core_profiles)
  # print('--> Write core_transport IDS')
  # output.put(core_transport)
  print('--> Write pf_active IDS')
  output.put(pf_active)
  print('--> Write pf_passive IDS')
  output.put(pf_passive)
  print('--> Write equilibrium IDS')
  output.put(equilibrium)
  print('--> Write wall IDS')
  output.put(wall)
  print('--> Write disruption IDS')
  output.put(disruption)
  print('Done.')

  output.close()




