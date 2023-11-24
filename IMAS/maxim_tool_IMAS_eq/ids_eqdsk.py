#system libraries
import sys
import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.use('TkAgg')

#UAL library
import imas

import geqdsk


def Save(filename, eq1, wall1):
  # IMAS uses COCOS11
  # Psi-related values in IMAS are inversed
  cocos_psi = -1
  
  # Plasma current in the ITER coordinates has negative value
  cur_dir = -1
  
  # Toroidal field in the ITER coordinates has negative value
  btor_dir = -1

  ## 2D
  psirz = eq1.time_slice[0].profiles_2d[0].psi # poloidal flux, Wb, [ir, iz]
  r = eq1.time_slice[0].profiles_2d[0].grid.dim1 # radial coordinate, m
  z = eq1.time_slice[0].profiles_2d[0].grid.dim2 # vertical coordinate, m
  
  nr = len(r)
  nz = len(z)
  rleft = r[0]
  zmid = (z[nz-1] + z[0])/2.0
  rdim = r[nr-1] - r[0]
  zdim = z[nz-1] - z[0]
  
  
  ## 1D
  rho = eq1.time_slice[0].profiles_1d.rho_tor_norm # uniform flux grid (normalized toroidal flux)
  pprime = eq1.time_slice[0].profiles_1d.dpressure_dpsi # Pa/Wb
  ffprime = eq1.time_slice[0].profiles_1d.f_df_dpsi # T²*m²/Wb
  pres = eq1.time_slice[0].profiles_1d.pressure # plasma pressure, Pa
  fpol = eq1.time_slice[0].profiles_1d.f # diamagnetic function f=r*Bt, T*m
  qpsi = eq1.time_slice[0].profiles_1d.q
  
  nrho = len(rho)
  ## 0D
  current = eq1.time_slice[0].global_quantities.ip # plasma current, A
  
  simag = eq1.time_slice[0].global_quantities.psi_axis # poloidal flux at magnetic axis, Wb
  sibry = eq1.time_slice[0].global_quantities.psi_boundary # poloidal flux at the plasma boundary, Wb
  
  rmaxis = eq1.time_slice[0].global_quantities.magnetic_axis.r # m
  zmaxis = eq1.time_slice[0].global_quantities.magnetic_axis.z # m
  
  rcentr = eq1.vacuum_toroidal_field.r0 # m
  bcentr = eq1.vacuum_toroidal_field.b0[0] # T
  
  # Precise elongation in this slice
  elong = eq1.time_slice[0].boundary.elongation
  
  ## Limiter outline
  limitr = len(wall1.description_2d[0].limiter.unit[0].outline.r)
  rlim = wall1.description_2d[0].limiter.unit[0].outline.r # m
  zlim = wall1.description_2d[0].limiter.unit[0].outline.z # m
  
  ## Plasma boundary
  nbbbs = len(eq1.time_slice[0].boundary.outline.r)
  rbbbs = eq1.time_slice[0].boundary.outline.r # m
  zbbbs = eq1.time_slice[0].boundary.outline.z # m
  
  
  
  ## Apply coordinate transformations
  
  # Eliminate COCOS11
  pprime *= cocos_psi
  ffprime *= cocos_psi
  psirz *= cocos_psi
  simag *= cocos_psi
  sibry *= cocos_psi
  
  # Plasma current direction
  current *= cur_dir
  pprime *= cur_dir
  ffprime *= cur_dir
  psirz *= cur_dir
  simag *= cur_dir
  sibry *= cur_dir
  
  # Toroidal field direction
  bcentr *= btor_dir
  fpol *= btor_dir
  
  
  ## Units transformations for geqdsk format
  pix2 = 2.0*np.pi
  
  psirz /= pix2
  simag /= pix2
  sibry /= pix2
  
  ffprime *= pix2
  pprime *= pix2
  
  # resample 1d from 1:nrho to 1:nr
  # use uniform grid as used in DINA
  rho_nw = np.zeros(nr)
  for i in range(nr):
    rho_nw[i] = np.sqrt(i/(nr-1))
  
  pprime_nw = np.interp(rho_nw, rho, pprime)
  ffprime_nw = np.interp(rho_nw, rho, ffprime)
  qpsi_nw = np.interp(rho_nw, rho, qpsi)
  fpol_nw = np.interp(rho_nw, rho, fpol)
  pres_nw = np.interp(rho_nw, rho, pres)
  
  
  ## Create dictionary of values and write to the file
  data = {'nr': nr, 'nz':nz,        # Number of horizontal and vertical points
          'r':r, 'z':z,                     # Location of the grid-points
          'rdim':rdim, 'zdim':zdim,         # Size of the domain in meters
          'rcentr':rcentr, 'bcentr':bcentr, # Reference vacuum toroidal field (m, T)
          'rleft':rleft,                  # R of left side of domain
          'zmid':zmid,                      # Z at the middle of the domain
          'rmaxis':rmaxis, 'zmaxis':zmaxis,     # Location of magnetic axis
          'simagx':simag, # Poloidal flux at the axis (Weber / rad)
          'sibdry':sibry, # Poloidal flux at plasma boundary (Weber / rad)
          'current':current, # Plasma current in (A)
          'psirz':psirz,    # Poloidal flux in Weber/rad on grid points
          'fpol':fpol_nw,  # Poloidal current function on uniform flux grid
          'pressure':pres_nw,  # Plasma pressure in Pa=nt/m^2 on uniform flux grid
          'ffprime':ffprime_nw, # 
          'pprime':pprime_nw, # 
          'qpsi':qpsi_nw,  # q values on uniform flux grid
          'nbbbs':nbbbs, 'rbbbs':rbbbs, 'zbbbs':zbbbs, # Plasma boundary
          'limitr':limitr, 'rlim':rlim, 'zlim':zlim} # Wall boundary
  
  #f = "DINA-IMAS.eqdsk"
  geqdsk.write(filename, data)
  
  return
  
  ## Draw picture
  fig_equil = plt.figure()
  plt.clf()
  
  # Axis
  data=plt.plot(rmaxis, zmaxis, 'red')
  #data.set_label('axis,  psi=' + "{:.2f}".format(simag))
  
  # Boundary
  data=plt.plot(rbbbs, zbbbs, 'red', linewidth=1.0)
  
  
  # Contours of psi
  n_levels = 9
  dpsi = (simag - sibry)/n_levels
  vmax = simag
  vmin = sibry
  levels = np.arange(vmin, vmax, dpsi)
  colors = 'blue'
  data0 = plt.contour(r, z, np.transpose(psirz), levels=levels,
    colors=colors, linewidths=0.5, linestyles='solid')
  
  
  # Limiter
  data=plt.plot(rlim, zlim, 'black', linewidth=0.5)
  
  time0 = eq1.time_slice[0].time
  plt.title('t = %f s, Ip = %f MA, kp = %f'%(time0, current*1.e-6, elong))
  plt.xlabel('R, m')
  plt.ylabel('Z, m')
  plt.legend()
  plt.axis('scaled')
  
  
  if False:
    fig_interp = plt.figure()
    plt.clf()
    
    data=plt.plot(rho, pres, 'black', linewidth=0.5)
    data=plt.plot(rho_nw, pres_nw, 'red', linewidth=0.5)
    
    
    fig_interp2 = plt.figure()
    plt.clf()
    
    data=plt.plot(range(nrho), rho, 'black', linewidth=0.5)
    data=plt.plot(range(nr), rho_nw, 'red', linewidth=0.5)
  
  
  plt.show()




