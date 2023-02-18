
import imas, os, copy 
import numpy as np
from scipy import interpolate
from scipy.interpolate import griddata
import matplotlib.pyplot as plt


def build_JOREK_boundary(tokamak_name):

  R_scale = 1

  if (tokamak_name == 'ITER'):
    #--------------------close fit to ITER wall
    ellip  = 2.0
    tria_u = 0.55
    tria_l = 0.65
    quad_u = -0.1
    quad_l = 0.15
    n_tht  = 257
    r0     = 6.2  * R_scale
    z0     = 0.1  * R_scale
    a0     = 2.25 * R_scale
  
    #-------------------- contour outside ITER wall
    ellip  = 2.1
    tria_u = 0.58
    tria_l = 0.65
    quad_u = -0.12
    quad_l = -0.
    n_tht  = 257
    r0     = 6.2   * R_scale
    z0     = -0.05 * R_scale
    a0     = 2.34  * R_scale
  
  elif (tokamak_name == 'JET'):
    
    #-------------------- contour outside JET wall
    # blue contour in https://www.jorek.eu/wiki/doku.php?id=eqdsk2jorek.f90
    ellip  = 1.85
    tria_u = 0.4
    tria_l = 0.4
    quad_u = -0.2
    quad_l = -0.2
    n_tht  = 257
    r0     = 2.9  * R_scale
    z0     = 0.1  * R_scale
    a0     = 1.08 * R_scale
  
    #-------------------- contour to avoid too long divertor legs
    # red contour in https://www.jorek.eu/wiki/doku.php?id=eqdsk2jorek.f90
    ellip  = 1.7
    tria_u = 0.4
    tria_l = 0.4
    quad_u = -0.4
    quad_l = -0.2
    n_tht  = 257
    r0     = 2.85 * R_scale
    z0     = 0.15 * R_scale
    a0     = 1.1  * R_scale
  
  elif (tokamak_name == 'DIII-D'):
  
    #-------------------- contour outside DIII-D wall
    ellip  = 1.85
    tria_u = 0.4
    tria_l = 0.4
    quad_u = -0.2
    quad_l = -0.2
    n_tht  = 257
    r0     = 1.7 * R_scale
    z0     = 0.  * R_scale
    a0     = 0.7 * R_scale
    
    #-------------------- Atomic physics JOREK/NIMROD/M3D-C1 benchmark case (paper by B. Lyons)
    ellip  = 1.35/0.7
    tria_u = 0.3
    tria_l = 0.3
    quad_u = 0.
    quad_l = 0.
    n_tht  = 257
    r0     = 1.7 * R_scale
    z0     = 0.  * R_scale
    a0     = 0.7 * R_scale
  
  else:
    print( 'Tokamak name not or wrongly specified, stopping' )
 
  r_bnd = np.zeros( n_tht )
  z_bnd = np.zeros( n_tht )

  n_mid = int(n_tht/2)
 
  for i in range(1, n_mid + 1):
    angle    = 2 * np.pi * float(i-1)/float(n_tht-1)
    r_bnd[i-1] = r0 + a0 * np.cos(angle + tria_u*np.sin(angle) + quad_u*np.sin(2*angle))
    z_bnd[i-1] = z0 + a0 * ellip * np.sin(angle)
  for i in range(n_mid + 1,n_tht+1):
    angle    = 2 * np.pi * float(i-1)/float(n_tht-1)
    r_bnd[i-1] = r0 + a0 * np.cos(angle + tria_l*np.sin(angle) + quad_l*np.sin(2*angle))
    z_bnd[i-1] = z0 + a0 * ellip * np.sin(angle)

  return r_bnd, z_bnd, r0, z0

# Import shot
shot = 105028 
run  = 1      
it   = 211   # Time slice 

# cocos factors
cocos_psi  = 1.0/(2*np.pi)       # Transform to COCOS convention 11 --> 8
cocos_curr = -1.0 
cocos_Bphi = -1.0

input = imas.ids(shot,run,0,0)
input.open_env('public','iterdb','3')
input.equilibrium.get()
input.pf_active.get()
input.core_profiles.get()

#i=0
#for time in input.equilibrium.time:
#    print(str(i)+" "+str(time))
#    i += 1

########## Create boundary of the JOREK domain ###############
boundary_line = build_JOREK_boundary('ITER')
R_bnd = boundary_line[0] 
Z_bnd = boundary_line[1] 
R_geo = boundary_line[2] 
Z_geo = boundary_line[3] 
###############################################################

# Constants
mu0      = 12.566370614359e-7
m_proton = 1.67262192e-27
e_ch     = 1.6021766e-19

# Read 0D parameters
a_min        = input.equilibrium.time_slice[it].boundary.minor_radius
eps          = a_min / R_geo
B_geo        = input.equilibrium.vacuum_toroidal_field.r0 * input.equilibrium.vacuum_toroidal_field.b0[0] / R_geo * cocos_Bphi
xip          = input.equilibrium.time_slice[it].global_quantities.ip           * cocos_curr
psi_axis     = input.equilibrium.time_slice[it].global_quantities.psi_axis     * cocos_psi
psi_boundary = input.equilibrium.time_slice[it].global_quantities.psi_boundary * cocos_psi

# Read 1D profiles  
psi_1d      = input.equilibrium.time_slice[it].profiles_1d.psi                 * cocos_psi
pressure_1d = input.equilibrium.time_slice[it].profiles_1d.pressure
ffprime_1d  = input.equilibrium.time_slice[it].profiles_1d.f_df_dpsi     *(-1.0)      / cocos_psi  
ions        = input.core_profiles.profiles_1d[it].ion

psi_norm    = (psi_1d - psi_axis) / (psi_boundary - psi_axis)

# Read 2D profiles  
R_2d   = input.equilibrium.time_slice[it].profiles_2d[0].r
Z_2d   = input.equilibrium.time_slice[it].profiles_2d[0].z
psi_2d = input.equilibrium.time_slice[it].profiles_2d[0].psi                   * cocos_psi

# Read PF coil currents
coils   = input.pf_active.coil

# Get poloidal flux at the JOREK boundary
Ra = R_2d.flatten()
Za = Z_2d.flatten()
pa = psi_2d.flatten()

psi_bnd = griddata((Ra, Za), pa, (R_bnd, Z_bnd))  

# Get ion density profile
n_tot    = np.zeros( len(ions[0].density))
rho_tot  = np.zeros( len(ions[0].density))

for ion in ions:
  n_tot   += ion.density
  rho_tot += ion.density*ion.element[0].a * m_proton

n0              = n_tot[0]
central_density = n0*1e-20
central_mass    = rho_tot[0] / (n0 * m_proton)
rho_unit        = rho_tot / rho_tot[0]

# Extend profiles into the SOL
n_sol          = 80
n_core         = len(psi_1d)
n_tot          = n_sol + n_core
psin_sol       = 2.0
psin_sep       = psi_norm[-1] + (psi_norm[-1] - psi_norm[-2])

psi_norm_ext    = np.zeros( n_tot )
ffprime_1d_ext  = np.zeros( n_tot )
pressure_1d_ext = np.zeros( n_tot )
rho_1d_ext      = np.zeros( n_tot )

psi_norm_ext[0:n_core]       = psi_norm
ffprime_1d_ext[0:n_core]     = ffprime_1d
pressure_1d_ext[0:n_core]    = pressure_1d
rho_1d_ext[0:n_core]         = rho_unit 
psi_norm_ext[n_core:n_tot]   = np.linspace(psin_sep,psin_sol,num=n_sol) 
ffprime_1d_ext[n_core:n_tot] = np.full( n_sol,  ffprime_1d[-1] ) 
pressure_1d_ext[n_core:n_tot]= np.full( n_sol, pressure_1d[-1] ) 
rho_1d_ext[n_core:n_tot]     = np.full( n_sol,    rho_unit[-1] ) 
rho_1d_ext   = np.full( n_tot,    1.0 ) 

ramp_down   = 0.5*(1 - np.tanh((psi_norm_ext - 1.01)/0.01) )
rho_bnd     = rho_unit[-1]
rho_1d_ext  = (rho_1d_ext-rho_bnd) * ramp_down + rho_bnd
temp_1d_ext = pressure_1d_ext*mu0/rho_1d_ext
T_bnd       = 2 * (e_ch * mu0 * n0 )  # 2 eVs
temp_1d_ext = (temp_1d_ext-T_bnd) * ramp_down + T_bnd

# Export profiles
np.savetxt( 'jorek_density',         np.transpose( [psi_norm_ext, rho_1d_ext]  ))
np.savetxt( 'jorek_temperature',     np.transpose( [psi_norm_ext, temp_1d_ext] ) )
np.savetxt( 'jorek_ffprime',         np.transpose( [psi_norm_ext, ffprime_1d_ext     ] ) )

# Write namelist files
namelist = open('jorek_namelist', 'w')
  
namelist.write( "**********************************************\n")
namelist.write( "* namelist from imas2jorek.py                 *\n")
namelist.write("***********************************************\n")
namelist.write(" B_geo = "+str(B_geo)+"\n")
namelist.write(" R_geo = "+str(R_geo)+"\n")
namelist.write("***********************************************\n")

namelist.write(" &in1"+"\n")

namelist.write("  tstep = 1."+"\n")
namelist.write("  nstep = 0"+"\n")
namelist.write("\n")
namelist.write("  freeboundary       = .f."+"\n")
namelist.write("  resistive_wall     = .f."+"\n")
namelist.write("  freeboundary_equil = .f."+"\n")
namelist.write("\n")
namelist.write("  psi_axis_init = "+str(psi_axis)+"\n")
namelist.write("  amix          = 0.d0"+"\n")
namelist.write("  amix_freeb    = 0.d0"+"\n")
namelist.write("\n")
for i in range(0, len(coils)):
 namelist.write("  pf_coils("+str(i+1)+")%current = "+str(coils[i].current.data[it]*cocos_curr)+"  !"+coils[i].name+"\n")
namelist.write("\n")
namelist.write("  n_R      = 0"+"\n")
namelist.write("  n_Z      = 0"+"\n")
namelist.write("  n_radial = 41"+"\n")
namelist.write("  n_pol    = 64"+"\n")
namelist.write("  n_flux   = 0"+"\n")
namelist.write("  n_tht    = 64"+"\n")
namelist.write("  n_open   = 15"+"\n")
namelist.write("  n_leg    = 15"+"\n")
namelist.write("  n_private = 9"+"\n")
namelist.write("  dPSI_open    = 0.04"+"\n")
namelist.write("  dPSI_private = 0.02"+"\n")
namelist.write(" "+"\n")
namelist.write("  rho_file     = 'jorek_density'"+"\n")
namelist.write("  T_file       = 'jorek_temperature'"+"\n")
namelist.write("  ffprime_file = 'jorek_ffprime'"+"\n")
namelist.write(" "+"\n")

namelist.write(" "+"\n")
namelist.write("  F0 = "+str(R_geo * B_geo)+"\n")
namelist.write("  xpoint = .t."+"\n")
namelist.write("\n")
namelist.write("  D_par     = 0.d0"+"\n")
namelist.write("  D_perp(1) = 1.d-6"+"\n")
namelist.write("  D_perp(2) = 0.85d0"+"\n")
namelist.write("  D_perp(3) = 0.d0"+"\n")
namelist.write("  D_perp(4) = 0.01d0"+"\n")
namelist.write("  D_perp(5) = 2.d0"+"\n")
namelist.write("  ZK_par     = 1.d0"+"\n")
namelist.write("  ZK_perp(1) = 1.d-6"+"\n")
namelist.write("  ZK_perp(2) = 0.85d0"+"\n")
namelist.write("  ZK_perp(3) = 0.d0"+"\n")
namelist.write("  ZK_perp(4) = 0.01d0"+"\n")
namelist.write("  ZK_perp(5) = 2.0d0"+"\n")
namelist.write(""+"\n")
namelist.write("  heatsource     = 0.d0"+"\n")
namelist.write("  particlesource = 0.d0"+"\n")
namelist.write(""+"\n")
namelist.write("  central_mass    = "+str(central_mass)    + "\n")
namelist.write("  central_density = "+str(central_density) + "\n")
namelist.write(""+"\n")
namelist.write("  R_geo = "+str(R_geo)+"\n")
namelist.write("  Z_geo = "+str(Z_geo)+"\n")
namelist.write("  mf = 0 \n")
namelist.write("  n_boundary = "+str(len(R_bnd))+"\n")

for i in range(0, len(R_bnd)):
 namelist.write("R_boundary("+str(i+1)+")="+str(R_bnd[i])+",Z_boundary("+str(i+1)+")="+str(Z_bnd[i])+",psi_boundary("+str(i+1)+")= " +str(psi_bnd[i])+"\n")

namelist.write( "/"+"\n")
namelist.close()


#fig = plt.figure()
#ax = fig.add_subplot(111)
#ax.plot(R_bnd,Z_bnd)
#cax = ax.pcolor(R_2d,Z_2d,psi_2d, cmap='rainbow')
#cbar = fig.colorbar(cax)
#plt.show()
