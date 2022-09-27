import numpy as np
from sys import exit
from interp_jorek import *
from raysect.primitive import  import_vtk, export_vtk, import_stl
from raysect.optical.observer import MeshCamera, PowerPipeline1D, MonoAdaptiveSampler1D
from raysect.optical import World
from raysect.optical.material import AbsorbingSurface, VolumeTransform
from raysect.core import translate
from cherab.tools.emitters import RadiationFunction
from raysect.primitive import Cylinder, Subtract
from raysect.core.math.function.float.function3d.interpolate import Discrete3DMesh
import imas
import logging
import time
import re

# Function to transform cylindrical to cartesian coordinates
def to_cart_coord(R,Z,phi):
    x =  R * np.cos(phi)
    y = -R * np.sin(phi)
    return np.array( [x,y,Z] )

# Main program
def main():

    # Input parameters
    username     = "artolaj"
    device       = "di_SPI"
    shot_list    = [111112]
    run          = 1
    N_phi        = 64    # Number of toroidal points for 3D mesh
    N_passes     = 60    # Number of camera render passes in RaySect
    N_pixel      = 1000  # Number of pixel samples
    rot_plasma   = 0     # Rotates the plasma by this angle 
    surface_file = "../Sector1.stl"  # Surface where radiation is calculated
    scale_vtk    = 0.001            # 0.001 for FW pannels
    flip_norm    = False            # Flips normals

    # Inputs for energy conservation check (< 0.5% relative error)
    #N_phi       = 64     # Number of toroidal points for 3D mesh
    #N_passes    = 30   # Number of camera render passes in RaySect
    #N_pixel     = 3000 # Number of pixel samples
    #rot_plasma  = -2.8   # Rotates the plasma by this angle 
    #surface_file= "./observer_torus.vtk"  # Surface where radiation is calculated
    #scale_vtk   = 1.16  # 0.001 for FW pannels, 1.16 for torus observer
    #flip_norm   = True  # Flips normals
   
    # Some hardcoded parameters
    N_vertex   = 4      # Number of vertices of each poloidal element
    N_loc_tet  = 5      # Every rectangular prism is divided into 5 tetrahedra 
    i_grid     = 0      # Time index of grid (In JOREK the grid does not vary over time)     

    # Check that the radiated surface file has the right format
    if ".stl" in surface_file:
        f_format = "stl"
        fname    = re.search('/(.*).stl', surface_file).group(1)
    elif ".vtk" in surface_file:
        f_format = "vtk"
        fname    = re.search('/(.*).vtk', surface_file).group(1)
    else:
        print( "File format for surface not known!")
        exit()
    
    logging.basicConfig(level=logging.DEBUG, filename='./output.log')

    # Print parameters
    logging.info('Input parameters')
    logging.info('   N_phi    = '+str(N_phi))
    logging.info('   N_passes = '+str(N_passes))
    logging.info('   N_pixel  = '+str(N_pixel))

    for shot in shot_list:
    
        logging.info("Reading shot = %s, run = %s from database = %s of user = %s "%(shot,run,device,username))
        x = imas.ids(shot, run)
        x.open_env(username, device, "3")
        x.radiation.get()

        # Get JOREK grid
        # Reconstruct JOREK element/node lists
        element_list  = reconstruct_element_list( x.radiation.grid_ggd[i_grid] )
        node_list     = reconstruct_node_list(    x.radiation.grid_ggd[i_grid] )
 
        # Definition of variables for Raysect
        r = [] # r coordinate -> input for raysect
        z = [] # z coordinate -> input for raysect

        # Get poloidal R and Z vectors              
        for node in node_list.nodes:
            r.append(  node.x[0,0]  )
            z.append(  node.x[0,1]  )
       
        N_elm       = element_list.n_elm   # Number of poloidal quadrilateral elements
        N_pol_nodes = node_list.n_nodes    # Number of nodes in the poloidal plane
 
        logging.info('N_pol_nodes = ' + str(N_pol_nodes))
        logging.info('N_elm       = ' + str(N_elm))
           
        n_tetra_nodes   = N_phi * N_pol_nodes                         # Number of 3D mesh nodes
        n_tetra         = N_phi * N_elm * N_loc_tet                   # Number of tetrahedra
        nodes_xyz       = np.zeros( shape=(n_tetra_nodes, 3) )        # A list with the coordinates of the 3D mesh nodes
        tetra_ind       = np.zeros( shape=(n_tetra, 4), dtype=int )   # Tetrahedra connectivity matrix
        tetra_val       = np.zeros( n_tetra )                         # Values at tetrahedra centres
        
        # Fill up global node list coordinates
        for i_pol_node in range(0, N_pol_nodes):
        
            Rnode = node_list.nodes[i_pol_node].x[0,0]
            Znode = node_list.nodes[i_pol_node].x[0,1]
        
            # Go toroidally
            for i_phi in range(0,N_phi):
        
                phi         = float(i_phi) / float(N_phi) * 2.0*np.pi + rot_plasma  # Toroidal angle of node 
                i_node_glob = i_phi * N_pol_nodes + i_pol_node          # Global 3D index of node
        
                nodes_xyz[i_node_glob] = to_cart_coord( Rnode, Znode, phi )

        # Time loop: Go over time slices in the IDS
        N_times = len(x.radiation.time)
        logging.info( "Number of time slices = " + str(N_times))

        for i_time in range(0, N_times):

            time_now = x.radiation.time[i_time]
            logging.info('Time  = ' + str(time_now) + ' s ')

            try:
                val_coeff     = x.radiation.process[0].ggd.array[i_time].ion[0].emissivity.array[0].coefficients
            except:
                logging.exception("Error reading radiation coefficients for time slice  "+str(i_time)+":")
                continue  # Skip this shot


            # Quantities needed for integration over elements
            n_dof    = element_list.n_dof   
            n_vertex = element_list.n_vertex
            n_tor    = len( val_coeff[:,0] )

            dV_gauss = np.zeros( shape=( ngauss,   ngauss) ) # volume element at gaussian points
            coeff    = np.zeros( shape=( n_vertex, n_dof, n_tor  ) )
            coeff_R  = np.zeros( shape=( n_vertex, n_dof ) )
            coeff_Z  = np.zeros( shape=( n_vertex, n_dof ) )
            val      = np.zeros(N_phi)

            # Get poloidal and toroidal basis at integration points
            basis_gauss   = get_basis_at_gaussian(n_vertex,n_dof)
            basis_tor_all = get_toroidal_basis(N_phi, n_tor, node_list.n_period)
 
            # Go over elements and form prisms and tetrahedra
            logging.info("Computing 3D tetrahedra mesh...")
            i_tetra = -1
            t_start = time.time()        
            for i_elm in range(0, N_elm):

                ######################### Find average emissitivity by integrating the element ######################### 
                sizes    = element_list.elems[i_elm].size

                for kv in range(0, element_list.n_vertex):  
                    iv             = element_list.elems[i_elm].vertex_ind[kv] - 1
                    coeff_R[kv,:]  = node_list.nodes[iv].x[:,0] * sizes[kv,:] 
                    coeff_Z[kv,:]  = node_list.nodes[iv].x[:,1] * sizes[kv,:]

                    for idof in range(0, n_dof):
                        icoeff            =  iv + idof*node_list.n_nodes 
                        coeff[kv,idof,:]  =  val_coeff[:,icoeff]

                # Integrate over the element to calculate average emissivity
                vol = 0.0;  val = 0.0;
                for ms in range(ngauss):
                    for mt in range(ngauss):
                        out             = interp_RZ_fast(sg[ms], sg[mt], coeff_R, coeff_Z)      # Get R, Z and derivatives for Jacobian
                        dV_gauss[ms,mt] = (out[1]*out[5] - out[2]*out[4])*out[0]*wg[ms]*wg[mt]  # Jacobian*R*gauss_wemshts
                        vol             = vol + dV_gauss[ms,mt] 
                        val = val + np.einsum('ijk,ij,ij,lk->l', coeff, sizes, basis_gauss[ms,mt], basis_tor_all) * dV_gauss[ms,mt]

                val = val / vol          
                ######################### End finding average emissitivity ############################################# 

                # Go toroidally to form tetrahedra
                for i_phi in range(0,N_phi):
            
                    # Global node indices
                    i0 = element_list.elems[i_elm].vertex_ind[0]-1 
                    i1 = element_list.elems[i_elm].vertex_ind[1]-1 
                    i2 = element_list.elems[i_elm].vertex_ind[2]-1 
                    i3 = element_list.elems[i_elm].vertex_ind[3]-1 
            
                    # Bottom face (quadrilateral element at phi)
                    B0_ind = i_phi *N_pol_nodes + i0 
                    B1_ind = i_phi *N_pol_nodes + i1 
                    B2_ind = i_phi *N_pol_nodes + i2
                    B3_ind = i_phi *N_pol_nodes + i3
            
                    i_phiT = i_phi +1 
                    if (i_phi == N_phi-1):
                        i_phiT = 0
            
                    # Top face (quadrilateral element at phi + dphi)
                    T0_ind = i_phiT*N_pol_nodes + i0
                    T1_ind = i_phiT*N_pol_nodes + i1
                    T2_ind = i_phiT*N_pol_nodes + i2
                    T3_ind = i_phiT*N_pol_nodes + i3
            
                    # Define the 5 tetrahedra and their indices
                    # See https://www.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/48509/versions/3/previews/COMP_GEOM_TLBX/html/Divide_hypercube_5_simplices_3D.html
                    # Tetrahedron 0 = 3,1,0,0'
                    i_tetra            = i_tetra + 1 
                    tetra_ind[i_tetra] = [B3_ind, B1_ind, B0_ind, T0_ind]     
                    tetra_val[i_tetra] = val[i_phi]  
            
                    # Tetrahedron 1 = 3,1,2,2'
                    i_tetra            = i_tetra + 1 
                    tetra_ind[i_tetra] = [B3_ind, B1_ind, B2_ind, T2_ind]    
                    tetra_val[i_tetra] = val[i_phi]  
            
                    # Tetrahedron 2 = 3,3',0',2'
                    i_tetra            = i_tetra + 1 
                    tetra_ind[i_tetra] = [B3_ind, T3_ind, T0_ind, T2_ind]     
                    tetra_val[i_tetra] = val[i_phi]  
            
                    # Tetrahedron 3 = 0',1,1',2'
                    i_tetra            = i_tetra + 1 
                    tetra_ind[i_tetra] = [T0_ind, B1_ind, T1_ind, T2_ind]     
                    tetra_val[i_tetra] = val[i_phi]  
            
                    # Tetrahedron 4 = 3, 1, 2', 0'
                    i_tetra            = i_tetra + 1 
                    tetra_ind[i_tetra] =  [B3_ind, B1_ind, T2_ind, T0_ind]     
                    tetra_val[i_tetra] = val[i_phi]  
           
         
            # Export 3D radiation to VTK 
            with open("jorek3d_radiation_shot"+str(shot)+"_run"+str(run)+"_timeslice_"+str(i_time)+".vtk", "w") as f:
                f.write("# vtk DataFile Version 2.0\n")
                f.write("vtk output\n")
                f.write("ASCII\n")
                f.write("\n")
                f.write("DATASET UNSTRUCTURED_GRID\n")
                f.write("\n")
                f.write("FIELD FieldData 1\n")
                f.write("TIME 1 1 double\n")
                f.write(str(time_now)+"\n")
                f.write("\n")
                f.write("POINTS "+str(n_tetra_nodes)+" float\n")
                for node in nodes_xyz:
                    f.write(str(node[0])+" "+str(node[1])+" "+str(node[2])+" \n") 
            
                f.write("\n")
                f.write("CELLS "+str(n_tetra)+" "+str(n_tetra*5)+"\n")
            
                for tetra in tetra_ind:
                    f.write("4 "+str(tetra[0])+" "+
                             str(tetra[1])+" "+
                             str(tetra[2])+" "+
                             str(tetra[3])+"\n")
            
                f.write("\n")
            
                f.write("CELL_TYPES "+str(n_tetra)+"\n")
                for cell in range(0, n_tetra):
                    f.write("10\n")
            
                f.write("\n")
                f.write("CELL_DATA "+str(n_tetra)+"\n")
                f.write("SCALARS Radiation[W/m^3] float 1\n")
                f.write("LOOKUP_TABLE default\n")
                for val in tetra_val:
                    f.write(str(val)+"\n")
             
            t_end = time.time()        
            logging.info("Time to compute tetrahedra and values = "+str(t_end-t_start))
            #############################
            
            
            # Define min and max values of grid (to later construct Raysect cylinders)
            r_min  = min(r)
            r_max  = max(r)
            z_min  = min(z)
            z_max  = max(z)
            height = max(z)-min(z)
            
            # Create 3D interpolation function in Raysect
            logging.info('Creating 3D interpolation function...')
            t_start = time.time()        
            radiation_interp = Discrete3DMesh(nodes_xyz, tetra_ind, tetra_val,
                                              limit=False, default_value=0)
            t_end = time.time()        
            logging.info("Time to compute 3D interpolation function = "+str(t_end-t_start))
            
            logging.info('Start RaySect...')
            radiation_emitter = VolumeTransform(RadiationFunction(radiation_interp))
            
            # Define world in Raysect
            world = World()
            
            # Define two cylinders
            cyl1 = Cylinder(r_max, height, transform=translate(0, 0, z_min), parent=world)
            cyl2 = Cylinder(r_min, height, transform=translate(0, 0, z_min), parent=world)
            
            # Substract cylinders and get plasma body
            emitter = Subtract(cyl1, cyl2, parent=world)
            
            # Append plasma radiation from SOLPS to Raysect plasma body 
            emitter.material = radiation_emitter
            
            # Import vtk mesh (here bolometer housing as example)
            if f_format == "vtk":
                bolometer1801 = import_vtk(surface_file, scaling=scale_vtk, parent=world, flip_normals=flip_norm)
            elif f_format == "stl":
                bolometer1801 = import_stl(surface_file, scaling=scale_vtk, parent=world, flip_normals=flip_norm)
            else:
                logging.info("Unknown format for surface file")
    
            bolometer1801.material = AbsorbingSurface()
            
            # Define Raysect power pipeline
            power = PowerPipeline1D()
            
            # Define Raysect camera
            sampler = MonoAdaptiveSampler1D(power, fraction=0.01, ratio=25.0, min_samples=500, cutoff=0.01)
            camera = MeshCamera(
                bolometer1801,
                surface_offset=1e-6,  # launch rays 1mm off surface to avoid
                                      # intersection with absorbing bolometer1801
                pipelines=[power],
                frame_sampler=sampler,
                parent=world,
                spectral_bins=1,
                min_wavelength=1,
                max_wavelength=2,
                pixel_samples=N_pixel) # number of rays to launch from each triangle
            
            # Render and compute
            logging.info('Start rendering...')
            t_start     = time.time()        
            render_pass = 0
            while (not camera.render_complete) and (render_pass < N_passes):
                render_pass += 1
                logging.info('Render pass {}:'.format(render_pass))
                camera.observe()
             
            # obtain frame from pipeline            
            frame = camera.pipelines[0].frame
            power_density = frame.mean / camera.collection_areas
            error = frame.errors() / camera.collection_areas
            triangle_data = {'PowerDensity': power_density, 'PowerDensityError': error}
            
            # Define output name and export result
            output_basename = "rad_power_" + fname +"_shot"+str(shot)+"_run"+str(run)+"_timeslice_"+str(i_time)
            export_vtk(bolometer1801, output_basename + '.vtk', triangle_data=triangle_data)

            t_end = time.time()        
            logging.info("Time in RaySect loop = "+str(t_end-t_start))

if __name__ == "__main__":
    main()

