import numpy as np
from sys import exit
from interp_jorek import *
from raysect.primitive import  import_vtk, export_vtk
from raysect.optical.observer import MeshCamera, PowerPipeline1D, MonoAdaptiveSampler1D
from raysect.optical import World
from raysect.optical.material import AbsorbingSurface, VolumeTransform
from raysect.core import translate
from cherab.tools.emitters import RadiationFunction
from raysect.primitive import Cylinder, Subtract
from raysect.core.math.function.float.function3d.interpolate import Discrete3DMesh
import imas
import logging

# Function to transfer from cylindrical to cartesian coordinates
def to_cart_coord(R,Z,phi):
    x =  R * np.cos(phi)
    y = -R * np.sin(phi)
    return np.array( [x,y,Z] )

# Main program
def main():

    # Input parameters
    username   = "artolaj"
    device     = "di_SPI"
    shot_list  = [111111]
    run        = 1
    N_phi      = 32     # Number of toroidal points for 3D mesh
    N_passes   = 10     # Number of camera render passes in RaySect
    N_pixel    = 10     # Number of pixel samples
    
    # Some hardcoded parameters
    N_vertex   = 4      # Number of vertices of each poloidal element
    N_loc_tet  = 5      # Every rectangular prism is divided into 5 tetrahedra 
    i_grid     = 0      # Time index of grid (In JOREK the grid does not vary over time)     
    
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
        
                phi         = float(i_phi) / float(N_phi) * 2.0*np.pi   # Toroidal angle of node 
                i_node_glob = i_phi * N_pol_nodes + i_pol_node          # Global 3D index of node
        
                nodes_xyz[i_node_glob] = to_cart_coord( Rnode, Znode, phi )

        # Time loop: Go over time slices in the IDS
        N_times = len(x.radiation.time)
        logging.info( "Number of time slices = " + str(N_times))

        for i_time in range(0, N_times):

            logging.info('Time  = ' + str(x.radiation.time[i_time]) + ' s ')

            try:
                val_coeff     = x.radiation.process[0].ggd.array[i_time].ion[0].emissivity.array[0].coefficients
            except:
                logging.exception("Error reading radiation coefficients for time slice  "+str(i_time)+":")
                continue  # Skip this shot

            # Go over elements and form prisms and tetrahedra
            logging.info("Computing 3D tetrahedra mesh...")
            i_tetra = -1
            
            for i_elm in range(0, N_elm):
            
                # Go toroidally
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
            
                    # Emissivity value at prism center
                    phi_mid = 0.5 * ( float(i_phi)/float(N_phi) + float(i_phi+1)/float(N_phi) ) * 2.0 * np.pi 
                    val     = interp_val(val_coeff, node_list, element_list, i_elm, 0.5, 0.5, phi_mid)
            
                    # Define the 5 tetrahedra and their indices
                    # See https://www.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/48509/versions/3/previews/COMP_GEOM_TLBX/html/Divide_hypercube_5_simplices_3D.html
                    # Tetrahedron 0 = 3,1,0,0'
                    i_tetra            = i_tetra + 1 
                    tetra_ind[i_tetra] = [B3_ind, B1_ind, B0_ind, T0_ind]     
                    tetra_val[i_tetra] = val  
            
                    # Tetrahedron 1 = 3,1,2,2'
                    i_tetra            = i_tetra + 1 
                    tetra_ind[i_tetra] = [B3_ind, B1_ind, B2_ind, T2_ind]    
                    tetra_val[i_tetra] = val  
            
                    # Tetrahedron 2 = 3,3',0',2'
                    i_tetra            = i_tetra + 1 
                    tetra_ind[i_tetra] = [B3_ind, T3_ind, T0_ind, T2_ind]     
                    tetra_val[i_tetra] = val  
            
                    # Tetrahedron 3 = 0',1,1',2'
                    i_tetra            = i_tetra + 1 
                    tetra_ind[i_tetra] = [T0_ind, B1_ind, T1_ind, T2_ind]     
                    tetra_val[i_tetra] = val  
            
                    # Tetrahedron 4 = 3, 1, 2', 0'
                    i_tetra            = i_tetra + 1 
                    tetra_ind[i_tetra] =  [B3_ind, B1_ind, T2_ind, T0_ind]     
                    tetra_val[i_tetra] = val  
           
         
            # Export 3D radiation to VTK 
            with open("jorek3d_radiation_shot"+str(shot)+"_run"+str(run)+"_timeslice_"+str(i_time)+".vtk", "w") as f:
                f.write("# vtk DataFile Version 2.0\n")
                f.write("vtk output\n")
                f.write("ASCII\n")
                f.write("\n")
                f.write("DATASET UNSTRUCTURED_GRID\n")
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
             
            logging.info("tetrahedra done")
            #############################
            
            
            # Define min and max values of grid (to later construct Raysect cylinders)
            r_min  = min(r)
            r_max  = max(r)
            z_min  = min(z)
            z_max  = max(z)
            height = max(z)-min(z)
            
            # Create 3D interpolation function in Raysect
            logging.info('Creating 3D interpolation function...')
            radiation_interp = Discrete3DMesh(nodes_xyz, tetra_ind, tetra_val,
                                              limit=False, default_value=0)
            
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
            bolometer1801 = import_vtk("./simplified_bolometer_dense_mesh.vtk",
                                      scaling=0.001, parent=world)
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
            output_basename = "bolometer_power_shot"+str(shot)+"_run"+str(run)+"_timeslice_"+str(i_time)+".vtk"
            export_vtk(bolometer1801, output_basename + '.vtk', triangle_data=triangle_data)


if __name__ == "__main__":
    main()

