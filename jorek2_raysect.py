import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from sys import exit
from interp_jorek import *
from raysect.primitive import Sphere, import_vtk, export_vtk, Box
from raysect.optical.observer import MeshPixel, MeshCamera, PowerPipeline0D, PowerPipeline1D, MonoAdaptiveSampler1D
from raysect.core import translate, rotate_x
from raysect.optical import World
from raysect.optical.material.emitter import UnityVolumeEmitter
from raysect.optical.material import AbsorbingSurface, VolumeTransform
from raysect.core import Point3D, translate
from scipy import interpolate
import scipy as sc
from cherab.core.math import sample2d, AxisymmetricMapper
from cherab.tools.emitters import RadiationFunction
from raysect.primitive import Cylinder, Subtract
from raysect.optical.library import Copper, Beryllium
from cherab.core.math import Constant3D, ConstantVector3D, sample3d, AxisymmetricMapper
from raysect.optical.material import Lambert, UniformSurfaceEmitter, Roughen
from raysect.optical.material import RoughConductor
from raysect.optical import InterpolatedSF
from raysect.core.math.function.float.function3d.interpolate import Discrete3DMesh

class elm_class:
    def __init__(self):
        self.vertex_ind = np.empty(0)
        self.RZ         = np.empty(0)
        self.parent     = 0

def quads_to_tris(quads, cell_values_quads):
    tris = [[None for j in range(3)] for i in range(2*len(quads))]
    cell_values_tris = [0 for i in range(2*len(cell_values_quads))]
    for i in range(len(quads)):
        j = 2*i
        cell_values_tris[j] = cell_values_quads[i]
        cell_values_tris[j+1] = cell_values_quads[i]
        n0 = quads[i][0]
        n1 = quads[i][1]
        n2 = quads[i][2]
        n3 = quads[i][3]
        tris[j][0] = n0
        tris[j][1] = n1
        tris[j][2] = n2
        tris[j + 1][0] = n2
        tris[j + 1][1] = n3
        tris[j + 1][2] = n0
    return tris, cell_values_tris


def to_cart_coord(R,Z,phi):
    x =  R * np.cos(phi)
    y = -R * np.sin(phi)

    return np.array( [x,y,Z] )


import imas
import logging


username = "artolaj"
device = "di_SPI"

shot_list = [111111]
run = 1

n_sub = 1

logging.basicConfig(level=logging.DEBUG, filename='./logfile.log')
shot_dict = {} 
for shot in shot_list:
    try:
        x = imas.ids(shot, run)
        x.open_env(username, device, "3")
        x.radiation.get()
        # 0d
        pts = []

        # Definition of variables for Raysect
        q_to_raysect = [] # values on quadrangles -> input for raysect
        r = [] # r coordinate -> input for raysect
        z = [] # z coordinate -> input for raysect
        quadrangles = [] # quadrangle node IDs -> input for raysect

        nr_nodes = len(x.radiation.grid_ggd[0].space[0].objects_per_dimension[0].object.array)
        nr_faces = len(x.radiation.grid_ggd[0].space[0].objects_per_dimension[2].object.array)

        for node_id in range(nr_nodes):
            pts.append(x.radiation.grid_ggd[0].space[0].objects_per_dimension[0].object[node_id].geometry)

        # 2d
        faces = []
        for face_id in range(nr_faces):
            faces.append(x.radiation.grid_ggd[0].space[0].objects_per_dimension[2].object[face_id].nodes)

        # JOREK type namelist
        element_list  = reconstruct_element_list( x.radiation.grid_ggd[0] )
        node_list     = reconstruct_node_list(    x.radiation.grid_ggd[0] )
        val_coeff     = x.radiation.process[0].ggd.array[0].ion[0].emissivity.array[0].coefficients
 
        # Calculate emissivity in element center
        q_to_raysect = np.zeros(nr_faces)
        for face_id in range(nr_faces):
            q_to_raysect[face_id] =  interp_val(val_coeff, node_list, element_list, face_id, 0.5, 0.5, 0.0)


        # Get R and Z vectors              
        for node in node_list.nodes:
            r.append(  node.x[0,0]  )
            z.append(  node.x[0,1]  )


        print("done")            
    
        with open("jorek_radiation_shot"+str(shot)+"_run"+str(run)+".vtk", "w") as f:
            f.write("# vtk DataFile Version 2.0\n")
            f.write("vtk output\n")
            f.write("ASCII\n")
            f.write("\n")
            f.write("DATASET UNSTRUCTURED_GRID\n")
            f.write("POINTS "+str(nr_nodes)+" float\n")
            for pt in enumerate(pts):
                f.write(str(pt[1][0])+" "+str(pt[1][1])+" 0\n")
                # Append to r and z for Raysect
                r.append(pt[1][0])
                z.append(pt[1][1])
        
            f.write("\n")
            f.write("CELLS "+str(nr_faces)+" "+str(nr_faces*5)+"\n")

            quadrangle_center = []
            for cell in faces:
                f.write("4 "+str(cell[0]-1)+" "+
                        str(cell[1]-1)+" "+
                        str(cell[2]-1)+" "+
                        str(cell[3]-1)+"\n")
                # Append quadrangle IDs for Raysect
                quadrangles.append([cell[0]-1, cell[1]-1, cell[2]-1, cell[3]-1])

            f.write("\n")

            f.write("CELL_TYPES "+str(nr_faces)+"\n")
            for cell in range(len(faces)):
                f.write("9\n")

            f.write("\n")
            f.write("CELL_DATA "+str(nr_faces)+"\n")
            f.write("SCALARS Q_process"+str(1)+"_ion_"+str(1)+"_grid_subset_"+str(1)+"[W/m^2] float 1\n")
            f.write("LOOKUP_TABLE default\n")
            for q in q_to_raysect:
                f.write(str(q)+"\n")
           
    except:
        logging.exception("Error reading shot "+str(shot)+":")


print("end")

# Calculate tetrahedra mesh
print("doing tetrahedra...")

N_phi       = 32                 # Number of toroidal points
N_elm       = element_list.n_elm   # Number of poloidal quadrilateral elements
N_pol_nodes = node_list.n_nodes    # Number of nodes in the poloidal plane
N_vertex    = 4                    # Number of vertices of each element
N_loc_tet   = 5                    # Every rectangular prism is divided into 5 tetrahedra 

print('Npol_nodes = ' + str(N_pol_nodes))
print('Nelem      = ' + str(N_elm))

n_tetra_nodes   = N_phi * N_pol_nodes
n_tetra         = N_phi * N_elm * N_loc_tet
nodes_xyz       = np.zeros( shape=(n_tetra_nodes, 3) )
tetra_ind       = np.zeros( shape=(n_tetra, 4), dtype=int )
tetra_val       = np.zeros( n_tetra )

# Fill global node list
for i_pol_node in range(0, N_pol_nodes):

    Rnode = node_list.nodes[i_pol_node].x[0,0]
    Znode = node_list.nodes[i_pol_node].x[0,1]

    # Go toroidally
    for i_phi in range(0,N_phi):

        phi = float(i_phi  ) / float(N_phi) * 2.0*np.pi 

        i_node_glob = i_phi * N_pol_nodes + i_pol_node 

        nodes_xyz[i_node_glob] = to_cart_coord( Rnode, Znode, phi )


# Go over elements and form prisms
i_tetra = -1

for i_elm in range(0, N_elm):

    # Go toroidally
    for i_phi in range(0,N_phi):

        # Global node indices
        i0 =element_list.elems[i_elm].vertex_ind[0]-1 
        i1 =element_list.elems[i_elm].vertex_ind[1]-1 
        i2 =element_list.elems[i_elm].vertex_ind[2]-1 
        i3 =element_list.elems[i_elm].vertex_ind[3]-1 


        # Bottom face
        B0_ind = i_phi *N_pol_nodes + i0 
        B1_ind = i_phi *N_pol_nodes + i1 
        B2_ind = i_phi *N_pol_nodes + i2
        B3_ind = i_phi *N_pol_nodes + i3

        if (i_phi == N_phi - 1):
            i_phiT = 0
        else:
            i_phiT = i_phi +1 

        # Top face
        T0_ind = i_phiT*N_pol_nodes + i0
        T1_ind = i_phiT*N_pol_nodes + i1
        T2_ind = i_phiT*N_pol_nodes + i2
        T3_ind = i_phiT*N_pol_nodes + i3

        # Emissivity value at prism center
        phi_mid = 0.5 * ( float(i_phi)/float(N_phi) + float(i_phi+1)/float(N_phi) ) * 2.0 * np.pi 
        val     = interp_val(val_coeff, node_list, element_list, i_elm, 0.5, 0.5, phi_mid)

        # Define the 5 tetrahedra and their indices
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


with open("jorek3d_radiation_shot"+str(shot)+"_run"+str(run)+".vtk", "w") as f:
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
    f.write("SCALARS Q_process"+str(1)+"_ion_"+str(1)+"_grid_subset_"+str(1)+"[W/m^3] float 1\n")
    f.write("LOOKUP_TABLE default\n")
    for val in tetra_val:
        f.write(str(val)+"\n")
 
print("tetrahedra done")
#############################



# Convert quadrangles into triangles (Discrete2DMesh takes triangles
# and triangle values)
triangles, tria_values = quads_to_tris(quadrangles, q_to_raysect)

# Define min and max values of grid (to later construct Raysect cylinders)
r_min = min(r)
r_max = max(r)
z_min = min(z)
z_max = max(z)
height = max(z)-min(z)

vertex_coords = np.array([[r[i], z[i]] for i in range(len(r))])

# Define r-z mesh in Raysect
print('Creating 3D interpolation function')
radiation_interp = Discrete3DMesh(nodes_xyz, tetra_ind, tetra_val,
                                  limit=False, default_value=0)

print('Start RaySect...')
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
    pixel_samples=1) # number of rays to launch from each triangle

# Render and compute
print('Observing the Sphere...')
render_pass = 0
while (not camera.render_complete) and (render_pass < 100):
    render_pass += 1
    print('Render pass {}:'.format(render_pass))
    camera.observe()
    # obtain frame from pipeline

frame = camera.pipelines[0].frame
power_density = frame.mean / camera.collection_areas
error = frame.errors() / camera.collection_areas
triangle_data = {'PowerDensity': power_density, 'PowerDensityError': error}

# Define output name and export result
output_basename = "JOREK3D_bolomer_power"
export_vtk(bolometer1801, output_basename + '.vtk', triangle_data=triangle_data)
