import numpy as np
import matplotlib.pyplot as plt
import matplotlib
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
from raysect.core.math.function.float.function2d.interpolate import Interpolator2DMesh, Discrete2DMesh

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
#    try:
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

        # subdivide the grid
#        nr_nodes_new = nr_nodes*n_sub**2
#        nr_faces_new = nr_faces*n_sub**2

#        sp = np.linspace(0.0, 1.0, n_sub)
#        tp = np.linspace(0.0, 1.0, n_sub)

#        ref_elm_list = []
#        nvertex = 4
#
#
#        # Create refined list of elements
#        for face_id in range(nr_faces):
#
#            for is1 in range(0,n_sub):
#                s1 = float(is1)  /float(n_sub)
#                s2 = float(is1+1)/float(n_sub)
#                for it1 in range(0,n_sub):
#                    t1 = float(it1  )/float(n_sub)
#                    t2 = float(it1+1)/float(n_sub)
#
#                    # fill sub elm
#                    elm        = elm_class()
#                    elm.parent = face_id
#
#                    RZ1 = interp_RZ(node_list, element_list, face_id, s1, t1)
#                    RZ2 = interp_RZ(node_list, element_list, face_id, s1, t2)
#                    RZ3 = interp_RZ(node_list, element_list, face_id, s2, t1)
#                    RZ4 = interp_RZ(node_list, element_list, face_id, s2, t2)
#
#                    elm.RZ = np.array( [RZ1, RZ2, RZ3, RZ4] )
#
#                    ref_elm_list.append( elm )
#
#        vertex_node_array = np.zeros(shape=(nvertex, nr_faces_new))
#
#        node_RZ = np.zeros(shape=(nr_nodes_new, 2)) 
#
#        print(nr_nodes_new, nr_faces_new, nvertex)
#        print("connect")            
        # Brute force creation of connectivity matrix
#        for inode in range(0, nr_nodes_new):
#            node_RZ[inode] = np.array([-1,-1])

#            for ielm in range(0,nr_faces_new):
#                elm = ref_elm_list[ielm]
#                for iv in range(0, nvertex):
#                    pass
#                    if (node_RZ[inode,0]< 0):
#                        node_RZ[inode] = elm.RZ[iv]
#                        vertex_node_array[iv, ielm] = inode 


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
#            q_to_raysect = 1.0
            f.write("SCALARS Q_process"+str(1)+"_ion_"+str(1)+"_grid_subset_"+str(1)+"[W/m^2] float 1\n")
            f.write("LOOKUP_TABLE default\n")
            for q in q_to_raysect:
                f.write(str(q)+"\n")
           
#    except:
#        logging.exception("Error reading shot "+str(shot)+":")


print("end")

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
radiation_interp = Discrete2DMesh(vertex_coords, triangles, tria_values,
                                  limit=False, default_value=0)

# Extend r-z mesh in toroidal direction
rad_function_3d = AxisymmetricMapper(radiation_interp)
radiation_emitter = VolumeTransform(RadiationFunction(rad_function_3d))

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
output_basename = "example"
export_vtk(bolometer1801, output_basename + '.vtk', triangle_data=triangle_data)

