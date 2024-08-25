import numpy as np
import pyvista as pv
import h5py

scale = 0.98
scale_geo = 1.0/1000.0 
dR        = 6.2 * (1-scale)

# Prepare wall that will be wetted
meshes = []
meshes.append( pv.read('Sector1.stl') )
meshes.append( pv.read('Sector2.stl') )
meshes.append( pv.read('Sector3.stl') )
meshes.append( pv.read('Sector4.stl') )
meshes.append( pv.read('Sector5.stl') )
meshes.append( pv.read('Sector6.stl') )
meshes.append( pv.read('Sector7.stl') )
meshes.append( pv.read('Sector8.stl') )
meshes.append( pv.read('Sector9.stl') )
for imesh in range(len(meshes)):
    meshes[imesh] = meshes[imesh].scale([scale_geo, scale_geo, scale_geo], inplace=False) 

meshes.append( pv.read('filter1_meters.stl') )
meshes.append( pv.read('filter2_meters.stl') )

#S cale but centering at the Rgeo
for imesh in range(len(meshes)):
    for ip in range(len(meshes[imesh].points)):
        
        xyz = meshes[imesh].points[ip][:]*scale 
        
        R = (xyz[0]**2 + xyz[1]**2)**0.5
        
        meshes[imesh].points[ip][:] = np.array([xyz[0] + dR*xyz[0]/R, xyz[1] + dR*xyz[1]/R, xyz[2]])


# If you have multiple files you can merge and scale like this (units are assumed to be in meters in the wall input)
mesh  = meshes[0].merge(meshes[1:])

mesh.save('mesh_to_load.stl')
faces = mesh.faces.reshape((-1, 4))[:, 1:]
nodes = mesh.points[faces]
ntriangle = nodes.shape[0]
nodes = np.transpose(nodes,(0,1,2)).ravel() 
 
with h5py.File('wall.h5','w') as h5:
    h5.create_dataset("ntriangle",  (1,),              data=ntriangle, dtype='i4')
    h5.create_dataset("nodes",      (ntriangle*3*3,),  data=nodes,     dtype='f8')
    h5.create_dataset("indices",       data=np.transpose(faces),   dtype='i4')
 
# with h5py.File('wall_to_load.h5','w') as h5:
#     h5.create_dataset("ntriangle",  (1,),              data=ntriangle, dtype='i4')
#     h5.create_dataset("nodes",      (ntriangle*3*3,),  data=nodes,     dtype='f8')
#     h5.create_dataset("indices",       data=np.transpose(faces),   dtype='i4')