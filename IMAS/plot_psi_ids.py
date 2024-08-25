import imas
import matplotlib
import matplotlib.pyplot as plt
import argparse
import numpy as np
from scipy.interpolate import interp2d

vv = np.loadtxt('/home/ITER/artolaj/ITER_components/ITER_vessel_inner.txt')
fw = np.loadtxt('/home/ITER/artolaj/ITER_components/1st_wall.txt')

# MANAGEMENT OF INPUT ARGUMENTS
# ------------------------------
parser = argparse.ArgumentParser(description=\
        '---- Display scenario')
parser.add_argument('-s','--shot',help='Shot number', required=True,type=int)
parser.add_argument('-r','--run',help='Run number',required=True,type=int)
parser.add_argument('-u','--user_or_path',help='User or absolute path name where the data-entry is located', required=False)
parser.add_argument('-d','--database',help='Database name where the data-entry is located', required=False)
parser.add_argument('-t','--time',help='Time', required=False,type=float)
parser.add_argument('-it','--it',help='Time index',required=True,type=int)

args = vars(parser.parse_args())

shot = args["shot"]
run  = args["run"]

# User or absolute path name
if args['user_or_path'] != None:
    user = args['user_or_path']
else:
    user = 'public'

# Database name
if args['database'] != None:
    database = args['database']
else:
    database = 'iter'

# Time
if args['time'] != None:
    time = args['time']
    SingleSlice = True
else:
    time = -99.
    SingleSlice = False

it = args["it"]

imas_entry_init = imas.DBEntry(imas.imasdef.MDSPLUS_BACKEND, database, shot, run, user, data_version = '3')
imas_entry_init.open()

idslist = {}

idslist['equilibrium'] = imas_entry_init.get('equilibrium')

# x = idslist['equilibrium'].time_slice[it].profiles_2d[0].grid.dim1
# y = idslist['equilibrium'].time_slice[it].profiles_2d[0].grid.dim2

      
# psi2d = np.transpose(idslist['equilibrium'].time_slice[it].profiles_2d[0].psi)

# psi_bnd = idslist['equilibrium'].time_slice[it].boundary.psi
# psi_sep = idslist['equilibrium'].time_slice[it].boundary_separatrix.psi
# psi_sep2 = idslist['equilibrium'].time_slice[it].boundary_secondary_separatrix.psi
# psi_axis = idslist['equilibrium'].time_slice[it].global_quantities.psi_axis

# psin = np.linspace(1.0, 1.1, num=40)
# contour_levels = sorted(psi_axis + psin * (psi_bnd - psi_axis ))

# # Create a contour plot
# contour_plot = plt.contour(x, y, psi2d, levels=contour_levels, cmap='viridis')
# plt.plot(vv[:,0]*1e-3,vv[:,1]*1e-3)
# plt.plot(fw[:,0],fw[:,1])
# print(idslist['equilibrium'].time[it])
# plt.gca().set_aspect('equal', adjustable='box')
# # Add labels and a colorbar

# plt.xlabel('R [m]')
# plt.ylabel('Z [m]')
# plt.title('Psi')
# plt.colorbar(contour_plot, label='Psi')

# # Show the plot
# plt.show()


for it in range(0,len(idslist['equilibrium'].time)):
    # Assuming x, y, and z are your data arrays
    x_values = idslist['equilibrium'].time_slice[it].profiles_2d[0].grid.dim1
    y_values = idslist['equilibrium'].time_slice[it].profiles_2d[0].grid.dim2
    psi_val  = np.transpose(idslist['equilibrium'].time_slice[it].profiles_2d[0].psi)


    psi_bnd  = idslist['equilibrium'].time_slice[it].boundary.psi
    psi_axis = idslist['equilibrium'].time_slice[it].global_quantities.psi_axis

    interp_function = interp2d(x_values, y_values, psi_val, kind='linear')

    # Take boundary at a given special location (psin=1)
    x1 = 3.9890909;  y1 = 4.63832;
    x1 = 5.0311395;  y1 = 5.13444;
    psi_bnd = interp_function(x1,y1)

    psi_norm = (psi_val-psi_axis)/(psi_bnd-psi_axis)

    # Create a meshgrid from x and y values
    x, y = np.meshgrid(x_values, y_values)

    # Flatten the arrays for writing to VTK
    x_flat = x.flatten()
    y_flat = y.flatten()
    field_values_flat = psi_norm.flatten()

    resolution =len(x_values)

    # Open the VTK file for writing
    with open(f"outputPP2{it}.vtk", "w") as vtk_file:
        # Write the header
        vtk_file.write("# vtk DataFile Version 3.0\n")
        vtk_file.write("Generated by Python\n")
        vtk_file.write("ASCII\n")
        vtk_file.write("DATASET STRUCTURED_POINTS\n")
        vtk_file.write(f"DIMENSIONS {resolution} 1 {resolution}\n")
        vtk_file.write(f"ORIGIN {x_values[0]} 0.0 {y_values[0]}\n")
        vtk_file.write(f"SPACING {(x_values[-1] - x_values[0]) / (resolution - 1)} 0.0 {(y_values[-1] - y_values[0]) / (resolution - 1)}\n")
        vtk_file.write("POINT_DATA " + str(resolution * resolution) + "\n")
        vtk_file.write("SCALARS field_values float 1\n")
        vtk_file.write("LOOKUP_TABLE default\n")

        # Write the field values
        for value in field_values_flat:
            vtk_file.write(f"{value}\n")

    print(f"VTK file {it} has been written successfully.")