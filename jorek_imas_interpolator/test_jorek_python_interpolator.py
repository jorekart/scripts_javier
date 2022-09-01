import numpy as np
import getpass
import argparse

from interp_jorek import *
from imas import imasdef
import imas

import logging
import sys

prec=np.float32

parser = argparse.ArgumentParser(description="Convert IMAS MHD IDS to VTK file",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-s", "--shot", type=int, default=1, help="Shot number")
parser.add_argument("-r", "--run", type=int, default=7, help="Run number")
parser.add_argument("-u", "--user", type=str, default=getpass.getuser(),
                    help="Location of ~$USER/public/imasdb")
parser.add_argument("-d", "--database", type=str, default="smiter", help="Database name under public/imasdb/")
parser.add_argument("-o", "--occurrence", type=int, default=0, help="Occurrence number")
parser.add_argument("-f", "--backend", type=int, default=imasdef.MDSPLUS_BACKEND,
                    help="Database format: 12=MDSPLUS, 13=HDF5")
parser.add_argument("vtkfile", metavar='jorek.vtk', nargs='?', help="Resulting VTK filename", default="jorek_ids.vtu")
parser.add_argument("-p", "--phi", type=list, default=[0, 90], help="Phi coordinate")
parser.add_argument("-n", "--n_plane", type=int, default=3, help="Number of planes")
parser.add_argument("-b", "--bezier", type=bool, default=True, help="Bezier grid")

args = parser.parse_args()

# Open IMAS DB entry
data_entry = imas.DBEntry(args.backend, args.database, args.shot, args.run,
                          user_name=args.user)
status, idx = data_entry.open()
if status != 0:
    logging.info('Creation of data entry FAILED! Exiting.')
    sys.exit(-1)
else:
    logging.info('Creation of data entry Ok!')

ids = data_entry.get("mhd")


def main():

    nlen = 10000
    Rvec = np.linspace(5.0,8.0,nlen)
    Zvec = np.linspace(1,1.0,nlen)
    pvec = np.linspace(0, 0.0,nlen)

#    f = open('grid_test.txt', 'w')
#    for R in Rvec:
#        for Z in Zvec:    
#            val = value_at_RZphi( ids, Rvec, Zvec, 0.0, 'rho' )
#            f.write( str(R) + '  ' + str(Z) + '  ' + str( val ) + '\n' )
#    f.close()

    # Time and grid slice
    i_grid  = 0
    t_slice = 0

    # Get GGD grid
    grid  = ids.grid_ggd[i_grid]

    # Get Bezier coefficients
    var_name = 'mass_density'
    if ( var_name == 'mass_density'):
        val_coeff = ids.ggd.array[t_slice].mass_density.array[0].coefficients
    else:
        print( "Error: variable with name " + var_name + " unkown")

    # Find points in FEM coordinates and get value there
    val = value_at_RZphi( grid, val_coeff, Rvec, Zvec, pvec )

    print('writing')
    f = open('grid_test.txt', 'w')
    for i in range(0,len(Rvec)):
        f.write( str(Rvec[i]) + '  ' + str(Zvec[i]) + '  ' + str( val[i] ) + '\n' )
    f.close()

if __name__ == "__main__":
    main()
