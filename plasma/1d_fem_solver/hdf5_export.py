import h5py
import numpy as np
from basis import basisfunctions1
from constants import *

def create_solution_file(filename, xvals):

    n_points = len(xvals)

    with h5py.File(filename, "w") as f:

        # fixed mesh coordinates
        f.create_dataset("xvals", data=xvals)

        # time axis
        f.create_dataset(
            "time",
            shape=(0,),
            maxshape=(None,),
            dtype=np.float64
        )

        # solution history
        f.create_dataset(
            "varvals",
            shape=(0, n_points),
            maxshape=(None, n_points),
            dtype=np.float64
        )

        f.create_dataset(
            "deltas",
            shape=(0, n_points),
            maxshape=(None, n_points),
            dtype=np.float64
        )


def store_solution(filename, time, varvals, deltas):

    with h5py.File(filename, "a") as f:

        it = f["time"].shape[0]

        #
        # append time
        #
        f["time"].resize(it + 1, axis=0)
        f["time"][it] = time

        #
        # append solution
        #
        f["varvals"].resize(it + 1, axis=0)
        f["varvals"][it, :] = varvals

        #
        # append delta
        #
        f["deltas"].resize(it + 1, axis=0)
        f["deltas"][it, :] = deltas


def solution_real_space(grid, n_sub):

    svec = np.linspace(0,1,num=n_sub)
    xvals = []
    deltas = []
    varvals = []

    for element in grid.elements:
        sizes = element.sizes

        for s in svec:
            H, H_s, H_ss = basisfunctions1(s)
      
            x = 0.0
            var = 0.0
            delt = 0.0
            for i in range(NODES_PER_ELEM):
                inode = element.vertex_glob_index[i]
                node = grid.nodes[inode]
                for j in range(N_DOF_NODE):
                        x    += H[i,j] * sizes[i,j] * node.coord[j]
                        var  += H[i,j] * sizes[i,j] * node.values[0,j]
                        delt += H[i,j] * sizes[i,j] * node.deltas[0,j]

            xvals.append(x)
            varvals.append(var)
            deltas.append(delt)

    return xvals, varvals, deltas