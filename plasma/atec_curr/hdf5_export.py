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

        # solution history of each variable (n, T)
        for name in VAR_NAMES:
            f.create_dataset(
                name,
                shape=(0, n_points),
                maxshape=(None, n_points),
                dtype=np.float64
            )


def store_solution(filename, time, varvals):
    """varvals : array of shape (N_VAR, n_points)"""

    with h5py.File(filename, "a") as f:

        it = f["time"].shape[0]

        f["time"].resize(it + 1, axis=0)
        f["time"][it] = time

        for i_var, name in enumerate(VAR_NAMES):
            f[name].resize(it + 1, axis=0)
            f[name][it, :] = varvals[i_var]


def solution_real_space(grid, n_sub, ref_values):
    """ref_values : reference values of each variable, to return physical (not normalized) values"""

    svec = np.linspace(0,1,num=n_sub)
    xvals = []
    varvals = [[] for _ in range(N_VAR)]

    for element in grid.elements:
        sizes = element.sizes

        for s in svec:
            H, H_s, H_ss = basisfunctions1(s)

            x = 0.0
            var = np.zeros(N_VAR)
            for i in range(NODES_PER_ELEM):
                inode = element.vertex_glob_index[i]
                node = grid.nodes[inode]
                for j in range(N_DOF_NODE):
                    x   += H[i,j] * sizes[i,j] * node.coord[j]
                    var += H[i,j] * sizes[i,j] * node.values[:,j]

            xvals.append(x)
            for i_var in range(N_VAR):
                varvals[i_var].append(var[i_var])

    return np.array(xvals), np.array(varvals) * np.array(ref_values)[:, None]
