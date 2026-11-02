import h5py
import numpy as np
from basis import basisfunctions1
from constants import *

def store_input_params(group, params):
    """Write the (nested) input parameters as attributes of group, one subgroup per section"""

    for key, value in params.items():
        if isinstance(value, dict):
            store_input_params(group.require_group(key), value)
        elif isinstance(value, (list, tuple)) and any(isinstance(v, str) for v in value):
            group.attrs[key] = ", ".join(str(v) for v in value)    # e.g. reflective: [left, right]
        elif value is None:
            group.attrs[key] = "None"
        else:
            group.attrs[key] = value


def create_solution_file(filename, xvals, input_text, input_params):
    """
    input_text   : content of the input file, stored as is (with comments) in the dataset "input_yaml"
    input_params : parsed input parameters, stored as attributes of the group "input" (one subgroup per section)
    """

    n_points = len(xvals)

    with h5py.File(filename, "w") as f:

        # input parameters
        f.create_dataset("input_yaml", data=input_text, dtype=h5py.string_dtype())
        store_input_params(f.create_group("input"), input_params)

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


def store_balance(filename, balance):
    """Append the particle balance (dict of scalars, see balance.py) to the group "balance", one value per stored time"""

    with h5py.File(filename, "a") as f:

        grp = f.require_group("balance")
        for key, value in balance.items():
            if key not in grp:
                grp.create_dataset(key, shape=(0,), maxshape=(None,), dtype=np.float64)
            it = grp[key].shape[0]
            grp[key].resize(it + 1, axis=0)
            grp[key][it] = value


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
