import numpy as np
from constants import *
from basis import basisfunctions_gauss
from element_contributions import get_element_contributions


def index_global(i_node, i_var, i_dof):

    return (i_node * N_VAR + i_var) * N_DOF_NODE + i_dof

def construct_matrix(grid, input_params):

    n_dofs_tot = grid.n_dofs_tot
    a_mat = np.zeros((n_dofs_tot, n_dofs_tot))
    rhs_vec = np.zeros(n_dofs_tot)

    # Initialize basis functions and gauss points
    H_gauss, H_s_gauss, H_ss_gauss = basisfunctions_gauss()

    for element in grid.elements:

        nodes = [grid.nodes[i] for i in element.vertex_glob_index]

        a_loc, rhs_loc = get_element_contributions(element, nodes, H_gauss, H_s_gauss, H_ss_gauss, input_params)

        # global indices in the local ordering (variable, vertex, dof)
        index_glob = np.array([index_global(element.vertex_glob_index[i_v], i_var, i_dof)
                               for i_var in range(N_VAR)
                               for i_v in range(NODES_PER_ELEM)
                               for i_dof in range(N_DOF_NODE)])

        rhs_vec[index_glob] += rhs_loc
        a_mat[np.ix_(index_glob, index_glob)] += a_loc

    # Apply Dirichlet boundary conditions (left and right ends) for n, T and nn.
    # The unknown is the change of the value, so the rhs is (target - current value).
    # Reflective (zero flux) conditions are the natural condition of the weak form, so nothing is
    # imposed for them: the sides listed in "reflective", and nn at both ends if it has no entry.
    bnd_conditions = input_params["bnd_conditions"]
    for i_var, name in ((I_N, "n"), (I_T, "T"), (I_NN, "nn")):
        if name not in bnd_conditions:
            continue
        ref = input_params["normalization"][REF_NAMES[name]]   # variables are normalized
        reflective = bnd_conditions[name].get("reflective", [])
        for i_node, side in ((0, "left"), (grid.n_nodes - 1, "right")):
            if side in reflective:
                continue
            index_i = index_global(i_node, i_var, 0)
            a_mat[index_i,:] = 0.0
            a_mat[index_i,index_i] = 1.0
            rhs_vec[index_i] = bnd_conditions[name][side] / ref - grid.nodes[i_node].values[i_var, 0]

    return a_mat, rhs_vec
