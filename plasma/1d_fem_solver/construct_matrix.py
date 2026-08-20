import numpy as np
from constants import *
from basis import basisfunctions_gauss
from element_contributions import get_element_contributions, index_local


def index_global(i_node, i_dof):

    return N_DOF_NODE * i_node + i_dof


def construct_matrix(grid):

    n_dofs_tot = grid.n_dofs_tot
    a_mat = np.zeros((n_dofs_tot, n_dofs_tot))
    rhs_vec = np.zeros(n_dofs_tot)

    # Initialize basis functions and gauss points
    H_gauss, H_s_gauss, H_ss_gauss = basisfunctions_gauss()

    print("Matrix construction and RHS vector assembly:")

    for element in grid.elements:
        
        inode1 = element.vertex_glob_index[0]
        inode2 = element.vertex_glob_index[1]

        nodes = [grid.nodes[inode1], grid.nodes[inode2]]

        a_loc, rhs_loc = get_element_contributions(element, nodes, H_gauss, H_s_gauss, H_ss_gauss)

        for i_v in range(NODES_PER_ELEM):
            for i_dof in range(N_DOF_NODE):
                i_node  = element.vertex_glob_index[i_v]
                index_i = index_global(i_node, i_dof)
                index_i_loc = index_local(i_v, i_dof)

                rhs_vec[index_i] += rhs_loc[index_i_loc]

                for j_v in range(NODES_PER_ELEM):
                    for j_dof in range(N_DOF_NODE):
                        j_node  = element.vertex_glob_index[j_v]
                        index_j = index_global(j_node, j_dof)
                        index_j_loc = index_local(j_v, j_dof)

                        a_mat[index_i, index_j] += a_loc[index_i_loc, index_j_loc]

    return a_mat, rhs_vec
