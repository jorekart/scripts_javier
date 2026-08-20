import numpy as np
from constants import *


def index_local(i_v, i_dof):

    return N_DOF_NODE*i_v +  i_dof  


def get_element_contributions(element, nodes,
                              H_gauss, H_s_gauss, H_ss_gauss, input_params):

    a_loc   = np.zeros((NODES_PER_ELEM * N_DOF_NODE, NODES_PER_ELEM * N_DOF_NODE))
    rhs_loc = np.zeros(NODES_PER_ELEM * N_DOF_NODE)    

    theta = input_params["time"]["theta"]
    zeta  = input_params["time"]["zeta"]
    tstep = input_params["time"]["tstep"]

    i_var = 0  # Assuming a single variable for simplicity

    # Get quantities at gaussian points
    x0_s  = np.zeros(N_GAUSS)
    f0    = np.zeros(N_GAUSS)
    f0_s  = np.zeros(N_GAUSS)
    x0    = np.zeros(N_GAUSS)
    sizes = element.sizes
    rhs_loc = np.zeros(NODES_PER_ELEM * N_DOF_NODE)
    a_loc = np.zeros((NODES_PER_ELEM * N_DOF_NODE, NODES_PER_ELEM * N_DOF_NODE))

    for igauss in range(N_GAUSS):
        H = H_gauss[igauss]
        H_s = H_s_gauss[igauss]
        H_ss = H_ss_gauss[igauss]

        for i_v in range(NODES_PER_ELEM):
            node = nodes[i_v]
            for i_dof in range(N_DOF_NODE):
                x0[igauss]    +=  H[i_v,i_dof]   * sizes[i_v,i_dof] * node.coord[i_dof]
                x0_s[igauss]  +=  H_s[i_v,i_dof] * sizes[i_v,i_dof] * node.coord[i_dof]
                f0[igauss]    +=  H[i_v,i_dof]   * sizes[i_v,i_dof] * node.values[i_var,i_dof]
                f0_s[igauss]  +=  H_s[i_v,i_dof] * sizes[i_v,i_dof] * node.values[i_var,i_dof]

    
    for igauss in range(N_GAUSS):
        x_jac = x0_s[igauss]
        var0  = f0[igauss]
        var0_x = f0_s[igauss] / x_jac
        wg    = W_GAUSS[igauss]
        H = H_gauss[igauss]
        H_s = H_s_gauss[igauss]
        H_ss = H_ss_gauss[igauss]

        for i_v in range(NODES_PER_ELEM):
            node = nodes[i_v]
            for i_dof in range(N_DOF_NODE):
                index_ij = index_local(i_v, i_dof)

                v_test = H[i_v,i_dof] * sizes[i_v,i_dof]
                v_x    = H_s[i_v,i_dof] * sizes[i_v,i_dof] / x_jac

                rhs_var  = - v_x * var0_x
                #rhs_var = - v_test * var0 

                rhs_loc[index_ij] +=  rhs_var * wg * x_jac * tstep

                for j_v in range(NODES_PER_ELEM):
                    for j_dof in range(N_DOF_NODE):

                        index_kl = index_local(j_v, j_dof)

                        bas = H[j_v,j_dof] * sizes[j_v,j_dof] 
                        bas_x = H_s[j_v,j_dof] * sizes[j_v,j_dof] / x_jac

                        amat_var = (+ v_test * (1+zeta)  * bas 
                                    #+ v_test             * bas * theta * tstep )
                                    + v_x * bas_x * theta * tstep )

                        a_loc[index_ij,index_kl] +=  amat_var * wg * x_jac

    return a_loc, rhs_loc
