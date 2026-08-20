# This code solves a 1D equation with finite elements
import numpy as np
from dataclasses import dataclass, field
import matplotlib.pyplot as plt
import h5py

# define basis functions for 1D Bezier elements
N_ORDER = 3
N_DOF_NODE = int((N_ORDER + 1)/2)
NODES_PER_ELEM = 2
N_GAUSS = 4
N_VAR = 1  # Number of variables in the system of equations

S_GAUSS = np.array([0.0694318442029735, 0.3300094782075720, 0.6699905217924280, 0.9305681557970265])  # Positions of Gaussian points
W_GAUSS = np.array([0.173927422568727,  0.326072577431273, 0.326072577431273,  0.173927422568727])    # Weights of Gaussian points

tstep = 0.01



def basisfunctions1(s):
    """
    1D mixed Bezier/Cubic finite element basis functions.

    Parameters
    ----------
    s : float
        Local coordinate in [0, 1].

    Returns
    -------
    H : ndarray, shape (NODES_PER_ELEM, N_DOF_NODE)
        Basis functions.
    H_s : ndarray, shape (NODES_PER_ELEM, N_DOF_NODE)
        First derivatives.
    H_ss : ndarray, shape (NODES_PER_ELEM, N_DOF_NODE)
        Second derivatives.
    """

    H = np.zeros((NODES_PER_ELEM, N_DOF_NODE))
    H_s = np.zeros((NODES_PER_ELEM, N_DOF_NODE))
    H_ss = np.zeros((NODES_PER_ELEM, N_DOF_NODE))

    # ---------------------------------------------------------- vertex (1)
    H[0, 0]    = (s - 1.0)**2 * (1.0 + 2.0*s)
    H_s[0, 0]  = 6.0 * (s - 1.0) * s
    H_ss[0, 0] = 6.0 * (2.0*s - 1.0)

    H[0, 1]    = 3.0 * (s - 1.0)**2 * s
    H_s[0, 1]  = 3.0 * (s - 1.0) * (-1.0 + 3.0*s)
    H_ss[0, 1] = 6.0 * (-2.0 + 3.0*s)

    # ---------------------------------------------------------- vertex (2)
    H[1, 0]    = -s**2 * (-3.0 + 2.0*s)
    H_s[1, 0]  = -6.0 * (s - 1.0) * s
    H_ss[1, 0] = -6.0 * (2.0*s - 1.0)

    H[1, 1]    = -3.0 * (s - 1.0) * s**2
    H_s[1, 1]  = -3.0 * s * (-2.0 + 3.0*s)
    H_ss[1, 1] = -6.0 * (-1.0 + 3.0*s)

    return H, H_s, H_ss

# basis at gaussian points
def basisfunctions_gauss():

    H_gauss = np.zeros((N_GAUSS, NODES_PER_ELEM, N_DOF_NODE))
    H_s_gauss = np.zeros((N_GAUSS, NODES_PER_ELEM, N_DOF_NODE))
    H_ss_gauss = np.zeros((N_GAUSS, NODES_PER_ELEM, N_DOF_NODE))

    for i in range(N_GAUSS):
        s = S_GAUSS[i]
        H_gauss[i], H_s_gauss[i], H_ss_gauss[i] = basisfunctions1(s)

    return H_gauss, H_s_gauss, H_ss_gauss

def meshac2(
    nr: int,
    xr1: float,
    xr2: float,
    sig1: float,
    sig2: float,
    bgf: float,
    fact: float,
    sg0: float = 0.0,
    nmax: int = 10001,
):
    """
    Construct a non-uniform 1D mesh using a Gaussian-based
    mesh density function.

    Parameters
    ----------
    nr : int
        Number of mesh points.

    xr1, xr2 : float
        Locations of the Gaussian clustering regions.

    sig1, sig2 : float
        Widths of the Gaussian clustering regions.

    bgf : float
        Background density level.

    fact : float
        Amplitude factor of the Gaussian contribution.

    sg0 : float, optional
        Starting coordinate of the mesh.
        Usually 0.0.

    nmax : int, optional
        Number of integration points used internally.

    Returns
    -------
    sg : ndarray, shape (nr,)
        Non-uniform mesh coordinates in [0,1].

    Notes
    -----
    The algorithm:

    1. Samples a density function rho(s)=FGAUS(...)
       on a fine background grid.
    2. Computes the cumulative integral using the
       trapezoidal rule.
    3. Normalizes the cumulative integral to [0,1].
    4. Inverts the cumulative mapping to generate
       a stretched mesh.
    """

    if nr <= 0:
        return np.array([])

    sg0 = np.clip(sg0, 0.0, 1.0)

    # Fine grid for numerical integration
    s1 = np.linspace(sg0, 1.0, nmax)

    # Evaluate density function
    rho = np.array([
        fgaus(s, bgf, xr1, xr2, sig1, sig2, fact)[0]
        for s in s1
    ])


    # Cumulative integral (trapezoidal rule)
    fsum = np.zeros(nmax)

    for i in range(1, nmax):
        ds = s1[i] - s1[i - 1]
        fsum[i] = (
            fsum[i - 1]
            + 0.5 * (rho[i - 1] + rho[i]) * ds
        )

    # Normalize cumulative integral
    fsum /= fsum[-1]

    # Uniform coordinates in cumulative space
    fi = np.linspace(0.0, 1.0, nr)

    # Inverse mapping
    sg = np.interp(fi, fsum, s1)

    sg[0] = 0.0
    sg[-1] = 1.0

    return sg


def fgaus(
    zs: float,
    bgf: float,
    xr1: float,
    xr2: float,
    sig1: float,
    sig2: float,
    fact: float,
):
    """
    Double-Gaussian mesh density function.

    Parameters
    ----------
    zs : float
        Evaluation coordinate.

    bgf : float
        Background factor (0 <= bgf <= 1).

    xr1, xr2 : float
        Centers of the two Gaussian peaks.

    sig1, sig2 : float
        Standard deviations of the two Gaussians.

    fact : float
        Weighting factor of the second Gaussian.

    Returns
    -------
    fgaus : float
        Mesh density function.

    dfgaus : float
        Derivative of the mesh density function.
    """

    znorm1 = 0.39894 / sig1
    znorm2 = 0.39894 / sig2

    zex1 = -0.5 * (zs - xr1)**2 / sig1**2
    zex2 = -0.5 * (zs - xr2)**2 / sig2**2

    dex1 = -(zs - xr1) / sig1**2
    dex2 = -(zs - xr2) / sig2**2

    f1 = znorm1 * np.exp(zex1)
    f2 = znorm2 * np.exp(zex2)

    df1 = znorm1 * dex1 * np.exp(zex1)
    df2 = znorm2 * dex2 * np.exp(zex2)

    fgaus = bgf + (1.0 - bgf) * (f1 + fact * f2) / fact

    dfgaus = (1.0 - bgf) * (df1 + fact * df2) / fact

    return fgaus, dfgaus

@dataclass
class element_1D_class:
    vertex_glob_index: np.ndarray   #  contains global indices of the nodes of the element
    sizes: np.ndarray             #  contains sizes of the element in each direction

@dataclass
class node_class:
    coord: np.ndarray  # contains coordinates of the node and their dofs
    values: np.ndarray  # contains values and dofs at the node
    deltas: np.ndarray  # contains deltas of the node and their dofs (changes in values)

@dataclass
class grid_1D_class:
    n_nodes: int       # number of nodes 
    n_elements: int    # number of elements
    nodes: list = field(default_factory=list) # list of node_class instances
    elements: list = field(default_factory=list) # list of element_1D_class instances
    n_dofs_tot: int = 0  # total number of degrees of freedom in the grid

@dataclass
class grid_input_params_class:
    x_min: float  # minimum x-coordinate of the grid
    x_max: float  # maximum x-coordinate of the grid
    n_elements: int  # number of elements in the grid
    x_acc1: float = 0.05  # accumulation factor for the first Gaussian
    x_acc2: float = 0.95   # accumulation factor for the second Gaussian
    x_sig1: float = 0.1     # Standard deviation for the first Gaussian
    x_sig2: float = 0.1     # Standard deviation for the second Gaussian


# Build 1D grid from input parameters
def build_grid_1D(grid_input_params):
    n_elements = grid_input_params.n_elements
    x_min = grid_input_params.x_min
    x_max = grid_input_params.x_max
    x_acc1 = grid_input_params.x_acc1
    x_acc2 = grid_input_params.x_acc2
    x_sig1 = grid_input_params.x_sig1
    x_sig2 = grid_input_params.x_sig2

    n_nodes = n_elements + 1
    n_dofs_tot = n_nodes * N_DOF_NODE * N_VAR
    grid = grid_1D_class(n_nodes=n_nodes, n_elements=n_elements, n_dofs_tot=n_dofs_tot)

    x_ac_unit = meshac2(n_nodes, x_acc1, x_acc2, x_sig1, x_sig2, 0.5, 1.0)
    x_coords = x_min + (x_max - x_min) * x_ac_unit

    
    for i in range(n_nodes):
        node = node_class(coord=np.zeros(N_DOF_NODE),
                        values=np.zeros((N_VAR, N_DOF_NODE)),
                        deltas=np.zeros((N_VAR, N_DOF_NODE)))

        element = element_1D_class(vertex_glob_index=np.zeros(NODES_PER_ELEM, dtype=int),
                                sizes=np.ones((NODES_PER_ELEM, N_DOF_NODE)))
        inode1 = i
        inode2 = i + 1 

        node.coord[0] = x_coords[inode1]
        node.coord[1] = 1.0

        if i < n_elements:
            x2 = x_coords[inode2]
            x1 = x_coords[inode1]
            h_elem = x2 - x1
            if (h_elem <= 0.0):
                raise ValueError(f"Element size must be positive. Found h_elem={h_elem} for element {i}.")

            element.vertex_glob_index[0] = inode1
            element.vertex_glob_index[1] = inode2 
            element.sizes[0, 1] = h_elem / float(N_ORDER)      # the size of the second node points at -x 
            element.sizes[1, 1] = -h_elem / float(N_ORDER)      # the size of the second node points at -x 
            
            grid.elements.append(element)
        
        grid.nodes.append(node)
        
    return grid

def index_local(i_v, i_dof):

    return N_DOF_NODE*i_v +  i_dof  

def index_global(i_node, i_dof):

    return N_DOF_NODE * i_node + i_dof


def get_element_contributions(element, nodes,H_gauss, H_s_gauss, H_ss_gauss):

    a_loc   = np.zeros((NODES_PER_ELEM * N_DOF_NODE, NODES_PER_ELEM * N_DOF_NODE))
    rhs_loc = np.zeros(NODES_PER_ELEM * N_DOF_NODE)    

    # This should be read from outside
    theta = 0.5
    zeta  = 0.0
    i_var = 0  # Assuming a single variable for simplicity

    # Get quantities at gaussian points
    x0_s  = np.zeros(N_GAUSS)
    f0    = np.zeros(N_GAUSS)
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

    
    for igauss in range(N_GAUSS):
        x_jac = x0_s[igauss]
        var0  = f0[igauss]
        wg    = W_GAUSS[igauss]
        H = H_gauss[igauss]
        H_s = H_s_gauss[igauss]
        H_ss = H_ss_gauss[igauss]

        for i_v in range(NODES_PER_ELEM):
            node = nodes[i_v]
            for i_dof in range(N_DOF_NODE):
                index_ij = index_local(i_v, i_dof)

                v_test =  H[i_v,i_dof] * sizes[i_v,i_dof] 
                rhs_var = - v_test * var0 

                rhs_loc[index_ij] +=  rhs_var * wg * x_jac * tstep

                for j_v in range(NODES_PER_ELEM):
                    for j_dof in range(N_DOF_NODE):

                        index_kl = index_local(j_v, j_dof)

                        bas = H[j_v,j_dof] * sizes[j_v,j_dof] 

                        amat_var = (+ v_test * (1+zeta)  * bas 
                                    + v_test             * bas * theta * tstep )

                        a_loc[index_ij,index_kl] +=  amat_var * wg * x_jac

    return a_loc, rhs_loc


def construct_matrix(grid, H_gauss, H_s_gauss, H_ss_gauss):

    n_dofs_tot = grid.n_dofs_tot
    a_mat = np.zeros((n_dofs_tot, n_dofs_tot))
    rhs_vec = np.zeros(n_dofs_tot)

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


def initial_conditions(grid):
    for node in grid.nodes:
        x = node.coord[0]
        node.values[0, 0] = 1.0
        node.values[0, 1] = 0.0 # Example initial condition for the second DOF

    return grid

def update_node_values(grid, sol):
    for i_node, node in enumerate(grid.nodes):
        for i_dof in range(N_DOF_NODE):
            index = index_global(i_node, i_dof)
            node.values[0, i_dof] += sol[index]
            node.deltas[0, i_dof] = sol[index]
    return grid


import h5py
import numpy as np

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
    


def main():
    print('This code solves a 1D equation with finite elements')

    # Hardcoded parameters/ read input file?

    # Initialize basis functions and gauss points
    H_gauss, H_s_gauss, H_ss_gauss = basisfunctions_gauss()

    # Build and initialize the grid
    grid_input_params = grid_input_params_class(x_min=0.0, 
                                                x_max=5.0, 
                                                x_acc1=0.0,
                                                x_acc2=1.0,
                                                x_sig1=0.3,
                                                x_sig2=0.3,
                                                n_elements=100)

    grid = build_grid_1D(grid_input_params)

    # Apply initial conditions to values and deltas at nodes
    grid = initial_conditions(grid)

    # Start time loop
    t = 0.0
    n_steps = 101
    n_out = 20
    filename = "solution.h5"
    n_sub = 4

    # Store initial grid and create file
    xvals, varvals, deltas = solution_real_space(grid, n_sub)
    create_solution_file(filename, xvals)
    store_solution(filename=filename, time=t, varvals=varvals, deltas=deltas)
            

    for step in range(n_steps):

        print(f"Time step {step+1}/{n_steps}, Time: {t:.4f}")

        # Collect element contributions from equations left hand side and right hand side
        a_mat, rhs_vec = construct_matrix(grid, H_gauss, H_s_gauss, H_ss_gauss)

        # Solve the system of equations
        sol = np.linalg.solve(a_mat, rhs_vec)

        # Update node values with the solution
        grid = update_node_values(grid, sol)

        # advance time
        t += tstep

        # store solution
        if (np.mod(step,n_out)==0):
            print(f"Storing solution for step = {step+1}")
            xvals, varvals, deltas = solution_real_space(grid, n_sub)
            store_solution(filename=filename, time=t, varvals=varvals, deltas=deltas)

    print('Time loop finished')



if __name__ == "__main__":
    main()


    # CHECK GRID!S
    # svec = np.linspace(0,1,num=4)
    # xvals = []
    # for element in grid.elements:
    #     sizes = element.sizes

    #     for s in svec:
    #         H, H_s, H_ss = basisfunctions1(s)
      
    #         x = 0.0
    #         for i in range(NODES_PER_ELEM):
    #             inode = element.vertex_glob_index[i]
    #             node = grid.nodes[inode]
    #             for j in range(N_DOF_NODE):
    #                     x = x + H[i,j] * sizes[i,j] * node.coord[j]

    #         xvals.append(x)

    # x_old = -1e10
    # for x in xvals:
    #     if (x < x_old):
    #         raise ValueError(f"x values are not monotonically increasing: {x} < {x_old}")
    #     x_old = x   
    # print(xvals)
    # plt.plot(xvals, xvals, 'o-')
    # plt.xlabel('s')
    # plt.ylabel('x')
    # plt.title('1D Grid')
    # plt.show()