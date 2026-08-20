import numpy as np
from dataclasses import dataclass, field
from constants import *
from construct_matrix import index_global

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


def update_node_values(grid, sol):
    for i_node, node in enumerate(grid.nodes):
        for i_dof in range(N_DOF_NODE):
            index = index_global(i_node, i_dof)
            node.values[0, i_dof] += sol[index]
            node.deltas[0, i_dof] = sol[index]
    return grid



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