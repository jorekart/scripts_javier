# This code solves a 1D equation with finite elements
import numpy as np
from dataclasses import dataclass, field
import matplotlib.pyplot as plt

# define basis functions for 1D Bezier elements
N_ORDER = 3
N_DOF_NODE = int((N_ORDER + 1)/2)
NODES_PER_ELEM = 2
N_GAUSS = 4
N_VAR = 1  # Number of variables in the system of equations

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
    grid = grid_1D_class(n_nodes=n_nodes, n_elements=n_elements)

    x_ac_unit = meshac2(n_nodes, x_acc1, x_acc2, x_sig1, x_sig2, 0.5, 1.0)
    x_coords = x_min + (x_max - x_min) * x_ac_unit

    print(f"x_coords: {x_coords}")
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

def main():
    print('This code solves a 1D equation with finite elements')

    # Hardcoded parameters/ read input file?

    # Initialize basis functions and gauss points

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

    # Collect element contributions from equations left hand side and right hand side

    # Solve the system of equations

    # Print





if __name__ == "__main__":
    main()


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