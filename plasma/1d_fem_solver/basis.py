import numpy as np
from constants import *


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