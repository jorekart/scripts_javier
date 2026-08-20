import numpy as np

# Related for FEM
N_ORDER = 3
N_DOF_NODE = int((N_ORDER + 1)/2)
NODES_PER_ELEM = 2

# Related to integration
N_GAUSS = 4
S_GAUSS = np.array([0.0694318442029735, 0.3300094782075720, 0.6699905217924280, 0.9305681557970265])  # Positions of Gaussian points
W_GAUSS = np.array([0.173927422568727,  0.326072577431273, 0.326072577431273,  0.173927422568727])    # Weights of Gaussian points

# Related to physics
N_VAR = 1  # Number of variables in the system of equations

tstep = 0.01