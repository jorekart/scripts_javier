import numpy as np

# Related for FEM
N_ORDER = 3
N_DOF_NODE = int((N_ORDER + 1)/2)
NODES_PER_ELEM = 2

# Related to integration
N_GAUSS = 4
S_GAUSS = np.array([0.0694318442029735, 0.3300094782075720, 0.6699905217924280, 0.9305681557970265])  # Positions of Gaussian points
W_GAUSS = np.array([0.173927422568727,  0.326072577431273, 0.326072577431273,  0.173927422568727])    # Weights of Gaussian points

# Variables of the system of equations
VAR_NAMES = ["n", "T", "nn"]   # particle density, temperature, neutral density
N_VAR = len(VAR_NAMES)
I_N  = 0
I_T  = 1
I_NN = 2

# Reference value used to normalize each variable (the neutral density is normalized with n_ref)
REF_NAMES = {"n": "n_ref", "T": "T_ref", "nn": "n_ref"}

# Related to physics
E_CHARGE = 1.60217662e-19
E_MASS   = 9.1093837e-31
ATOMIC_MASS = 2.01  # D as reference
ION_MASS = 1.67e-27 * ATOMIC_MASS
KB       = 1.380649 * 1e-23
