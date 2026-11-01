import numpy as np
from constants import *


def initial_conditions(grid, input_params):
    """
    Parabolic profile from the left to the right boundary value, for each variable:
    f(x) = right + (left - right) * (1 - s^2),  s = (x - x_min) / L
    i.e. zero slope at the left end and the boundary values at both ends.
    """

    x_min = input_params["mesh"]["x_min"]
    x_max = input_params["mesh"]["x_max"]
    L0  = x_max - x_min

    for i_var, name in enumerate(VAR_NAMES):
        left  = input_params["bnd_conditions"][name]["left"]
        right = input_params["bnd_conditions"][name]["right"]

        # normalized variables
        ref   = input_params["normalization"][name + "_ref"]
        left  = left / ref
        right = right / ref

        for node in grid.nodes:
            x = node.coord[0]
            x_s = node.coord[1]
            s = (x - x_min) / L0
            node.values[i_var, 0] = right + (left - right) * (1.0 - s**2)
            node.values[i_var, 1] = - (left - right) * 2.0 * s / L0 * x_s

    return grid
