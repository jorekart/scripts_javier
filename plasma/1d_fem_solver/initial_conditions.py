import numpy as np
from constants import *

def initial_conditions(grid):
    for node in grid.nodes:
        x = node.coord[0]
        node.values[0, 0] = 1.0
        node.values[0, 1] = 0.0 # Example initial condition for the second DOF

    return grid