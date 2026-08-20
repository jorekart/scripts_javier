import numpy as np
from constants import *
from basis import basisfunctions1


def f1(x,x0,w):
    f = np.exp(-((x-x0)/w)**2)
    f_x = - f * ((x-x0)/w) * 2.0 / w
    return f, f_x

def initial_conditions(grid, input_params):

    i_var = 0

    for inode, node in enumerate(grid.nodes):
        x = node.coord[0] 
        x_s = node.coord[1]
        f, f_x = f1(x,2.5,0.8)
        node.values[i_var, 0] = f
        node.values[i_var, 1] = f_x * x_s 

    return grid