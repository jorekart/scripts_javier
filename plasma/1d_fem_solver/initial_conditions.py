import numpy as np
from constants import *
from basis import basisfunctions1


def f1(x,x0,w):
    f = np.exp(-((x-x0)/w)**2)
    f_x = - f * ((x-x0)/w) * 2.0 / w
    return f, f_x

def initial_conditions(grid, input_params):

    i_var = 0

    Te0 = input_params["physics"]["Te0"]
    Te0_bnd = input_params["physics"]["Te0_bnd"]
    x_min = input_params["mesh"]["x_min"]
    x_max = input_params["mesh"]["x_max"]
    L0  = x_max - x_min
    Lmid = L0*0.5
    width = Lmid*0.3

    for inode, node in enumerate(grid.nodes):
        x = node.coord[0] 
        x_s = node.coord[1]
        f, f_x = f1(x,Lmid,width)
        node.values[i_var, 0] = f * (Te0-Te0_bnd) + Te0_bnd 
        node.values[i_var, 1] = f_x * x_s * (Te0-Te0_bnd)

    return grid