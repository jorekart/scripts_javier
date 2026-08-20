# This code solves a 1D equation with finite elements
import numpy as np

from constants import *
from grid import update_node_values, grid_input_params_class, build_grid_1D
from construct_matrix import construct_matrix
from initial_conditions import initial_conditions
from hdf5_export import solution_real_space, store_solution, create_solution_file

def main():
    print('This code solves a 1D equation with finite elements')

    # Hardcoded parameters/ read input file?


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
        a_mat, rhs_vec = construct_matrix(grid)

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

