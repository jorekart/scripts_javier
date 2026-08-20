# This code solves a 1D equation with finite elements
import numpy as np
import yaml

from constants import *
from grid import update_node_values, grid_input_params_class, build_grid_1D
from construct_matrix import construct_matrix
from initial_conditions import initial_conditions
from hdf5_export import solution_real_space, store_solution, create_solution_file

def main():
    print('This code solves a 1D equation with finite elements')

    # Read input file parameters
    with open("input.yaml", "r") as f:
        input_params = yaml.safe_load(f)

    time_params = input_params["time"]
    print(time_params)

    # Create grid
    grid_params = grid_input_params_class(**input_params["mesh"])
    grid = build_grid_1D(grid_params)

    # Apply initial conditions to values and deltas at nodes
    grid = initial_conditions(grid, input_params)

    # Start time loop
    t = time_params["t_0"]
    n_steps = time_params["n_steps"]
    n_out = time_params["n_out"]
    tstep = time_params["tstep"]
    filename = "solution.h5"
    n_sub = 4

    # Store initial grid and create file
    xvals, varvals, deltas = solution_real_space(grid, n_sub)
    create_solution_file(filename, xvals)
    store_solution(filename=filename, time=t, varvals=varvals, deltas=deltas)
            

    for step in range(n_steps):

        print(f"Time step {step+1}/{n_steps}, Time: {t:.4f}")

        # Collect element contributions from equations left hand side and right hand side
        a_mat, rhs_vec = construct_matrix(grid, input_params)

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

