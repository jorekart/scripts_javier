# This code solves two 1D diffusion equations with finite elements:
#   dn/dt = d/dx( D_perp   dn/dx )   particle density
#   dT/dt = d/dx( chi_perp dT/dx )   temperature
import sys
import numpy as np
import yaml

from constants import *
from grid import update_node_values, grid_input_params_class, build_grid_1D
from construct_matrix import construct_matrix
from initial_conditions import initial_conditions
from hdf5_export import solution_real_space, store_solution, create_solution_file


def as_list(value, dtype):
    """Accept a scalar, a YAML list or a comma-separated string, and return a list of dtype"""
    if isinstance(value, str):
        value = value.split(",")
    if not isinstance(value, (list, tuple)):
        value = [value]
    return [dtype(v) for v in value]


def main(input_filename="input.yaml", filename="solution.h5"):
    print('This code solves the 1D diffusion equations for n and T with finite elements')

    # Read input file parameters
    with open(input_filename, "r") as f:
        input_params = yaml.safe_load(f)

    time_params = input_params["time"]
    print(time_params)

    # Create grid
    grid_params = grid_input_params_class(**input_params["mesh"])
    grid = build_grid_1D(grid_params)

    # Apply initial conditions to values at nodes
    grid = initial_conditions(grid, input_params)

    # Start time loop
    t = time_params["t_0"]
    n_out = time_params["n_out"]

    # Variable time step: tstep[i] is used during n_steps[i] steps
    tsteps  = as_list(time_params["tstep"], float)
    n_steps_list = as_list(time_params["n_steps"], int)
    if len(tsteps) != len(n_steps_list):
        raise ValueError(f"tstep and n_steps need the same number of values, found {len(tsteps)} and {len(n_steps_list)}")
    step_sizes = np.repeat(tsteps, n_steps_list)    # time step of every step
    n_steps = len(step_sizes)
    n_sub = 4
    ref_values = [input_params["normalization"][name + "_ref"] for name in VAR_NAMES]

    # Store initial grid and create file
    xvals, varvals = solution_real_space(grid, n_sub, ref_values)
    create_solution_file(filename, xvals)
    store_solution(filename=filename, time=t, varvals=varvals)

    # Time loop
    for step in range(n_steps):

        # time step of this step, read by the element contributions
        tstep = step_sizes[step]
        time_params["tstep"] = tstep
        time_params["t"] = t            # time at the start of the step

        print(f"Time step {step+1}/{n_steps}, Time: {t:.4e}")

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
            xvals, varvals = solution_real_space(grid, n_sub, ref_values)
            store_solution(filename=filename, time=t, varvals=varvals)

    print('Time loop finished')



if __name__ == "__main__":
    main(*sys.argv[1:2])
