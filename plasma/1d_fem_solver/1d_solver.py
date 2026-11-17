# This code solves a 1D equation with finite elements
import os
import numpy as np
import yaml

from constants import *
from grid import update_node_values, grid_input_params_class, build_grid_1D
from construct_matrix import construct_matrix
from initial_conditions import initial_conditions
from hdf5_export import solution_real_space, store_solution, create_solution_file


def read_jphi_profile(filename, profile_index=0):
    """Read distance/Jphi/Te/B rows from a selected-field-line fort.444 file."""
    profile_blocks = []
    profile_rows = []
    current_block = None
    has_headers = False

    with open(filename) as profile_file:
        for line_number, line in enumerate(profile_file, start=1):
            stripped = line.strip()
            if stripped.startswith("# fieldline"):
                has_headers = True
                if current_block:
                    profile_blocks.append(current_block)
                current_block = []
                continue
            if not stripped or stripped.startswith("#"):
                continue

            values = stripped.split()
            if len(values) < 3:
                raise ValueError(
                    f"Expected distance, Jphi and Te on line {line_number} of {filename}"
                )
            row = [float(values[0]), float(values[1])]
            if len(values) >= 4:
                row.append(float(values[3]))
            if has_headers:
                current_block.append(row)
            else:
                profile_rows.append(row)

    if current_block:
        profile_blocks.append(current_block)

    if has_headers:
        if profile_index < 0 or profile_index >= len(profile_blocks):
            raise ValueError(
                f"Profile index {profile_index} is out of range for {filename}"
            )
        profile_rows = profile_blocks[profile_index]

    if not profile_rows:
        raise ValueError(f"No profile data found in {filename}")

    profile = np.asarray(profile_rows)
    order = np.argsort(profile[:, 0])
    distance = profile[order, 0]
    # Poincare exports Jphi in MA/m^2; the 1D solver uses A/m^2.
    current = profile[order, 1] * 1.0e6

    magnetic_field = None
    if profile.shape[1] >= 3:
        magnetic_field = profile[order, 2]

    # np.interp requires increasing, unique x values.
    distance, unique_indices = np.unique(distance, return_index=True)
    current = current[unique_indices]
    if magnetic_field is not None:
        magnetic_field = magnetic_field[unique_indices]
    return distance, current, magnetic_field


def main():
    print('This code solves a 1D equation with finite elements')

    # Read input file parameters
    with open("input.yaml", "r") as f:
        input_params = yaml.safe_load(f)

    time_params = input_params["time"]
    print(time_params)

    jphi_filename = input_params["physics"].get("Jphi_profile", "fort.444")
    jphi_profile = None
    if os.path.isfile(jphi_filename):
        jphi_profile_index = int(input_params["physics"].get("Jphi_profile_index", 0))
        jphi_x, jphi_values, B_values = read_jphi_profile(jphi_filename, jphi_profile_index)
        jphi_profile = (jphi_x, jphi_values)
        B_profile = None if B_values is None else (jphi_x, B_values)
        input_params["mesh"]["x_max"] = float(jphi_x[-1])
        print(f"Using Jphi profile {jphi_profile_index} from {jphi_filename} in A/m2")
        if B_profile is not None:
            print(f"Using total magnetic-field profile from {jphi_filename} in T")
        else:
            print(f"No B profile found in {jphi_filename}; continuing without B")
        print(f"Setting mesh x_max from Jphi profile: {input_params['mesh']['x_max']:.8g} m")
    else:
        B_profile = None
        print(f"Jphi profile {jphi_filename} not found; using constant Jphi and YAML x_max")

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

    # Time loop
    for step in range(n_steps):

        print(f"Time step {step+1}/{n_steps}, Time: {t:.4f}")

        # Collect element contributions from equations left hand side and right hand side
        a_mat, rhs_vec = construct_matrix(grid, input_params, jphi_profile, B_profile)

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
