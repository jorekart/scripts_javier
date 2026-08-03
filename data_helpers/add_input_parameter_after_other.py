import os

def add_input_param(file_path):
    try:
        with open(file_path, 'r') as file:
            lines = file.readlines()

        with open(file_path, 'w') as file:
            added_input_param = False
            for line in lines:
                file.write(line)
                if "F0" in line and not added_input_param:        # Put here the parameter for reference
                    file.write(" eta_coul_log_dep = .false.\n")   # Put here your input parameter
                    added_input_param = True

        if added_input_param:
            print(f"Added string after 'F0' in {file_path}") 
        else:
            print(f"String 'F0' not found in {file_path}")

    except IOError:
        print(f"Error processing file: {file_path}")


def add_input_param_to_input_files(starting_directory):
    for root, dirs, files in os.walk(starting_directory):
        # Filter directories to exclude those that don't contain "600" or "750" in their names
        dirs[:] = [d for d in dirs if "600" in d or "750" in d]

        for file in files:
            if file == "input":
                file_path = os.path.join(root, file)
                add_input_param(file_path)

# Directory where the search for 'input' files starts
starting_directory = "."

add_input_param_to_input_files(starting_directory)