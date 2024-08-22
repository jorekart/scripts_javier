# This script goes over a selected set of simulation folders (listed on a text file)
# and exports it to an IMAS DB (specified in imas.nml)
# At the moment, the code uses yaml files to describe simulation phases in which the 
# JOREK input parameters vary
import os
import shutil
import h5py
import yaml
import fileinput

# A class to define a simulation phase
class sim_phase:
    def __init__(self, name, input_file, t_start, t_end, h5_start, h5_end, n_slices):
        self.name       = name
        self.input_file = input_file
        self.t_start    = float(t_start)
        self.t_end      = float(t_end)
        self.h5_start   = int(h5_start)
        self.h5_end     = int(h5_end)
        self.n_slices   = int(n_slices)

    @classmethod
    def from_dict(cls, data):
        return cls(data['name'], data['input_file'],data['t_start'],data['t_end'],data['h5_start'],data['h5_end'],data['n_slices'])
    
# A routine to copy h5 files with conditions
def copy_hdf5_files(source_folder, destination_folder, dt, min_time, max_time):
    if not os.path.exists(destination_folder):
        os.makedirs(destination_folder)

    # List all files in the source folder
    files = [f for f in os.listdir(source_folder) if f.endswith('.h5') and f != 'jorek_restart.h5']

    # Sort files based on numerical part of the filename
    files.sort(key=lambda x: int(''.join(filter(str.isdigit, x))))

    selected_files = []
    prev_time = None

    for file in files:
        file_path = os.path.join(source_folder, file)

        with h5py.File(file_path, 'r') as f:
            time = f['t_now'][()]  # Get the 'time' parameter from the HDF5 file

            if (prev_time is None or ((time - prev_time) >= dt)) and (time >= min_time) and (time <= max_time):
                selected_files.append(file)
                prev_time = time

    # Copy selected files to the destination folder
    for file in selected_files:
        shutil.copy2(os.path.join(source_folder, file), os.path.join(destination_folder, file))

    print(f"Selected and copied {len(selected_files)} files.")

# A routine to replace a parameter in a file
def replace_parameter(file_name,var,value):
    with fileinput.FileInput(file_name, inplace=True, backup='.bak') as file:
        for line in file:                         # go over lines
            if var in line:                       # check if the variable is in the line
                split_line = line.split("=",1)    # divide in 2 the line by "=" delimiter
                print(line.replace(split_line[1]," " + value + "\n"), end='')  # replace value on the RHS
            else: 
                print(line.replace("1e26275462165", value), end='')  # absurd string to search, not elegant, but needed

################# INPUTS ################################################################################################               
skip_export      = True                   # If True, don't write to any database (just to print shots and folders)
file_folder_list = 'folder_list.txt'      # Path of the text file listing the simulation folders
phase_file       = "phase_file.yaml"      # Name of the yaml file describing the simulation phases
imas_file        = "imas.nml"             # Input file for jorek2_IDS
job_IDS          = "job_IDS"              # empty template for job_script
jorek2_IDS       = "./jorek2_IDS"         # Executable for jorek2_IDS
exclude_file     = "/home/ITER/artolaj/scripts_hub/scripts_javier/exclude_list"   # List with files to exclude when copying folder
out_shot_list    = "shot_list.txt"        # Write shot list with folder names here

###    IMPORTANT: put here the right initial shot number, this could re-write old shots!!!
shot_num         = 111111    # initial shot number
################# END INPUTS #############################################################################################


with open(file_folder_list, 'r') as file:
    # Read the content of the file line by line into a list and strip newline characters
    folder_list = [line.strip() for line in file.readlines()]



job_dir   = os.getcwd()
shot_list = []

for fold in folder_list:

    folder = job_dir + "/" + fold

    print(fold + "     shot = " + str(shot_num))
    shot_list.append(str(shot_num) + ' ' + fold + '\n')

    if (not skip_export):
        # Reading sim_phase instances from the YAML file
        phase_path = folder + "/" + phase_file
        with open(phase_path, 'r') as file:
            loaded_data = yaml.safe_load(file)

        phases = [sim_phase.from_dict(data) for data in loaded_data] 

        # Create temporary folder to copy selected files
        destination_folder = folder + "_tmp/"
        rsync_command = "rsync -av " + folder + "/ " + destination_folder  + " --exclude-from=\'" + exclude_file + "\'"
        os.system(rsync_command)

        # Copy required files to folder
        shutil.copy2(imas_file,  destination_folder + imas_file)
        shutil.copy2(jorek2_IDS, destination_folder + jorek2_IDS)
        shutil.copy2(job_IDS,    destination_folder + job_IDS)

        # Move to the directory
        os.chdir(destination_folder)

        # Add the correct shot number
        replace_parameter(imas_file,"shot_number", str(shot_num))

        first_step = True
        for phase in phases:

            # Create dedicated jorek2_IDS input file to this phase
            phase_nml = "imas_" + phase.name
            shutil.copy2(imas_file, phase_nml)

            # Overwrite shot/entry for the first phase
            if first_step:
                replace_parameter(phase_nml,"overwrite_entry", ".true.")
                first_step = False
            else:
                replace_parameter(phase_nml,"overwrite_entry", ".false.")

            # Adapt h5 range for jorek2_IDS
            replace_parameter(phase_nml,"i_begin", str(phase.h5_start))
            replace_parameter(phase_nml,"i_end",   str(phase.h5_end  ))

            time_interval = (phase.t_end - phase.t_start) / float(phase.n_slices)

            # Copy the h5 files for this phase
            copy_hdf5_files(folder+"/", destination_folder, time_interval, phase.t_start, phase.t_end)

            # Append commands to job script
            run_command = "mpirun -n 1 ./jorek2_IDS < " + phase.input_file + " > output_IDS_" + phase.name
            cp_command  = f"cp {phase_nml} imas.nml"
            with open(job_IDS, 'a') as file:
                file.write(cp_command +"\n")
                file.write(run_command+"\n")

        
        # Launch jobscript
        os.system("sbatch "+ job_IDS)

        # Go back to launching directory
        os.chdir(job_dir)

    shot_num += 1


with open(out_shot_list, 'w') as file_shot_list:
    for shot in shot_list:
        # Write each line to the file
        file_shot_list.write(shot)  # Add a newline character if needed