#!/usr/bin/env python3
import fileinput
import sys

# This script changes the value of a paramer in a file
# It replaces the RHS (after the = sign) with the value scaled
# Usage
#   python replace_param.py file scaling_factor

def replace_scale_parameter(file_name,var,value):
    with fileinput.FileInput(file_name, inplace=True, backup='.bak') as file:
        for line in file:                         # go over lines
            if var.lower() in line.lower():      # check if the variable is in the line
                split_line = line.split("=",1)    # divide in 2 the line by "=" delimiter
                str_old = split_line[1].split()[0]
                if 'd' in str_old:
                    val_old = float(str_old.replace('d','e'))
                elif 'D' in str_old:
                    val_old = float(str_old.replace('D','e'))
                else: 
                    val_old = float(str_old)
                print(line.replace(split_line[1]," " + str(val_old*float(scale_fact)) + "\n"), end='')  # replace value on the RHS
            else: 
                print(line.replace("1e26275462165", value), end='')  # absurd string to search, not elegant, but needed

fname      = sys.argv[1]
scale_fact = sys.argv[2]

replace_scale_parameter(fname, 'wall_resistivity_fact', scale_fact)
replace_scale_parameter(fname, 'eta ',                  scale_fact)
replace_scale_parameter(fname, 'eta_ohmic',             scale_fact)
replace_scale_parameter(fname, 'D_par',                 scale_fact)
replace_scale_parameter(fname, 'D_perp(1)',             scale_fact)
replace_scale_parameter(fname, 'ZK_par ',               scale_fact)
replace_scale_parameter(fname, 'ZK_par_max',            scale_fact)
replace_scale_parameter(fname, 'ZK_par_min',            scale_fact)
replace_scale_parameter(fname, 'ZK_perp(1)',            scale_fact)
