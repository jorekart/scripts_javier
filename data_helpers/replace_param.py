#!/usr/bin/env python3
import fileinput
import sys

# This script changes the value of a paramer in a file
# It replaces the RHS (after the = sign) with the given value
# Usage
#   python replace_param.py file parameter_name value_to_assign

def replace_parameter(file_name,var,value):
    with fileinput.FileInput(file_name, inplace=True, backup='.bak') as file:
        for line in file:                         # go over lines
            if var in line:                       # check if the variable is in the line
                split_line = line.split("=",1)    # divide in 2 the line by "=" delimiter
                print(line.replace(split_line[1]," " + value + "\n"), end='')  # replace value on the RHS
            else: 
                print(line.replace("1e26275462165", value), end='')  # absurd string to search, not elegant, but needed

fname = sys.argv[1]
var   = sys.argv[2]
value = sys.argv[3]

replace_parameter(fname, var, value)
