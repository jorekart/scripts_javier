#!/usr/bin/env python3
import fileinput
import sys
import re
import numpy as np

# This script changes the value of a paramer in a file
# It replaces the RHS (after the = sign) with the value scaled
# Usage
#   python replace_param.py file 

def replace_scale_parameter(file_name,var, scale_fact):
    p1 = r'\s*=\s*'   # Pattern for = sign
    p2 = r'([-+]?\d+(?:\.\d+)?(?:[eEdD][-+]?\d+)?)' # Pattern for number
    pattern = var + p1 + p2

    with fileinput.FileInput(file_name, inplace=True, backup='.bak') as file:
        for line in file:                         # go over lines
            match = re.search(pattern, line, flags=re.IGNORECASE )
            if match:     
                # Get pattern with updated value
                found_pattern = match.group()
                split_line    = found_pattern.split("=",1)
                str_old       = split_line[1]
                str_old       = str_old.replace('d', 'e')
                str_old       = str_old.replace('D', 'e')
                new_val       = float(str_old)*scale_fact
                new_pattern   = split_line[0] + " = %s"%( '{:18.10e}'.format(new_val))
            
                # Update line
                spn    = line.split(found_pattern)
                line2  = spn[0] + new_pattern + spn[1]
                print(line2, end='')  # replace value on the RHS
            else:
                print(line.replace("1e26275462165", str(scale_fact)), end='')  # absurd string to search, not elegant, but needed

    try: 
        print(var + ': ' + str_old + '--> ' + str(new_val) )
    except:
        print(var+ ' not found!!!!!!!!! ') 

fname = sys.argv[1]

scale_R = 0.5

pattern_Rbnd = r"R_boundary\s*\(\s*(\d+)\s*\)"
pattern_Zbnd = r"Z_boundary\s*\(\s*(\d+)\s*\)"
pattern_psi  = r"psi_boundary\s*\(\s*(\d+)\s*\)"
pattern_pf   = r"pf_coils\s*\(\s*(\d+)\s*\)\s*\%\s*current"
pattern_wall1= r"rc_w\s*\(\s*(\d+)\s*\)"
pattern_wall2= r"rs_w\s*\(\s*(\d+)\s*\)"
pattern_wall3= r"zc_w\s*\(\s*(\d+)\s*\)"
pattern_wall4= r"zs_w\s*\(\s*(\d+)\s*\)"

# Scale grid parameters
replace_scale_parameter(fname,      'R_geo', scale_R)
replace_scale_parameter(fname,      'Z_geo', scale_R)
replace_scale_parameter(fname, pattern_Rbnd, scale_R)
replace_scale_parameter(fname, pattern_Zbnd, scale_R)

# Scale psi
scale_psi = scale_R**3
replace_scale_parameter(fname,     pattern_psi, scale_psi)
replace_scale_parameter(fname, 'psi_axis_init', scale_psi)

# Scale F0
scale_F0  = scale_R**2
replace_scale_parameter(fname, 'F0 ', scale_F0)

# Scale rho
scale_rho = scale_R**2
replace_scale_parameter(fname,     'rho_0', scale_rho)
replace_scale_parameter(fname,     'rho_1', scale_rho)

# Scale conduction
scale_ZK = scale_rho
replace_scale_parameter(fname,     'ZK_perp\(1\)', scale_ZK)
replace_scale_parameter(fname,         'ZK_par', scale_ZK)
replace_scale_parameter(fname,     'ZK_par_max', scale_ZK)
replace_scale_parameter(fname,     'ZK_par_min', scale_ZK)

# Scale viscosity (due to R factor and rho change)
scale_visco = scale_R**2 * scale_rho
replace_scale_parameter(fname,         'visco', scale_ZK)

# Scale PF coil currents
scale_pf = scale_R**2
replace_scale_parameter(fname,   pattern_pf, scale_pf)

# Scale wall
replace_scale_parameter(fname,   pattern_wall1, scale_R)
replace_scale_parameter(fname,   pattern_wall2, scale_R)
replace_scale_parameter(fname,   pattern_wall3, scale_R)
replace_scale_parameter(fname,   pattern_wall4, scale_R)
replace_scale_parameter(fname,'eta_thin_w', 1.0/scale_R)




# Scale FFprime
scale_ffp = scale_R
fdat = np.loadtxt('jorek_ffprime')
np.savetxt('jorek_ffprime_scaled', np.transpose([ fdat[:,0], fdat[:,1]*scale_ffp] ))


