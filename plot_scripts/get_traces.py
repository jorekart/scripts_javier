import os
import argparse
from datetime import datetime


######## README ##############################################
#
#  This code does the following
#
#       - creates a job script to run jorek2_postproc
#       - Submits a job for postproc
#       - copies the result into the  plot_traces folder
#       - creates there a README file with some info
#
##############################################################

parser = argparse.ArgumentParser(description="Collect 0D data",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-q", "--queue", type=str, default='gen10_ib', help="Properties of the machine (gen10_ib, gen9_ib")
parser.add_argument("-p_file", "--p_file", type=str, default='nonsense', help="jorek2_postproc executable")
parser.add_argument("-jinp", "--jorek_input", type=str, default='nonsense', help="JOREK input file")
parser.add_argument("-rj", "--restart_jump", type=int, default=20, help="jump restart files by this step")

args = parser.parse_args()

# General directory of scripts
dir_scripts ="~/scripts_hub/scripts_javier/" 

# output folder
dir_out = "get_traces"

# Create directory to save the data
try:
    os.system("rm -r ./"+dir_out)
except:
    pass

os.mkdir(dir_out)

# Create 0D file
f2 = open("pinp_0D", "w")
f2.write( "namelist "+ args.jorek_input +"\n")
f2.write( "si-units" +"\n")
f2.write( "for step 0 to 99999 by "+str(args.restart_jump)+ " do" +"\n")
f2.write( "zeroD_quantities" + "\n")
f2.write( "done"+"\n" )
f2.close()

# postproc command
run_postproc = args.p_file + " < " + "pinp_0D > output_postproc \n"

# cp postproc into get_traces folder
run_cp  = "cp postproc/zeroD_quantities_s00000..99999.dat "+ dir_out + "/0D.dat \n"
run_cp2 = "mv output_postproc "+ dir_out +"\n"

# export info of running directory and time
cwd       = os.getcwd()
now       = datetime.now()
dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
f1 = open(dir_out+"/README", "w")
f1.write( "Time of calling postproc "+ dt_string +"\n")
f1.write( "Calling postproc from folder "+ cwd +"\n")
f1.close()

# Put last restart file into README
run_last_restart = "python "+dir_scripts+"print_last_restart_in_0D.py -zD " + dir_out+"/0D.dat -Rf "+dir_out+"/README "

# all commands together
command =  "\""+ run_postproc + run_cp + run_cp2 + run_last_restart + "\""

# create job script
create_job = "python " +dir_scripts+"job_creator.py -t openmp -p " + args.queue + " -c " + command
os.system(create_job)

# submit job_script
os.system("sbatch job_template.sh")
