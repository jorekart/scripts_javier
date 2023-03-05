import numpy as np
import argparse

######## README ##############################################
#
#  This code creates a SLURM job script template
#
##############################################################

parser = argparse.ArgumentParser(description="Create a job script template",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-p", "--profile", type=str, default='gen10_ib', help="Properties of the machine (gen10_ib, gen9_ib")
parser.add_argument("-t", "--type", type=str, default='hybrid', help="Type of parallelization (MPI, OpenMP, hybrid)")
parser.add_argument("-name", "--name", type=str, default='jorek', help="Name of the job")
parser.add_argument("-n", "--nodes", type=int, default=1, help="Number of requested nodes")
parser.add_argument("-c", "--command", type=str, default=" ", help="Command to be run")
parser.add_argument("-mpi_pn", "--mpi_per_node", type=int, default=1, help="Number of mpi tasks per node")

args = parser.parse_args()


# Available cluster profiles
queue_name="gen10_ib"
cpus=36
exclude_nodes="c07n[113-120]"

config =   "export SLURM_CPU_BIND=\"rank\""  + "\n" \
         + "export I_MPI_PIN_MODE=lib"       + "\n"      \
         + "export I_MPI_PIN_DOMAIN=auto"    + "\n"      \
         + "export I_MPI_THREAD_LEVEL_DEFAULT=multiple" + "\n"  \
         + "export FI_PROVIDER=verbs"        + "\n"  
        
if (args.profile == 'gen9_ib'):
    queue_name="gen9_ib_centos8"
    cpus=28
    exclude_nodes="none"

# Properties depending on paralellization
if (args.type.lower() == 'hybrid'):
    mpi_tasks_pn  = args.mpi_per_node
    cpus_per_task = int(cpus/args.mpi_per_node)
    omp_threads   = cpus_per_task       

elif (args.type.lower() == 'openmp'):
    mpi_tasks_pn  = 1 
    cpus_per_task = cpus 
    omp_threads   = cpus       

elif (args.type.lower() == 'mpi'):
    mpi_tasks_pn  = cpus 
    cpus_per_task = 1 
    omp_threads   = 1       

# Write the job script
f2 = open("job_template.sh", "w")
f2.write( "#!/bin/bash -l                                        "+"\n")
f2.write( "#SBATCH --time=730:00:00                              "+"\n")
f2.write( "#SBATCH --job-name=" + args.name                       +"\n")
f2.write( "#SBATCH --exclusive                                   "+"\n")

if (exclude_nodes != "none"):
    f2.write( "#SBATCH --exclude=" + exclude_nodes                    +"\n")

f2.write( "#SBATCH --ntasks-per-node=" + str(mpi_tasks_pn)        +"\n")
f2.write( "#SBATCH --cpus-per-task=" + str(cpus_per_task)         +"\n")
f2.write( "#SBATCH --nodes=" + str(args.nodes)                    +"\n")
f2.write( "#SBATCH --partition=" + queue_name                     +"\n")
f2.write( "#SBATCH --error myjob.err                             "+"\n")
f2.write( "                                                      "+"\n")
f2.write( "export OMP_NUM_THREADS=" + str(omp_threads)            +"\n")
f2.write( "                                                      "+"\n")
f2.write( config )
f2.write( "                                                      "+"\n")
f2.write( args.command                                            +"\n")

f2.close()  

