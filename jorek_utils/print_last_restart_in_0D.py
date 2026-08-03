import argparse
import numpy as np

parser = argparse.ArgumentParser(description="Prints last restart file in 0D file",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-zD", "--zD_file", type=str, default='nonsense', help="0D file")
parser.add_argument("-Rf", "--readme_file", type=str, default='nonsense', help="README file to append output")

args = parser.parse_args()

data = np.loadtxt( args.zD_file ) 
f_restart = data[:,1][ 0]
l_restart = data[:,1][-1]

file = open(args.readme_file, 'a') 
file.write('\n ') 
file.write('First  restart = ' + str(f_restart) + '\n') 
file.write('Last   restart = ' + str(l_restart)+ '\n')  
file.close() 

