import numpy as np
import argparse
from scipy import integrate


parser = argparse.ArgumentParser(description="Integrate data with x and y from file columns",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-x", "--x", type=int, default=1, help="Column for x", required=True)
parser.add_argument("-y", "--y", type=int, default=2, help="Column for y", required=True)
parser.add_argument("-f", "--file", type=str, default='none', help="File", required=True)
parser.add_argument("-ofact", "--factor", type=float, default=1, help="Multiply output by this factor")
parser.add_argument("-t", "--time", type=str, default='no', help="Export integral as a function of time (yes/no)")

args = parser.parse_args()

data = np.loadtxt(args.file)

x = data[:,args.x-1]
y = data[:,args.y-1]

print( 'The total integral is = ' + str( integrate.simps(y, x) * args.factor ) )
 
if args.time == 'yes':

    n_times = len(x)

    # Start integral to integrate after the first 3 steps
    i_start = 2

    integral_t = np.zeros_like(x)

    # Loop over time
    for i_time in range(n_times):

        # Skip first times
        if (i_time <= i_start): continue

        x_t = x[0:i_time] 
        y_t = y[0:i_time] 

        integral_t[i_time] = integrate.simps(y_t, x_t)  

np.savetxt("integral_vs_time.txt", np.transpose([x, integral_t]) )
 
