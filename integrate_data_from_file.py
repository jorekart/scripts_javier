import numpy as np
import argparse
from scipy import integrate


parser = argparse.ArgumentParser(description="Integrate data with x and y from file columns",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-x", "--x", type=int, default=1, help="Column for x", required=True)
parser.add_argument("-y", "--y", type=int, default=2, help="Column for y", required=True)
parser.add_argument("-f", "--file", type=str, default='none', help="File", required=True)
parser.add_argument("-ofact", "--factor", type=float, default=1, help="Multiply output by this factor")

args = parser.parse_args()

data = np.loadtxt(args.file)

x = data[:,args.x-1]
y = data[:,args.y-1]

print( 'The integral is = ' + str( integrate.simps(y, x) * args.factor ) )
 
