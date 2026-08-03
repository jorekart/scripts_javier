#  Gets a file from the command line and calculates a decay time
#        Usage:   python get_decay_time.py -f current.dat -options 
import numpy as np
import argparse

def find_value_in_list(lista, value):

    tol         = 0.03
    diff        = abs( lista - value )
    min_index   = np.argmin(diff)
    found_value = lista[min_index]
    err         = (found_value-value)/value
    if (err > tol):
        print('Warning, found value differs by more than ' + str(tol*100) + '%') 
    return min_index 

def get_tau(times, qtty, pcent_thresholds, use_max_val): 

    percent_top = pcent_thresholds[0];    percent_bot = pcent_thresholds[1]
    
    if use_max_val:
        qtty_ref      = np.max(qtty)
        index_max_val = np.argmax(qtty)
    else: 
        qtty_ref      = qtty[0]
        index_max_val = 0

    qtty_bot = qtty_ref * percent_bot 
    qtty_top = qtty_ref * percent_top

    times_final = times[index_max_val:]
    qtty_final  = qtty[index_max_val:]
    
    min_index = find_value_in_list(qtty_final, qtty_top)
    t_top     = times_final[min_index]
    qtty_top  = qtty_final[min_index]
    print(' ')
    print(f'  Time and qtty at {percent_top}% : {t_top} {qtty_top}')
    
    min_index = find_value_in_list(qtty_final, qtty_bot)
    t_bot     = times_final[min_index]
    qtty_bot  = qtty_final[min_index]
    print(f'  Time and qtty at {percent_bot}% : {t_bot} {qtty_bot}')
    print(' ')

    tau = (t_bot-t_top)/(percent_top - percent_bot) 
 
    return tau


parser = argparse.ArgumentParser(description="Script to get liner decay time of quantity",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("-f", "--file", type=str, default="current.dat", help="Data file")
parser.add_argument("-tn", "--tnorm", type=float, default=1, help="Time normalization factor")
parser.add_argument("-ts", "--tstart", type=float, default=0, help="Only consider data after this time")
parser.add_argument("-tc", "--tcolumn", type=int, default=0, help="Column index for time arrary")
parser.add_argument("-qc", "--qcolumn", type=int, default=1, help="Column index for quantity arrary")
parser.add_argument("-f1", "--tfrac1", type=float, default=0.8, help="Quantity percentage to calculate t1")
parser.add_argument("-f2", "--tfrac2", type=float, default=0.2, help="Quantity percentage to calculate t2")
parser.add_argument("-qm", "--maxq", type=bool, default=False, help="Use max qqty value as reference value?")

args = parser.parse_args()

# Extract the columns from the data file
fname = args.file
data  = np.loadtxt(fname, skiprows=1)  # Assuming comma-separated values and skipping header row
time  = data[:, args.tcolumn] * args.tnorm 
qtty  = data[:, args.qcolumn]

use_max_val      = args.maxq
pcent_thresholds = [args.tfrac1, args.tfrac2] 

print('')
print(f'  (t{args.tfrac1}-t{args.tfrac2})/{args.tfrac1-args.tfrac2}   = ' + str(get_tau(time, qtty, pcent_thresholds, use_max_val))+ ' ')
print('')

