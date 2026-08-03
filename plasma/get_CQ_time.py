# Gets a file from the command line (0D) and calculates the CQ time
import numpy as np
import sys

fname = sys.argv[1]

def find_value_in_list(lista, value):

    tol         = 0.03
    diff        = abs( lista - value )
    min_index   = np.argmin(diff)
    found_value = lista[min_index]
    err         = (found_value-value)/value
    if (err > tol):
        print('Warning, found value differs by more than ' + str(tol*100) + '%') 
    return min_index 

def get_tau_CQ(fname): 

    Ip0      = 15e6
    col_time = 0;    col_Ip   = 62
    
    try:
    #    data = np.loadtxt(fname, skiprows=1)
        data = np.loadtxt(fname)
    except:
        print('Input file not found')
        quit()
    
    percent_top = 0.8;    percent_bot = 0.2
    
    Ip_bot = Ip0 * percent_bot 
    Ip_top = Ip0 * percent_top
    
    times  = data[:,col_time]
    Ips    = data[:,col_Ip]
    
    min_index    = find_value_in_list(Ips, Ip_top)
    t_top        = data[min_index, col_time]
    Ip_top_found = data[min_index, col_Ip]
    
    min_index    = find_value_in_list(Ips, Ip_bot)
    t_bot        = data[min_index, col_time]
    Ip_bot_found = data[min_index, col_Ip]
    
    tau = (t_bot-t_top)/(percent_top - percent_bot) * 1e3 # to convert to ms
 
    return tau

print('')
print(' tau_CQ = ' + str(get_tau_CQ(fname))+ ' ms')
print('')

