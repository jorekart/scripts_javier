#  Gets a file from the command line (0D) and calculates the TQ time
#        Usage:   python get_tau_TQ.py thermal_energy.dat  Wth0  t_norm
import numpy as np
import sys

fname = sys.argv[1]
Wth0  = float( sys.argv[2] )
norm  = float( sys.argv[3] )

def find_value_in_list(lista, value):

    tol         = 0.03
    diff        = abs( lista - value )
    min_index   = np.argmin(diff)
    found_value = lista[min_index]
    err         = (found_value-value)/value
    if (err > tol):
        print('Warning, found value differs by more than ' + str(tol*100) + '%') 
    return min_index 

def get_tau_TQ(fname): 

    col_time = 0;    col_Wth   = 1  
    
    try:
        data = np.loadtxt(fname, skiprows=1)
#        data = np.loadtxt(fname)
    except:
        print('Input file not found')
        quit()
    
    percent_top = 0.9;    percent_bot = 0.2
    
    Wth_bot = Wth0 * percent_bot 
    Wth_top = Wth0 * percent_top
    
    times  = data[:,col_time]
    Wths    = data[:,col_Wth]
    
    min_index    = find_value_in_list(Wths, Wth_top)
    t_top        = data[min_index, col_time]
    Wth_top_found = data[min_index, col_Wth]
    
    min_index    = find_value_in_list(Wths, Wth_bot)
    t_bot        = data[min_index, col_time]
    Wth_bot_found = data[min_index, col_Wth]
    
    tau = (t_bot-t_top)/(percent_top - percent_bot) 
 
    return tau

def compute_derivative(x, y):
    # Compute the derivative using numpy's gradient function
    dy_dx = np.gradient(y, x)
    return dy_dx

def normalize_derivative(x, y):
    dy_dx = -compute_derivative(x, y)
    
    # Find the maximum derivative and its corresponding x value
    max_derivative = np.max(dy_dx)
    max_derivative_index = np.argmax(dy_dx)
    max_derivative_x = x[max_derivative_index]
    
    return max_derivative, max_derivative_x

data = np.loadtxt(fname, skiprows=1)  # Assuming comma-separated values and skipping header row

# Extract the columns
x = data[:, 0] *norm 
y = data[:, 1]

# Compute and print the normalized derivative
max_derivative, max_derivative_x = normalize_derivative(x, y)

print('')
print(' (t90-20)/0.7   = ' + str(get_tau_TQ(fname)*norm)+ ' ms')
print(' W0/(dW/dt)_max = ' + str(Wth0/max_derivative)+ ' ms')
print(' t@ (dW/dt)_max = ' + str(max_derivative_x)+ ' ms')
print('')

