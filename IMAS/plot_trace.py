import imas
import matplotlib.pyplot as plt
import argparse
import numpy as np


# MANAGEMENT OF INPUT ARGUMENTS
# ------------------------------
parser = argparse.ArgumentParser(description=\
        '---- Display scenario')
parser.add_argument('-s','--shot',help='Shot number', required=True,type=int)
parser.add_argument('-r','--run',help='Run number',required=True,type=int)
parser.add_argument('-q','--quantity',help='Quantity (Ip/Zaxis)',required=True,type=str)
parser.add_argument('-u','--user_or_path',help='User or absolute path name where the data-entry is located', required=False)
parser.add_argument('-d','--database',help='Database name where the data-entry is located', required=False)

args = vars(parser.parse_args())

shot = args["shot"]
run  = args["run"]
qtty = args["quantity"]

# User or absolute path name
if args['user_or_path'] != None:
    user = args['user_or_path']
else:
    user = 'public'

# Database name
if args['database'] != None:
    database = args['database']
else:
    database = 'iter'
    

imas_entry_init = imas.DBEntry(imas.imasdef.MDSPLUS_BACKEND, database, shot, run, user, data_version = '3')
imas_entry_init.open()

idslist = {}

summary = imas_entry_init.get('summary')

time = summary.time*1e3



if (qtty=="Ip"):
    value = summary.global_quantities.ip.value * 1e-6
    qtty_name = 'Ip [MA]'
elif (qtty=="Zaxis"):
    value = summary.local.magnetic_axis.position.z
    qtty_name = 'Zaxis [m]'

      
plt.plot(time,value,'*-')

plt.ylabel(qtty_name)
plt.xlabel('Time [ms]')

# Show the plot
plt.show()