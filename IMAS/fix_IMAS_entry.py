import imas
import sys
sys.path.append('/home/ITER/artolaj/scripts_hub/scripts_javier/IMAS') 
from imas_custom_utils import parse_and_open_imas_entry
import numpy as np

entry_and_args = parse_and_open_imas_entry() 

imas_entry = entry_and_args["entry"]
args       = entry_and_args["args"]

summary = imas_entry.get('summary')

#summary.global_quantities.ip.value[0] = 444

#imas_entry.put(summary)

imas_entry.close()