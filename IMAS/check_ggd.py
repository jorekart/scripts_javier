import imas
from imas import imasdef
import matplotlib.pyplot as plt
import sys
sys.path.append('/home/ITER/artolaj/scripts_hub/scripts_javier/IMAS') 
from imas_custom_utils import parse_and_open_imas_entry
import numpy as np

entry_and_args = parse_and_open_imas_entry() 

imas_entry = entry_and_args["entry"]

#ids = imas_entry.get_slice('plasma_profiles',0.05,imasdef.PREVIOUS_INTERP,occurrence=1)

ids = imas_entry.partial_get('plasma_profiles','ggd(5)/r_j_total_phi',occurrence=1)

ids = imas_entry.partial_get('plasma_profiles','grid_ggd(0)',occurrence=1)

print(ids)


