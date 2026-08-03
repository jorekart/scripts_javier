import imas
import sys
sys.path.append('/home/ITER/artolaj/scripts_hub/scripts_javier/IMAS') 
from imas_custom_utils import parse_and_open_imas_entry
import numpy as np

entry_and_args = parse_and_open_imas_entry() 

imas_entry = entry_and_args["entry"]
args       = entry_and_args["args"]

summary = imas_entry.get('summary')


time   = summary.time 
Te_avg = summary.volume_average.t_e.value
import numpy as np

def tq_bounds(time, Te, end_frac=0.10, tau=1e-3, rate_frac=0.3):
    """
    End: first time Te < end_frac * Te0
    Start: first time dTe/dt <= -rate_frac *Te0/tau for consecutive samples
    """
    time = np.asarray(time, float)
    Te   = np.asarray(Te,   float)

    Te0 = Te[0]
    end_thr = end_frac * Te0

    ends = np.where(Te < end_thr)[0]
    if ends.size == 0:
        return None, None
    i_end = int(ends[0])

    # dTe/dt (works for non-uniform time grids)
    dTedt = np.gradient(Te, time)

    # expected linear-drop slope magnitude over tau
    rate_thr = -rate_frac * Te0 / tau  # negative

    # Count backward while "fast drop" holds
    i = i_end
    while i >= 0:
        if dTedt[i] > rate_thr:          
            i_start = i + 1
            break
        i -= 1
    else:
        # We never exited; treat start as beginning of array
        i_start = 0

    # Safety clamp
    i_start = int(max(0, min(i_start, i_end)))

    return float(time[i_start]), float(time[i_end])
  
print(time)
print(Te_avg)

t_start, t_end = tq_bounds(time, Te_avg)
print("TQ start:", t_start, "TQ end:", t_end)

imas_entry.close()