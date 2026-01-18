import numpy as np
import matplotlib.pyplot as plt
import re
from pathlib import Path

def extract_from_name(filename):
    """
    Extract Te [keV] and ne (in your filename units) from:
    <elem>_Te<X>keV_ne<Y>_time<Z>ms
    """
    Te = None
    ne = None

    mTe = re.search(r"_Te([\d\.]+)keV", filename)
    if mTe:
        Te = float(mTe.group(1))

    mne = re.search(r"_ne([\d\.]+)", filename)
    if mne:
        ne = float(mne.group(1))

    return Te, ne

def load_output_file(path):
    """
    File format: two columns [t_ms, ratio]
    """
    data = np.loadtxt(path)
    t_ms = data[:, 0]
    ratio = data[:, 1]
    return t_ms, ratio

# ----------------------------
# Provide your 5 files here
# ----------------------------

elem = "xenon"
files = [
    f"{elem}_Te0.50keV_ne0.35_time10.00ms",
    f"{elem}_Te1.00keV_ne0.35_time10.00ms",
    f"{elem}_Te2.50keV_ne0.35_time10.00ms",
    f"{elem}_Te5.00keV_ne0.35_time10.00ms",
    f"{elem}_Te10.00keV_ne0.35_time10.00ms",
]

plt.figure()

for f in files:
    t_ms, ratio = load_output_file(f)

    Te, ne = extract_from_name(Path(f).name)
    if Te is None:
        label = Path(f).name
    else:
        # label only by Te (keV); add ne if you like
        label = fr"$T_e$={Te:.2f} keV"
        # label = fr"$T_e$={Te:.2f} keV, $n_e$={ne:.2f}"  # optional

    plt.plot(t_ms, ratio, label=label)

plt.xlabel("Time [ms]")
plt.ylabel("Radiation / Coronal radiation")
plt.yscale("log")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()
