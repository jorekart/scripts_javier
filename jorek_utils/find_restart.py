#!/usr/bin/env python3
import sys
from pathlib import Path

import numpy as np


def main() -> int:
    if len(sys.argv) != 4:
        print("Usage: python find_restart.py <time_file> <time_in_ms> <time_normalization_s>")
        return 1

    fname = Path(sys.argv[1])

    try:
        tsearch = float(sys.argv[2])
        tnorm = float(sys.argv[3])
    except ValueError:
        print("Error: <time_in_ms> and <time_normalization_s> must be numbers.")
        return 1

    if not fname.is_file():
        print(f"Error: file not found: {fname}")
        return 1

    try:
        data = np.loadtxt(fname, skiprows=1)
    except Exception as exc:
        print(f"Error reading {fname}: {exc}")
        return 1

    if data.ndim == 1:
        data = data.reshape(1, -1)

    if data.shape[1] < 2:
        print("Error: input file must contain at least two columns.")
        return 1

    restart_ids = data[:, 0]
    times = data[:, 1] * tnorm * 1e3

    diff = np.abs(times - tsearch)
    min_index = np.argmin(diff)

    print()
    print(
        f"t = {tsearch} corresponds to restart file {int(restart_ids[min_index])} "
        f"with t = {times[min_index]}"
    )
    print()

    return 0


if __name__ == "__main__":
    raise SystemExit(main())