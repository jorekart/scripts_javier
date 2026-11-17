"""Scan of L_parallel (x_max) and constant Jphi with 1d_solver.py.

Each case runs in its own directory <scan_dir>/L<x_max>_J<Jphi>/ with a copy of
the base input.yaml in which only mesh.x_max and physics.Jphi are changed.
No fort.444 is copied, so the solver uses the constant Jphi and the YAML x_max.

Usage:
    python scan_Lpar_Jphi.py [--scan-dir DIR] [--base-input input.yaml]
"""
import argparse
import copy
import os
import subprocess
import sys

import yaml

L_PAR = [30.0, 60.0, 120.0, 240.0]
JPHI = [0.5e6, 1.0e6, 1.5e6, 2.0e6]

HERE = os.path.dirname(os.path.abspath(__file__))


def case_name(L, J):
    return f"L{L:g}_J{J:.2e}"


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--scan-dir", default="/scratch/users/artolaj/1d_fem_Lpar_Jphi_scan")
    parser.add_argument("--base-input", default=os.path.join(HERE, "input.yaml"))
    args = parser.parse_args()

    with open(args.base_input) as f:
        base = yaml.safe_load(f)

    os.makedirs(args.scan_dir, exist_ok=True)
    solver = os.path.join(HERE, "1d_solver.py")

    for L in L_PAR:
        for J in JPHI:
            run_dir = os.path.join(args.scan_dir, case_name(L, J))
            os.makedirs(run_dir, exist_ok=True)

            params = copy.deepcopy(base)
            params["mesh"]["x_max"] = float(L)
            params["physics"]["Jphi"] = float(J)
            params["physics"].pop("Jphi_profile", None)
            with open(os.path.join(run_dir, "input.yaml"), "w") as f:
                yaml.safe_dump(params, f, sort_keys=False)

            print(f"Running {case_name(L, J)} in {run_dir}", flush=True)
            with open(os.path.join(run_dir, "solver.log"), "w") as log:
                subprocess.run([sys.executable, solver], cwd=run_dir,
                               stdout=log, stderr=subprocess.STDOUT, check=True)


if __name__ == "__main__":
    main()
