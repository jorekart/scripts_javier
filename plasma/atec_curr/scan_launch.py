"""
Scan of one input parameter: writes one input file per value and runs the simulations in series
(meant to be called inside a batch job). Each run writes sol_<label>.h5 and log_<label>.txt,
with <label> = <parameter>_<value>. scan.yaml describes the scan for scan_collect.py.

Example:
    python scan_launch.py input.yaml physics.D_n 1 5 10 50 100 500 --outdir /scratch/users/$USER/atec_scans/D_n
"""
import argparse
import subprocess
import sys
import time
from pathlib import Path

import yaml

SOLVER = Path(__file__).resolve().parent / "1d_solver.py"


def parse_value(text):
    """'50' -> 50, '1.0e+19' or '1e19' -> float (YAML would read '1e19' as a string)"""
    try:
        return int(text)
    except ValueError:
        pass
    try:
        return float(text)
    except ValueError:
        sys.exit(f"scan value '{text}' is not a number (all arguments after the parameter name are scan values)")


def set_param(params, path, value):
    """Set params["physics"]["D_n"] = value for path = "physics.D_n" (the parameter must exist)"""
    keys = path.split(".")
    section = params
    for key in keys[:-1]:
        section = section[key]
    if keys[-1] not in section:
        raise KeyError(f"parameter {path} not found in the input file")
    section[keys[-1]] = value


def run_label(path, value):
    return f"{path.split('.')[-1]}_{value:g}"


def main():
    parser = argparse.ArgumentParser(description="Run a scan of one input parameter, simulations in series")
    parser.add_argument("base_input", help="base input file (yaml)")
    parser.add_argument("param", help="scanned parameter as section.name, e.g. physics.D_n")
    parser.add_argument("values", nargs="+", help="values of the scanned parameter")
    parser.add_argument("--outdir", default=".", help="scan directory for inputs, solutions and logs (default: current directory)")
    parser.add_argument("--overwrite", action="store_true", help="allow existing sol_*.h5 of this scan to be overwritten")
    args = parser.parse_args()

    with open(args.base_input) as f:
        base_params = yaml.safe_load(f)

    outdir = Path(args.outdir).resolve()
    values = [parse_value(v) for v in args.values]
    labels = [run_label(args.param, v) for v in values]
    if len(set(labels)) != len(labels):
        raise ValueError(f"repeated values in the scan: {labels}")

    existing = [label for label in labels if (outdir / f"sol_{label}.h5").exists()]
    if existing and not args.overwrite:
        sys.exit(f"solutions already exist in {outdir} for {existing}; use --overwrite or another --outdir")

    outdir.mkdir(parents=True, exist_ok=True)

    # one input file per value
    for value, label in zip(values, labels):
        params = yaml.safe_load(yaml.safe_dump(base_params))      # deep copy
        set_param(params, args.param, value)
        with open(outdir / f"input_{label}.yaml", "w") as f:
            yaml.safe_dump(params, f, sort_keys=False)

    # description of the scan, read by scan_collect.py
    scan_info = {"param": args.param, "values": values, "labels": labels,
                 "files": [f"sol_{label}.h5" for label in labels],
                 "base_input": str(Path(args.base_input).resolve())}
    with open(outdir / "scan.yaml", "w") as f:
        yaml.safe_dump(scan_info, f, sort_keys=False)

    print(f"Scan of {args.param} over {values} in {outdir}", flush=True)

    # run the simulations in series
    failed = []
    for label in labels:
        print(f"{time.strftime('%F %T')} start {label}", flush=True)
        with open(outdir / f"log_{label}.txt", "w") as log:
            result = subprocess.run([sys.executable, str(SOLVER), f"input_{label}.yaml", f"sol_{label}.h5"],
                                    cwd=outdir, stdout=log, stderr=subprocess.STDOUT)
        if result.returncode == 0:
            last_line = (outdir / f"log_{label}.txt").read_text().strip().splitlines()[-1]
            print(f"{time.strftime('%F %T')} done  {label}: {last_line}", flush=True)
        else:
            failed.append(label)
            print(f"{time.strftime('%F %T')} FAILED {label}, see log_{label}.txt", flush=True)

    print(f"Scan finished: {len(labels) - len(failed)}/{len(labels)} runs completed"
          + (f", failed: {failed}" if failed else ""))


if __name__ == "__main__":
    main()
