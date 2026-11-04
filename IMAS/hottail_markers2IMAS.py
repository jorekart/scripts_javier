#!/usr/bin/env python3
"""
Convert a source HDF5 marker file to IMAS distributions IDS (markers representation).

Expected source datasets
------------------------
Required:
    t       : shape (nt,)
    R       : shape (nt, nm)
    Z       : shape (nt, nm)
    phi     : shape (nt, nm)
    p       : shape (nt, nm)   # momentum norm
    xi      : shape (nt, nm)   # cos(pitch angle)
    weight  : shape (nt, nm)

Optional:
    i_elm   : shape (nt, nm)   # currently not stored in standard distributions/markers

Output
------
Writes one IMAS `distributions` IDS containing:
    distribution[0].markers[it].time
    distribution[0].markers[it].coordinate_identifier[:]
    distribution[0].markers[it].weights[:]
    distribution[0].markers[it].positions[:, :]

Coordinates used in positions are:
    [R, Z, phi, p, xi]
with IMAS coordinate identifiers:
    [r, z, phi, momentum, pitch]

A minimal `summary` IDS (description, b0, r0, ip) is written alongside, as
required by the database. Defaults are those of mdcube_imas.py and can be
changed with --summary-description, --b0, --r0 and --ip.

Notes
-----
- If your source `p` is normalized by m_e*c, set `--p-normalized-by-mec`.
  Then the script converts it to SI momentum [kg*m/s].
- If you prefer to keep the normalized value unchanged, omit that flag.
- Species metadata is minimal and can be adjusted to your use case.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import h5py
import numpy as np

import imas


# Physical constants (SI)
M_E = 9.1093837015e-31      # electron mass [kg]
C = 299792458.0             # speed of light [m/s]

# Minimal plasma context required by the database summary IDS (defaults from mdcube_imas.py).
SUMMARY_DESCRIPTION = "ITER hot-tail runaway electron markers case"
SUMMARY_B0_T = -5.3
SUMMARY_R0_M = 6.2
SUMMARY_IP_A = -15.0e6


def _read_dataset(h5: h5py.File, name: str) -> np.ndarray:
    """Read a dataset and return it as a NumPy array."""
    if name not in h5:
        raise KeyError(f"Missing dataset '{name}' in input file.")
    return np.asarray(h5[name])


def _validate_shapes(
    t: np.ndarray,
    R: np.ndarray,
    Z: np.ndarray,
    phi: np.ndarray,
    p: np.ndarray,
    xi: np.ndarray,
    weight: np.ndarray,
) -> tuple[int, int]:
    """Validate source array shapes and return (nt, nm)."""
    if t.ndim != 1:
        raise ValueError(f"'t' must be 1D, got shape {t.shape}")

    arrays_2d = {
        "R": R,
        "Z": Z,
        "phi": phi,
        "p": p,
        "xi": xi,
        "weight": weight,
    }

    shapes = {name: arr.shape for name, arr in arrays_2d.items()}
    unique_shapes = set(shapes.values())
    if len(unique_shapes) != 1:
        raise ValueError(f"2D arrays do not all have the same shape: {shapes}")

    nt, nm = next(iter(unique_shapes))
    if t.shape[0] != nt:
        raise ValueError(
            f"Length of t ({t.shape[0]}) does not match time dimension of marker arrays ({nt})."
        )

    return nt, nm


def _set_identifier(node, value) -> None:
    """
    Assign an IMAS identifier.
    IMAS-Python supports assigning the identifier enum name as a string.
    """
    node = value  # kept for clarity; not used directly


def _fill_coordinate_identifiers(marker_slice) -> None:
    """
    Set the coordinate identifiers for columns of positions[:, :].
    Order must match the columns we build later.
    """
    names = ["r", "z", "phi", "momentum", "pitch", "energy_kinetic"]
    marker_slice.coordinate_identifier.resize(len(names))
    for i, name in enumerate(names):
        # IMAS-Python smart assignment fills name/index/description
        marker_slice.coordinate_identifier[i] = name



def build_distributions_ids(
    input_h5: Path,
    comment: str = "",
    dd_version: str | None = None,
):
    """Build and return an IMAS distributions IDS from the source HDF5 file."""
    # Create IDS factory
    factory = imas.IDSFactory(dd_version) if dd_version else imas.IDSFactory()
    ids = factory.distributions()

    # markers(itime) has its own time coordinate, so heterogeneous time is appropriate
    ids.ids_properties.homogeneous_time = imas.ids_defs.IDS_TIME_MODE_HOMOGENEOUS
    ids.ids_properties.comment = comment or "Converted from external HDF5 marker file"

    with h5py.File(input_h5, "r") as h5:
        group_path = "/groups/001"
        if group_path not in h5:
            raise KeyError(f"Missing group '{group_path}' in input file.")

        grp = h5[group_path]

        t = np.asarray(grp["t"], dtype=float)
        R = np.asarray(grp["R"], dtype=float)
        Z = np.asarray(grp["Z"], dtype=float)
        phi = np.asarray(grp["phi"], dtype=float)
        p = np.asarray(grp["p"], dtype=float)
        xi = np.asarray(grp["xi"], dtype=float)
        weight = np.asarray(grp["weight"], dtype=float)

        nt, nm = _validate_shapes(t, R, Z, phi, p, xi, weight)

        ids.time.resize(nt)
        ids.time[:] = t
        # Create one distribution object
        ids.distribution.resize(1)
        dist = ids.distribution[0]

        # Minimal metadata
        dist.gyro_type = 2  # 2 = gyro centre, per DD description
        dist.species.type.name = "electron"
        dist.species.type.index = 1  # per DD description, 1 = electrons
        
        # One markers AoS entry per time
        dist.markers.resize(nt)

        for it in range(nt):
            ms = dist.markers[it]
            ms.coordinate_identifier.resize(6)  # [r, z, phi, |momentum|, pitch, kinetic energy]

            ms.coordinate_identifier[0].index = 4    # r
            ms.coordinate_identifier[1].index = 3    # z
            ms.coordinate_identifier[2].index = 5    # phi
            ms.coordinate_identifier[3].index = 200  # momentum
            ms.coordinate_identifier[4].index = 403  # pitch
            ms.coordinate_identifier[5].index = 301  # kinetic energy (optional extra column for reference)

            _fill_coordinate_identifiers(ms)

            # weights: [n_markers]
            ms.weights = np.ascontiguousarray(weight[it, :], dtype=np.float64)

            # positions: [n_markers, n_coordinates]
            # Column order must match coordinate_identifier:
            # [r, z, phi, momentum, pitch]
            positions = np.column_stack(
                [
                    R[it, :],
                    Z[it, :],
                    phi[it, :],
                    p[it, :] * M_E*C,
                    xi[it, :],
                    (np.sqrt(1+p[it,:]**2) - 1)*M_E*C**2,  # kinetic energy in SI units, optional extra column for reference
                ]
            ).astype(np.float64, copy=False)

            if positions.shape != (nm, 6):
                raise RuntimeError(
                    f"Internal error: positions has shape {positions.shape}, expected {(nm, 6)}"
                )

            ms.positions = np.ascontiguousarray(positions)

    return ids


def build_summary_ids(
    time: np.ndarray,
    comment: str = "",
    dd_version: str | None = None,
    description: str = SUMMARY_DESCRIPTION,
    b0: float = SUMMARY_B0_T,
    r0: float = SUMMARY_R0_M,
    ip: float = SUMMARY_IP_A,
):
    """Build the minimal summary IDS required by the database, on the given time base."""
    factory = imas.IDSFactory(dd_version) if dd_version else imas.IDSFactory()
    summary = factory.summary()

    summary.ids_properties.homogeneous_time = imas.ids_defs.IDS_TIME_MODE_HOMOGENEOUS
    summary.ids_properties.comment = comment or "Converted from external HDF5 marker file"

    time = np.asarray(time, dtype=float)
    summary.description = description
    summary.time = time
    # b0 and ip are time-dependent in the DD, taken constant over the time base;
    # r0 is a scalar (FLT_0D)
    summary.global_quantities.b0.value = np.full(time.shape, b0, dtype=float)
    summary.global_quantities.r0.value = float(r0)
    summary.global_quantities.ip.value = np.full(time.shape, ip, dtype=float)
    return summary


def write_ids(ids_list, output_uri: str) -> None:
    """
    Write the IDSs to an IMAS data entry.
    Examples:
        imas:hdf5?path=./my_imas_entry
        output.nc
    """
    dbentry = imas.DBEntry(output_uri, "w")
    try:
        for ids in ids_list:
            # Optional but useful: catch common schema issues before writing
            ids.validate()
            dbentry.put(ids)
    finally:
        # Some backends expose close(), some are fine without it; be defensive.
        close = getattr(dbentry, "close", None)
        if callable(close):
            close()


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Convert a marker HDF5 file to IMAS distributions IDS."
    )
    parser.add_argument("input_h5", type=Path, help="Input HDF5 file")
    parser.add_argument(
        "output",
        help=(
            "Output IMAS target. "
            "Examples: 'imas:hdf5?path=./imas_entry' or 'distributions.nc'"
        ),
    )
    parser.add_argument(
        "--dd-version",
        default=None,
        help="Optional IMAS Data Dictionary version, e.g. 4.1.1",
    )
    parser.add_argument(
        "--comment",
        default="Converted from external HDF5 marker file",
        help="Comment stored in ids_properties.comment",
    )
    parser.add_argument(
        "--summary-description",
        default=SUMMARY_DESCRIPTION,
        help="Description stored in summary.description",
    )
    parser.add_argument(
        "--b0",
        type=float,
        default=SUMMARY_B0_T,
        help=f"Vacuum toroidal field at r0 for the summary IDS, in T (default: {SUMMARY_B0_T})",
    )
    parser.add_argument(
        "--r0",
        type=float,
        default=SUMMARY_R0_M,
        help=f"Reference major radius for the summary IDS, in m (default: {SUMMARY_R0_M})",
    )
    parser.add_argument(
        "--ip",
        type=float,
        default=SUMMARY_IP_A,
        help=f"Plasma current for the summary IDS, in A (default: {SUMMARY_IP_A})",
    )

    return parser.parse_args(argv)


def main(argv: list[str]) -> int:
    args = parse_args(argv)

    if not args.input_h5.exists():
        print(f"ERROR: input file does not exist: {args.input_h5}", file=sys.stderr)
        return 2

    try:
        ids = build_distributions_ids(
            input_h5=args.input_h5,
            comment=args.comment,
            dd_version=args.dd_version,
        )
        summary = build_summary_ids(
            time=ids.time,
            comment=args.comment,
            dd_version=args.dd_version,
            description=args.summary_description,
            b0=args.b0,
            r0=args.r0,
            ip=args.ip,
        )
        write_ids([summary, ids], args.output)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    print("Conversion completed successfully.")
    print(f"Input : {args.input_h5}")
    print(f"Output: {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))