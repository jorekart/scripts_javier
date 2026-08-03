#!/usr/bin/env python3
"""Convert MDCUBE SPI fragment output to the IMAS SPI IDS."""

from __future__ import annotations

import argparse
import time
from pathlib import Path

import numpy as np


# Hard-coded pellet/core composition data not present in the MDCUBE file.
# Set PELLET_CORE_ATOMS_N to the original pellet atom inventory.
PELLET_CORE_SPECIES_NAME = "H"
PELLET_CORE_SPECIES_Z_N = 1.0
PELLET_CORE_SPECIES_A = 1.0
PELLET_CORE_SPECIES_DENSITY = 5.2e28  # m^-3 (atom density not molecular)
PELLET_DIAMETER = 19.0* 1e-3  # m
PELLET_L_OVER_D = 2.0  # ratio of pellet length to diameter

PELLET_VELOCITY_SHATTER = 300  # m/s
INJECTOR_SHATTERING_ANGLE = 15.5  # degrees

PELLET_LENGTH = PELLET_L_OVER_D * PELLET_DIAMETER  # m
PELLET_RADIUS = PELLET_DIAMETER / 2.0
PELLET_CORE_ATOMS_N = np.pi * PELLET_RADIUS**2 * PELLET_LENGTH * PELLET_CORE_SPECIES_DENSITY

# The origin of the MD-CUBE fragment positions (the exit of the DFW or shattering point into the torus)
# https://user.iter.org/default.aspx?uid=4T7TMR
# We take EQ_02L_4_1 as representative
INJECTOR_NAME = "EQ_02L_4"
ORIGIN_GLOBAL_XYZ_M = np.array([7511.084535, 3880.419966, 685.5], dtype=float) * 1e-3  # m
ORIGIN_EXTRA_R_OFFSET_M = 0.3  # m, positive moves the local origin outward in global +R, this shifts the local MDCUBE coordinates to the global origin (recessed from the DFW exit point)

# Local MDCUBE axes expressed in the global cylindrical basis (R, phi, Z).
# This right-handed default means:
#   local x points inward, in the -R direction
#   local y points upward, in the +Z direction
#   local z points toroidally, in the +phi direction
LOCAL_X_IN_GLOBAL_RPHIZ = np.array([-1.0, 0.0, 0.0])
LOCAL_Y_IN_GLOBAL_RPHIZ = np.array([0.0, 0.0, 1.0])
LOCAL_Z_IN_GLOBAL_RPHIZ = np.array([0.0, 1.0, 0.0])

# Code metadata for the IMAS IDS
CODE_NAME = "MD-CUBE"
CODE_DESCRIPTION = "Calculates realistic 3D fragment distributions from shattered pellets, see [IDM_5KX6U2]"
RESULT_REPORT = ", see [IDM_5KX6U2]"

# Minimal hard-coded plasma context required by the database summary IDS.
SUMMARY_DESCRIPTION = "ITER SPI MDCUBE fragment case"
SUMMARY_B0_T = -5.3
SUMMARY_R0_M = 6.2
SUMMARY_IP_A = -15.0e6


def read_mdcube_fragments(
    filename: str | Path,
    origin_global_xyz_m=ORIGIN_GLOBAL_XYZ_M,
    volume_fraction=1.0,
    max_fragments=None,
):
    """Read MDCUBE rows and return global positions in m, velocities in m/s and volumes in m^3."""
    data = np.loadtxt(filename, comments="#")
    data = np.atleast_2d(data)
    if data.shape[1] < 9:
        raise ValueError("Expected at least 9 columns: Frag_ID #Partls Equiv_dia COM_x COM_y COM_z COM_vx COM_vy COM_vz")

    diameter_m = data[:, 2] * 1.0e-3
    volume_m3 = (4.0 / 3.0) * np.pi * (0.5 * diameter_m) ** 3
    data, volume_m3 = filter_fragments_by_volume(data, volume_m3, volume_fraction, max_fragments)

    local_xyz_m = data[:, 3:6] * 1.0e-3
    local_vxyz_ms = data[:, 6:9]
    global_xyz_m, global_vxyz_ms = local_to_global_cartesian(
        local_xyz_m,
        local_vxyz_ms,
        origin_global_xyz_m,
    )
    return global_xyz_m, global_vxyz_ms, volume_m3


def filter_fragments_by_volume(data, volume_m3, volume_fraction=1.0, max_fragments=None):
    """Keep the largest fragments that account for the requested volume fraction."""
    if not 0.0 < volume_fraction <= 1.0:
        raise ValueError("--volume-fraction must be > 0 and <= 1")

    order = np.argsort(volume_m3)[::-1]
    sorted_volume = volume_m3[order]
    total_volume = float(np.sum(sorted_volume))

    if total_volume > 0.0 and volume_fraction < 1.0:
        cumulative_fraction = np.cumsum(sorted_volume) / total_volume
        keep_count = int(np.searchsorted(cumulative_fraction, volume_fraction, side="left") + 1)
    else:
        keep_count = len(order)

    if max_fragments is not None:
        keep_count = min(keep_count, int(max_fragments))

    keep = order[:keep_count]
    keep = np.sort(keep)
    kept_fraction = float(np.sum(volume_m3[keep]) / total_volume) if total_volume > 0.0 else 0.0
    print(f"Keeping {keep_count} of {len(data)} fragments ({kept_fraction:.6f} of condensed volume)")
    return data[keep], volume_m3[keep]


def local_basis_in_global_cartesian(origin_global_xyz_m):
    """Return local MDCUBE basis vectors expressed in global Cartesian XYZ."""
    origin_global_xyz_m = np.asarray(origin_global_xyz_m, dtype=float)
    phi0 = np.arctan2(origin_global_xyz_m[1], origin_global_xyz_m[0])
    e_r = np.array([np.cos(phi0), np.sin(phi0), 0.0])
    e_phi = np.array([-np.sin(phi0), np.cos(phi0), 0.0])
    e_z = np.array([0.0, 0.0, 1.0])
    global_cyl_basis = np.vstack((e_r, e_phi, e_z))
    local_axes_rphiz = np.vstack((
        LOCAL_X_IN_GLOBAL_RPHIZ,
        LOCAL_Y_IN_GLOBAL_RPHIZ,
        LOCAL_Z_IN_GLOBAL_RPHIZ,
    ))
    local_axes_xyz = local_axes_rphiz @ global_cyl_basis
    handedness = float(np.linalg.det(local_axes_xyz))
    if handedness < 0.0:
        raise ValueError("The configured local-to-global axes are left-handed. Flip one axis to restore right-handedness.")
    return local_axes_xyz


def shifted_origin_global_xyz(origin_global_xyz_m):
    """Move the hard-coded origin by an extra global radial offset."""
    origin_global_xyz_m = np.asarray(origin_global_xyz_m, dtype=float)
    phi0 = np.arctan2(origin_global_xyz_m[1], origin_global_xyz_m[0])
    e_r = np.array([np.cos(phi0), np.sin(phi0), 0.0])
    return origin_global_xyz_m + ORIGIN_EXTRA_R_OFFSET_M * e_r


def local_to_global_cartesian(local_xyz_m, local_vxyz_ms, origin_global_xyz_m):
    """Shift local MDCUBE coordinates to the global origin and rotate axes into global XYZ."""
    origin_global_xyz_m = shifted_origin_global_xyz(origin_global_xyz_m)
    local_axes_xyz = local_basis_in_global_cartesian(origin_global_xyz_m)
    global_xyz_m = origin_global_xyz_m + np.asarray(local_xyz_m, dtype=float) @ local_axes_xyz
    global_vxyz_ms = np.asarray(local_vxyz_ms, dtype=float) @ local_axes_xyz
    return global_xyz_m, global_vxyz_ms


def cartesian_to_cylindrical(xyz, vxyz):
    """Convert a right-handed Cartesian basis to IMAS cylindrical components.

    Assumes phi = atan2(y, x), with e_phi = -sin(phi) e_x + cos(phi) e_y.
    """
    x, y, z = xyz.T
    vx, vy, vz = vxyz.T
    r = np.hypot(x, y)
    r_safe = np.where(r > 0.0, r, 1.0)
    phi = np.arctan2(y, x)
    velocity_r = (x * vx + y * vy) / r_safe
    velocity_phi = (x * vy - y * vx) / r_safe
    return r, phi, z, velocity_r, velocity_phi, vz


def weighted_mean(values, weights, axis=0):
    """Return a finite weighted mean, requiring at least one positive weight."""
    weights = np.asarray(weights, dtype=float)
    weight_sum = float(np.sum(weights))
    if weight_sum <= 0.0:
        raise ValueError("Cannot calculate a centre-of-mass quantity with zero total fragment volume")
    return np.average(np.asarray(values, dtype=float), axis=axis, weights=weights)


def fill_species(species_parent, include_density=False):
    species_parent.species.resize(1)
    species = species_parent.species[0]
    species.name = PELLET_CORE_SPECIES_NAME
    if include_density:
        species.density = PELLET_CORE_SPECIES_DENSITY
    species.z_n = PELLET_CORE_SPECIES_Z_N
    species.a = PELLET_CORE_SPECIES_A


def fill_ids_metadata(ids, mdcube_file):
    ids.ids_properties.comment = Path(mdcube_file).name + RESULT_REPORT
    ids.code.name = CODE_NAME
    ids.code.description = CODE_DESCRIPTION


def build_spi_ids(
    mdcube_file,
    dd_version=None,
    origin_global_xyz_m=ORIGIN_GLOBAL_XYZ_M,
    time_s=0.0,
    volume_fraction=1.0,
    max_fragments=None,
):
    import imas

    factory = imas.IDSFactory(version=dd_version) if dd_version else imas.IDSFactory()
    spi = factory.spi()
    try:
        spi.ids_properties.homogeneous_time = imas.ids_defs.IDS_TIME_MODE_HOMOGENEOUS
    except Exception:
        spi.ids_properties.homogeneous_time = 1
    fill_ids_metadata(spi, mdcube_file)

    xyz_m, vxyz_ms, volume_m3 = read_mdcube_fragments(
        mdcube_file,
        origin_global_xyz_m=origin_global_xyz_m,
        volume_fraction=volume_fraction,
        max_fragments=max_fragments,
    )
    r, phi, z, velocity_r, velocity_phi, velocity_z = cartesian_to_cylindrical(xyz_m, vxyz_ms)
    time = np.asarray([time_s], dtype=float)

    spi.time = time
    spi.injector.resize(1)
    injector = spi.injector[0]
    injector.name = INJECTOR_NAME
    injector.shattering_angle = INJECTOR_SHATTERING_ANGLE
    injector.pellet.velocity_shatter = PELLET_VELOCITY_SHATTER
    injector.pellet.diameter = PELLET_DIAMETER
    injector.pellet.length = PELLET_LENGTH

    condensed_atoms_n = float(PELLET_CORE_SPECIES_DENSITY * np.sum(volume_m3))
    pellet_atoms_n = float(PELLET_CORE_ATOMS_N) if PELLET_CORE_ATOMS_N is not None else condensed_atoms_n
    gas_atoms_n = max(pellet_atoms_n - condensed_atoms_n, 0.0)

    print(f"Pellet atoms = {PELLET_CORE_ATOMS_N:.3e}")
    print(f"Fragmentation gas atoms = {gas_atoms_n:.3e}")

    injector.pellet.core.atoms_n = pellet_atoms_n
    fill_species(injector.pellet.core, include_density=True)

    injector.fragmentation_gas.atoms_n = gas_atoms_n
    fill_species(injector.fragmentation_gas)

    injector.velocity_mass_centre_fragments_r = weighted_mean(velocity_r, volume_m3)
    injector.velocity_mass_centre_fragments_phi = weighted_mean(velocity_phi, volume_m3)
    injector.velocity_mass_centre_fragments_z = weighted_mean(velocity_z, volume_m3)

    injector.fragment.resize(len(r))
    for i, fragment in enumerate(injector.fragment):
        fragment.position.r = np.asarray([r[i]])
        fragment.position.phi = np.asarray([phi[i]])
        fragment.position.z = np.asarray([z[i]])
        fragment.velocity_r = np.asarray([velocity_r[i]])
        fragment.velocity_phi = np.asarray([velocity_phi[i]])
        fragment.velocity_z = np.asarray([velocity_z[i]])
        fragment.volume = np.asarray([volume_m3[i]])

    return spi


def build_summary_ids(mdcube_file, dd_version=None, time_s=0.0):
    import imas

    factory = imas.IDSFactory(version=dd_version) if dd_version else imas.IDSFactory()
    summary = factory.summary()
    try:
        summary.ids_properties.homogeneous_time = imas.ids_defs.IDS_TIME_MODE_HOMOGENEOUS
    except Exception:
        summary.ids_properties.homogeneous_time = 1
    fill_ids_metadata(summary, mdcube_file)

    summary.description = SUMMARY_DESCRIPTION
    summary.time = np.asarray([time_s], dtype=float)
    summary.global_quantities.b0.value = np.asarray([SUMMARY_B0_T], dtype=float)
    summary.global_quantities.r0.value = np.asarray([SUMMARY_R0_M], dtype=float)
    summary.global_quantities.ip.value = np.asarray([SUMMARY_IP_A], dtype=float)
    return summary


def write_ids(ids_list, uri):
    import imas

    entry = imas.DBEntry(uri, "w")
    try:
        print(f"Writing IMAS IDSs to {uri} (put)")
        start = time.perf_counter()
        for ids in ids_list:
            entry.put(ids)
        print(f"put completed in {time.perf_counter() - start:.1f} s")
    finally:
        entry.close()


def parse_args():
    parser = argparse.ArgumentParser(description="Convert MDCUBE SPI fragments to IMAS SPI IDS.")
    parser.add_argument("mdcube_file", help="MDCUBE fragment text file")
    parser.add_argument("--uri", default="imas:hdf5?path=./mdcube_spi", help="Output IMAS URI")
    parser.add_argument("--dd-version", default=None, help="Optional IMAS Data Dictionary version")
    parser.add_argument(
        "--origin-global-xyz-m",
        nargs=3,
        type=float,
        default=ORIGIN_GLOBAL_XYZ_M,
        metavar=("X", "Y", "Z"),
        help="Global Cartesian XYZ position of the MDCUBE local origin, in m",
    )
    parser.add_argument("--time", type=float, default=0.0, help="Snapshot time in seconds")
    parser.add_argument(
        "--volume-fraction",
        type=float,
        default=1.0,
        help="Keep largest fragments until this fraction of total condensed volume is represented",
    )
    parser.add_argument("--max-fragments", type=int, default=None, help="Keep at most this many largest fragments")
    parser.add_argument("--dry-run", action="store_true", help="Build the IDS without writing it")
    return parser.parse_args()


def main():
    args = parse_args()
    print(f"Reading MDCUBE fragments from {args.mdcube_file}")
    spi = build_spi_ids(
        args.mdcube_file,
        args.dd_version,
        args.origin_global_xyz_m,
        args.time,
        args.volume_fraction,
        args.max_fragments,
    )
    summary = build_summary_ids(args.mdcube_file, args.dd_version, args.time)
    print(f"Writing IMAS summary and SPI IDSs to {args.uri}")
    if not args.dry_run:
        write_ids([summary, spi], args.uri)
    print(f"Converted {len(spi.injector[0].fragment)} MDCUBE fragments to SPI IDS")


if __name__ == "__main__":
    main()
