#!/usr/bin/env python3
"""
dream_to_imas.py

Map a DREAM HDF5 output file into selected IMAS IDSs with IMAS-Python.

The converter focuses on DREAM data normally useful for disruption/runaway
studies and writes, when available:

  * runaway_electrons IDS:
      /eqsys/n_re, /eqsys/j_re, /other/fluid/runawayRate,
      /other/fluid/gammaDreicer, /gammaAvalanche, /gammaTritium,
      /gammaCompton, /pCrit, /pCritHottail, /EDreic, /Eceff, /Ecfree,
      /Ectot, /other/fluid/W_re

  * plasma_profiles IDS:
      /eqsys/T_cold, /eqsys/n_cold, /eqsys/n_tot, /eqsys/n_i,
      /eqsys/j_tot, /eqsys/j_ohm, /eqsys/E_field,
      /other/fluid/Zeff, /conductivity, /grid/r, /grid/t

  * equilibrium IDS:
      /eqsys/I_p, /eqsys/psi_p, /eqsys/psi_edge,
      /eqsys/psi_wall, /grid/R0, /grid/r, /grid/eq/*

  * summary IDS, when supported by the installed Data Dictionary:
      /eqsys/I_p, /eqsys/V_loop_w, /eqsys/W_cold, /eqsys/psi_edge

Design choices
--------------
DREAM and IMAS do not have a perfect one-to-one mapping. This script therefore:

  1. Maps high-confidence physical quantities to documented IMAS fields.
  2. Uses a defensive `set_path()` helper so the same script can run across
     Data Dictionary versions. If a target node does not exist in your DD,
     it is skipped and recorded in the report instead of crashing.
  3. Writes a JSON mapping report next to the output, including successful,
     skipped, missing-source, and uncertain mappings.
  4. Keeps values in DREAM units where these match IMAS nodes. DREAM T_cold is
     eV, densities are m^-3, current densities are A m^-2, E_field is V m^-1,
     time is seconds, and r/R0 are metres.

Requirements
------------
    pip install h5py numpy imas-python

For writing to an IMAS HDF5 data entry you also need an IMAS-Core installation
with the hdf5 backend available. Without IMAS-Core, use `--backend netcdf`.

Examples
--------
Write to a local IMAS HDF5 data-entry directory:

    python dream_to_imas.py output_CQ_S6.h5 --uri 'imas:hdf5?path=./imas_dream'

Write to a portable netCDF file:

    python dream_to_imas.py output_CQ_S6.h5 --uri dream_imas.nc

Only build IDSs and produce a report, without writing:

    python dream_to_imas.py output_CQ_S6.h5 --dry-run
"""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Optional

import h5py
import numpy as np

try:
    import imas
except Exception as exc:  # pragma: no cover - depends on user environment
    imas = None  # type: ignore[assignment]
    IMAS_IMPORT_ERROR = exc
else:
    IMAS_IMPORT_ERROR = None

# ------------------------ normalization factors---------------------------
psi_cocos = 2.0*np.pi # good
curr_dir  = 1.0       # 
gamma_adb = 5.0/3.0
c_light    = 299792458.0
m_electron = 9.10938356e-31
p_norm = m_electron * c_light


# ----------------------------- reporting ---------------------------------

@dataclass
class MappingReport:
    source_file: str
    written_uri: Optional[str] = None
    dd_version_requested: Optional[str] = None
    created_ids: list[str] = field(default_factory=list)
    set_nodes: list[dict[str, str]] = field(default_factory=list)
    skipped_nodes: list[dict[str, str]] = field(default_factory=list)
    missing_sources: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)

    def ok(self, ids: str, target: str, source: str) -> None:
        self.set_nodes.append({"ids": ids, "target": target, "source": source})

    def skip(self, ids: str, target: str, reason: str, source: str = "") -> None:
        self.skipped_nodes.append(
            {"ids": ids, "target": target, "source": source, "reason": reason}
        )

    def missing(self, source: str) -> None:
        if source not in self.missing_sources:
            self.missing_sources.append(source)

    def warn(self, message: str) -> None:
        self.warnings.append(message)

    def write(self, path: Path) -> None:
        path.write_text(json.dumps(self.__dict__, indent=2), encoding="utf-8")


# ----------------------------- HDF5 helpers -------------------------------

class DreamH5:
    def __init__(self, filename: str | Path):
        self.filename = str(filename)
        self.h5 = h5py.File(self.filename, "r")

    def close(self) -> None:
        self.h5.close()

    def exists(self, path: str) -> bool:
        return path in self.h5

    def arr(self, path: str, report: MappingReport | None = None) -> Optional[np.ndarray]:
        if path not in self.h5:
            if report is not None:
                report.missing(path)
            return None
        data = self.h5[path][()]
        return decode_h5_value(data)

    def scalar(self, path: str, default: float | None = None) -> Optional[float]:
        if path not in self.h5:
            return default
        value = decode_h5_value(self.h5[path][()])
        arr = np.asarray(value)
        if arr.size == 0:
            return default
        return float(arr.reshape(-1)[0])


def decode_h5_value(value: Any) -> Any:
    """Decode bytes/character arrays that appear in DREAM HDF5 files."""
    arr = np.asarray(value)
    if arr.dtype.kind in {"S", "U"}:
        try:
            if arr.ndim == 0:
                return arr.item().decode() if isinstance(arr.item(), bytes) else str(arr.item())
            # DREAM often stores strings as arrays of single characters.
            chars = []
            for x in arr.reshape(-1):
                if isinstance(x, bytes):
                    chars.append(x.decode(errors="ignore"))
                else:
                    chars.append(str(x))
            return "".join(chars).rstrip("\x00")
        except Exception:
            return arr.astype(str)
    return arr


def flatten_1d(x: Optional[np.ndarray]) -> Optional[np.ndarray]:
    if x is None:
        return None
    a = np.asarray(x)
    if a.size == 0:
        return a.astype(float)
    return np.asarray(a, dtype=float).reshape(-1)


def time_aligned(data: Optional[np.ndarray], nt: int) -> Optional[np.ndarray]:
    """Return data aligned to a time dimension, padding edge values if needed."""
    if data is None:
        return None
    a = np.asarray(data, dtype=float)
    if a.ndim == 0:
        return np.full((nt,), float(a))
    if a.shape[0] == nt:
        return a
    if a.shape[0] == nt - 1:
        pad = np.take(a, [-1], axis=0)
        return np.concatenate([a, pad], axis=0)
    if a.shape[0] > nt:
        return a[:nt]
    if a.shape[0] < nt and a.shape[0] > 0:
        reps = [1] * a.ndim
        reps[0] = nt - a.shape[0]
        pad = np.repeat(np.take(a, [-1], axis=0), reps[0], axis=0)
        return np.concatenate([a, pad], axis=0)
    return a


def normalized_radius(r: np.ndarray) -> np.ndarray:
    r = np.asarray(r, dtype=float)
    if r.size == 0:
        return r
    r0 = float(r[0])
    denom = float(r[-1] - r0)
    if abs(denom) < 1e-30:
        return np.zeros_like(r)
    return (r - r0) / denom

# ----------------------------- IMAS helpers -------------------------------

def ensure_imas() -> None:
    if imas is None:
        raise RuntimeError(
            "Could not import imas. Install IMAS-Python / load your IMAS environment. "
            f"Original import error: {IMAS_IMPORT_ERROR!r}"
        )


def make_factory(dd_version: str | None = None):
    ensure_imas()
    if dd_version:
        return imas.IDSFactory(version=dd_version)
    return imas.IDSFactory()


def make_ids(factory: Any, name: str, report: MappingReport):
    if not hasattr(factory, name):
        report.warn(f"Installed IMAS Data Dictionary does not provide IDS '{name}'.")
        return None
    ids = getattr(factory, name)()
    report.created_ids.append(name)
    set_homogeneous_time(ids)
    try:
        ids.ids_properties.comment = (
            "Generated from DREAM HDF5 output by dream_to_imas.py. "
            "Mappings are conservative and reported in the companion JSON file."
        )
    except Exception:
        pass
    return ids


def set_homogeneous_time(ids: Any) -> None:
    try:
        ids.ids_properties.homogeneous_time = imas.ids_defs.IDS_TIME_MODE_HOMOGENEOUS
    except Exception:
        try:
            ids.ids_properties.homogeneous_time = 1
        except Exception:
            pass


def resize_aos(aos: Any, n: int) -> bool:
    try:
        aos.resize(int(n))
        return True
    except Exception:
        return False


def set_path(root: Any, path: str, value: Any, report: MappingReport, ids_name: str, source: str) -> bool:
    """Set a slash-separated IMAS path if it exists in this DD version.

    Example: set_path(pp, 'profiles_1d/0/electrons/temperature', value, ...)

    Array-of-structures elements are handled by numeric path components.
    """
    parts = [p for p in path.split("/") if p]
    if not parts:
        report.skip(ids_name, path, "empty target path", source)
        return False
    obj = root
    try:
        for part in parts[:-1]:
            if part.isdigit():
                obj = obj[int(part)]
            else:
                obj = getattr(obj, part)
        leaf = parts[-1]
        if not hasattr(obj, leaf):
            report.skip(ids_name, path, "target node is not present in this DD version", source)
            return False
        setattr(obj, leaf, sanitize_for_imas(value))
        report.ok(ids_name, path, source)
        return True
    except Exception as exc:
        report.skip(ids_name, path, f"assignment failed: {type(exc).__name__}: {exc}", source)
        return False


def sanitize_for_imas(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        if value.dtype.kind in {"S", "U", "O"}:
            return value.astype(str)
        return np.asarray(value)
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    return value


def resize_child_aos(parent: Any, child_name: str, n: int) -> bool:
    if not hasattr(parent, child_name):
        return False
    return resize_aos(getattr(parent, child_name), n)


# ----------------------------- mappings -----------------------------------

def common_grids(dream: DreamH5, report: MappingReport) -> dict[str, Any]:
    time = flatten_1d(dream.arr("/grid/t", report))

    # Minor radius R-Raxis in [m]
    r  = flatten_1d(dream.arr("/grid/r", report))
    dr = flatten_1d(dream.arr("/grid/dr", report))

    # R at magnetic axis in [m]
    R0 = dream.scalar("/grid/R0", dream.scalar("/grid/eq/R0", None))

    # Flux averaged toroidal field in [T]
    Bphi_avg = flatten_1d(dream.arr("/grid/geometry/GR0", report))

    # Getting vacuum field from the toroidal field at the outermost radius
    if Bphi_avg is not None and Bphi_avg.size > 0:
        B0 = float(Bphi_avg[-1])
        pass
    
    # Toroidal flux in [Wb]
    phi_tor = flatten_1d(dream.arr("/grid/geometry/toroidalFlux", report))

    rho_tor = np.sqrt(phi_tor / (np.pi * B0)) # minor radius like
    rho_tor_norm = normalized_radius(rho_tor) # normalized to 1 at the outermost radius

    # Other geometric quantities
    vpvol  = flatten_1d(dream.arr("/grid/VpVol", report))
    G      = flatten_1d(dream.arr("/grid/geometry/GR0", report))
    R2inv  = flatten_1d(dream.arr("/grid/geometry/FSA_R02OverR2", report))
    Bmin   = flatten_1d(dream.arr("/grid/geometry/Bmin", report))
    Bmax   = flatten_1d(dream.arr("/grid/geometry/Bmax", report))

    weight_int_area = dr * vpvol * G * R2inv / (2*np.pi*Bmin)
    weight_int_vol  = dr * vpvol * R0
    dVdr =  vpvol * R0

    if time is None or time.size == 0:
        raise ValueError("DREAM file does not contain a valid /grid/t array")
    if r is None or r.size == 0:
        raise ValueError("DREAM file does not contain a valid /grid/r array")

    R_outboard = R0 + r

    return {
        "time": np.asarray(time, dtype=float),
        "r": np.asarray(r, dtype=float),
        "dr": np.asarray(dr, dtype=float),
        "R_outboard": np.asarray(R_outboard, dtype=float),
        "phi_tor": phi_tor,  
        "rho_tor": rho_tor,
        "rho_tor_norm": rho_tor_norm,
        "R0": R0,
        "a": dream.scalar("/grid/a", None),
        "B0": B0,
        "dVdr": dVdr,
        "weight_int_area": weight_int_area,
        "weight_int_vol": weight_int_vol
    }


def fill_vacuum_toroidal_field(ids: Any, ids_name: str, grids: dict[str, Any], report: MappingReport) -> None:
    # DREAM output structure shown here does not include B0 explicitly. We fill R0
    # where possible and leave b0 empty unless the user supplies a post-processing hook.
    R0 = grids.get("R0")
    B0 = grids.get("B0")
    if R0 is not None:
        set_path(ids, "vacuum_toroidal_field/r0", R0, report, ids_name, "/grid/R0")
    else:
        report.warn(
            f"{ids_name}: vacuum_toroidal_field/r0 was not filled because no explicit R0 "
            "dataset was present in the provided DREAM HDF5 structure."
        )
    if B0 is not None:
        set_path(ids, "vacuum_toroidal_field/b0", B0, report, ids_name, "/grid/B0")
    else:
        report.warn(
            f"{ids_name}: vacuum_toroidal_field/b0 was not filled because no explicit B0 "
            "dataset was present in the provided DREAM HDF5 structure."
        )


def fill_1d_grid(p: Any, grids: dict[str, Any], dream: DreamH5, nt: int, it: int, report: MappingReport, ids_name: str) -> None:
    """Fill grid data for a single time step in profiles_1d."""
    rho_tor = grids["rho_tor"]
    rho_tor_norm = grids["rho_tor_norm"]
    set_path(p, "grid/rho_tor", rho_tor, report, ids_name, "/grid/geometry/toroidalFlux")
    set_path(p, "grid/rho_tor_norm", rho_tor_norm, report, ids_name, "derived from rho_tor")

    psi_p = time_aligned(dream.arr("/eqsys/psi_p"), nt)*psi_cocos
    if psi_p is not None:
        psi_arr  = np.asarray(psi_p[it], dtype=float)
        psi_bnd  = float(psi_arr[-1]) if psi_arr.size > 0 else 0.0
        psi_axis = float(psi_arr[0]) if psi_arr.size > 0 else 0.0
        psi_norm = normalized_radius(psi_arr) if psi_arr is not None else None

        set_path(p, "grid/psi", psi_arr, report, ids_name, "/eqsys/psi_p")
        set_path(p, "grid/psi_magnetic_axis", psi_axis, report, ids_name, "/eqsys/psi_p[:,0]")
        set_path(p, "grid/psi_boundary", psi_bnd, report, ids_name, "/eqsys/psi_p[:,-1]")
        
        set_path(p, "grid/rho_pol_norm", psi_norm, report, ids_name, "derived")


def map_plasma_profiles(factory: Any, dream: DreamH5, grids: dict[str, Any], report: MappingReport):
    ids_name = "plasma_profiles"
    pp = make_ids(factory, ids_name, report)
    if pp is None:
        return None

    time = grids["time"]
    nt = len(time)

    set_path(pp, "time", time, report, ids_name, "/grid/t")
    fill_vacuum_toroidal_field(pp, ids_name, grids, report)

    if not resize_child_aos(pp, "profiles_1d", nt):
        report.skip(ids_name, "profiles_1d", "could not resize AoS", "")
        return pp

    # DREAM time-dependent profiles.
    profiles = {
        "electrons/temperature": ("/eqsys/T_cold", dream.arr("/eqsys/T_cold", report)),
        "electrons/density": ("/eqsys/n_tot", dream.arr("/eqsys/n_tot")),
        "electrons/density_thermal": ("/eqsys/n_cold", dream.arr("/eqsys/n_cold", report)),
       # "pressure_ion_total": ("/eqsys/W_i", dream.arr("/eqsys/W_i", report)/(gamma_adb-1.0)), 
        "conductivity_parallel": ("/other/fluid/conductivity", dream.arr("/other/fluid/conductivity")),
        "e_field/parallel": ("/eqsys/E_field", dream.arr("/eqsys/E_field", report)),
        "zeff": ("/other/fluid/Zeff", dream.arr("/other/fluid/Zeff")),
    }

    # n_tot is total electron density; if absent, use n_cold as the main density.
    if profiles["electrons/density"][1] is None:
        profiles["electrons/density"] = ("/eqsys/n_cold", dream.arr("/eqsys/n_cold", report))

    j_tot = dream.arr("/eqsys/j_tot", report)
    j_ohm = dream.arr("/eqsys/j_ohm", report)

    for it in range(nt):
        p = pp.profiles_1d[it]
        set_path(p, "time", time[it], report, ids_name, "/grid/t")

        fill_1d_grid(p, grids, dream, nt, it, report, ids_name)

        for target, (source, data) in profiles.items():
            aligned = time_aligned(data, nt)
            if aligned is not None:
                set_path(p, target, aligned[it], report, ids_name, source)
            elif source:
                report.missing(source)

        if j_tot is not None:
            jt = time_aligned(j_tot, nt)[it]
            set_path(p, "j_total", jt, report, ids_name, "/eqsys/j_tot")
        if j_ohm is not None:
            jo = time_aligned(j_ohm, nt)[it]
            set_path(p, "j_ohmic", jo, report, ids_name, "/eqsys/j_ohm")

        # Ion densities: DREAM /eqsys/n_i has shape (time, charge_state_flat, r).
        fill_ion_profiles_1d(p, dream, it, nt, report, ids_name)

    return pp


def dream_arr_first(dream: Any, *paths: str) -> Any:
    """
    Return the first available DREAM HDF5 dataset among `paths`.

    This lets us prefer /settings/eqsys/n_i/*, while falling back to /ionmeta/*.
    """
    for path in paths:
        try:
            value = dream.arr(path)
        except Exception:
            value = None
        if value is not None:
            return value
    return None


def parse_dream_string_list(value: Any) -> list[str]:
    """
    Parse DREAM string-list datasets.

    DREAM often stores lists as one character array, for example:
        D;T;D_inj_stage_1;Ne_inj_stage_1;

    This returns:
        ["D", "T", "D_inj_stage_1", "Ne_inj_stage_1"]
    """
    if value is None:
        return []

    if isinstance(value, str):
        text = value

    elif isinstance(value, (bytes, bytearray)):
        text = bytes(value).decode("utf-8", errors="ignore")

    else:
        arr = np.asarray(value)

        if arr.dtype.kind in {"S", "U", "O"}:
            chars = []
            for x in arr.ravel():
                if isinstance(x, bytes):
                    chars.append(x.decode("utf-8", errors="ignore"))
                else:
                    chars.append(str(x))
            text = "".join(chars)

        elif arr.dtype.kind in {"u", "i"}:
            text = bytes(arr.ravel().astype(np.uint8)).decode("utf-8", errors="ignore")

        else:
            text = str(value)

    text = text.replace("\x00", "").strip()

    # DREAM convention: semicolon-separated list.
    if ";" in text:
        return [p.strip() for p in text.split(";") if p.strip()]

    # Fallback for already-separated strings.
    for sep in [",", "\n"]:
        if sep in text:
            return [p.strip() for p in text.split(sep) if p.strip()]

    return [text] if text else []


def parse_dream_ion_names(value: Any, n: int | None = None) -> list[str]:
    """
    Backwards-compatible wrapper around parse_dream_string_list().
    """
    names = parse_dream_string_list(value)
    if n is not None and len(names) > n:
        # Do not blindly truncate here. The reconciliation routine below will decide
        # which names are compatible with Z/A.
        return names
    return names


def base_species_label(label: str) -> str:
    """
    Extract the physical species prefix from DREAM labels.

    Examples:
        D                -> D
        T                -> T
        D_inj_stage_1    -> D
        Ne_inj_stage_2   -> Ne
    """
    head = label.split("_", 1)[0]
    if len(head) >= 2 and head[:2].istitle():
        return head[:2]
    return head[:1]


def infer_isotope_mass_number(
    label: str,
    z: int,
    isotope: int | None,
    hydrogen_names: set[str],
    tritium_names: set[str],
) -> int | None:
    """
    Infer isotope mass number A.

    Priority:
      1. /settings/eqsys/n_i/isotopes, if present and positive
      2. explicit tritiumnames / hydrogennames
      3. label prefix H/D/T
      4. DREAM convention: Z=1 defaults to deuterium
    """
    if isotope is not None and isotope > 0:
        return int(isotope)

    if z == 1:
        if label in hydrogen_names:
            return 1
        if label in tritium_names:
            return 3

        base = base_species_label(label)
        if base == "H":
            return 1
        if base == "D":
            return 2
        if base == "T":
            return 3

        # DREAM default for Z=1 if not explicitly marked otherwise.
        return 2
    
    if z == 10:
        return 20   # Neon
    if z == 18:
        return 40   # Argon

    return None


def expected_z_from_label(label: str) -> int | None:
    """
    Return nuclear charge inferred from the label prefix, when obvious.

    This is used only to avoid assigning e.g. a D label to a Ne block.
    """
    element_z = {
        "H": 1,
        "D": 1,
        "T": 1,
        "He": 2,
        "Li": 3,
        "Be": 4,
        "B": 5,
        "C": 6,
        "N": 7,
        "O": 8,
        "F": 9,
        "Ne": 10,
        "Ar": 18,
        "W": 74,
    }

    base = base_species_label(label)
    return element_z.get(base)


def reconcile_dream_ion_names(
    raw_names: Any,
    z_list: list[int],
    isotope_list: list[int | None],
    report: Any | None = None,
) -> list[str]:
    """
    Return exactly len(z_list) labels.

    /settings/eqsys/n_i/Z and /eqsys/n_i define the real species blocks.
    /settings/eqsys/n_i/names or /ionmeta/names provide labels, but in some
    injection cases there can be extra labels. We therefore only accept labels
    compatible with the species Z.
    """
    parsed = parse_dream_ion_names(raw_names)
    out: list[str] = []
    used: set[int] = set()

    for i, z in enumerate(z_list):
        chosen = None

        # First try same-position label.
        if i < len(parsed):
            z_label = expected_z_from_label(parsed[i])
            if z_label is None or z_label == z:
                chosen = parsed[i]
                used.add(i)

        # If same-position failed, scan unused labels for a compatible one.
        if chosen is None:
            for j, label in enumerate(parsed):
                if j in used:
                    continue
                z_label = expected_z_from_label(label)
                if z_label is None or z_label == z:
                    chosen = label
                    used.add(j)
                    break

        if chosen is None:
            A = isotope_list[i] if i < len(isotope_list) else None
            chosen = f"species_{i}_Z{z}" if A is None else f"species_{i}_Z{z}_A{A}"
            if report is not None:
                report.warn(
                    f"Could not find a DREAM ion name compatible with species {i}, "
                    f"Z={z}. Parsed names were {parsed}. Using '{chosen}'."
                )

        out.append(chosen)

    if report is not None and len(parsed) != len(z_list):
        report.warn(
            f"Parsed {len(parsed)} DREAM ion names but found {len(z_list)} ion species "
            f"from Z/isotopes. Names were reconciled by compatibility with Z."
        )

    if report is not None:
        unused = [label for j, label in enumerate(parsed) if j not in used]
        if unused:
            report.warn(f"Unused DREAM ion labels not mapped to n_i species blocks: {unused}")

    return out


def set_element_atomic_properties(
    parent: Any,
    label: str,
    z: int,
    A: int | None,
    report: Any,
    ids_name: str,
    source: str,
) -> None:
    """
    Store atomic properties in parent.element[0], where available.

    Intended for both:
        profiles_1d.ion[i].element[0]
        profiles_1d.neutral[i].element[0]
    """
    if not hasattr(parent, "element") or not resize_aos(parent.element, 1):
        report.skip(ids_name, "element", "element AoS not available/resizable", source)
        return

    element = parent.element[0]

    # IMAS commonly uses z_n for nuclear charge and a for mass number.
    set_path(element, "z_n", float(z), report, ids_name, source)

    if A is not None:
        set_path(element, "a", float(A), report, ids_name, source)

    # Useful when available in the DD version.
    set_path(element, "label", label, report, ids_name, source)


def fill_ion_profiles_1d(
    p: Any,
    dream: Any,
    it: int,
    nt: int,
    report: Any,
    ids_name: str,
) -> None:
    """
    Map DREAM ion densities to IMAS plasma_profiles.profiles_1d.

    DREAM:
        /eqsys/n_i has shape (time, charge_state_flat, radius)

    DREAM ordering:
        species 0, Z0=0
        species 0, Z0=1
        ...
        species 0, Z0=Z

        species 1, Z0=0
        species 1, Z0=1
        ...

    Mapping:
        Z0=0      -> p.neutral[i].density
        Z0=1..Z   -> p.ion[i].state[Z0-1].density
        ion.density = sum over charged states only

    Atomic properties:
        p.ion[i].element[0].z_n     = Z
        p.ion[i].element[0].a       = A
        p.neutral[i].element[0].z_n = Z
        p.neutral[i].element[0].a   = A
    """
    n_i = time_aligned(dream_arr_first(dream, "/eqsys/n_i"), nt)

    Z_raw = dream_arr_first(
        dream,
        "/settings/eqsys/n_i/Z",
        "/ionmeta/Z",
    )

    isotopes_raw = dream_arr_first(
        dream,
        "/settings/eqsys/n_i/isotopes",
    )

    names_raw = dream_arr_first(
        dream,
        "/settings/eqsys/n_i/names",
        "/ionmeta/names",
    )

    hydrogen_names_raw = dream_arr_first(
        dream,
        "/settings/eqsys/n_i/hydrogennames",
    )

    tritium_names_raw = dream_arr_first(
        dream,
        "/settings/eqsys/n_i/tritiumnames",
    )

    if n_i is None or Z_raw is None:
        return

    z_arr = flatten_1d(Z_raw)
    if z_arr is None:
        return

    z_list = [int(z) for z in z_arr]

    isotope_arr = flatten_1d(isotopes_raw) if isotopes_raw is not None else None
    isotope_list: list[int | None] = []
    for i in range(len(z_list)):
        if isotope_arr is not None and i < len(isotope_arr):
            A = int(isotope_arr[i])
            isotope_list.append(A if A > 0 else None)
        else:
            isotope_list.append(None)

    labels = reconcile_dream_ion_names(
        names_raw,
        z_list=z_list,
        isotope_list=isotope_list,
        report=report,
    )

    hydrogen_names = set(parse_dream_string_list(hydrogen_names_raw))
    tritium_names = set(parse_dream_string_list(tritium_names_raw))

    A_list = [
        infer_isotope_mass_number(
            label=labels[i],
            z=z_list[i],
            isotope=isotope_list[i],
            hydrogen_names=hydrogen_names,
            tritium_names=tritium_names,
        )
        for i in range(len(z_list))
    ]

    expected = sum(z + 1 for z in z_list)
    if n_i.shape[1] < expected:
        report.warn(
            f"/eqsys/n_i has {n_i.shape[1]} charge-state channels, but "
            f"Z implies {expected}; mapping the available prefix only."
        )
    elif n_i.shape[1] > expected:
        report.warn(
            f"/eqsys/n_i has {n_i.shape[1]} charge-state channels, but "
            f"Z implies {expected}; extra trailing channels will not be mapped."
        )

    has_ion = hasattr(p, "ion") and resize_aos(p.ion, len(z_list))
    if not has_ion:
        report.skip(ids_name, "profiles_1d/ion", "ion AoS not available/resizable", "/eqsys/n_i")

    has_neutral = hasattr(p, "neutral") and resize_aos(p.neutral, len(z_list))
    if not has_neutral:
        report.skip(ids_name, "profiles_1d/neutral", "neutral AoS not available/resizable", "/eqsys/n_i")

    if not has_ion and not has_neutral:
        return

    idx = 0

    for iion, z in enumerate(z_list):
        label = labels[iion]
        A = A_list[iion]
        nstates = z + 1

        block = n_i[it, idx:min(idx + nstates, n_i.shape[1]), :]
        idx += nstates

        if block.size == 0:
            continue

        if has_ion and has_neutral:
            set_path(p.ion[iion], "neutral_index", iion, report, ids_name, "DREAM species index")
            set_path(p.neutral[iion], "ion_index", iion, report, ids_name, "DREAM species index") 

        # ------------------------------------------------------------
        # DREAM Z0=0: neutral density
        # ------------------------------------------------------------
        neutral_density = block[0, :]

        if has_neutral:
            neutral = p.neutral[iion]

            set_path(neutral, "label", label, report, ids_name, "/settings/eqsys/n_i/names")
            set_element_atomic_properties(
                neutral,
                label=label,
                z=z,
                A=A,
                report=report,
                ids_name=ids_name,
                source="/settings/eqsys/n_i/Z,/settings/eqsys/n_i/isotopes",
            )

            # Depending on DD version, neutral density may be direct or under state.
            set_path(neutral, "density", neutral_density, report, ids_name, "/eqsys/n_i Z0=0")

            if hasattr(neutral, "state") and resize_aos(neutral.state, 1):
                nstate = neutral.state[0]
                set_path(nstate, "label", f"{label} neutral", report, ids_name, "/ionmeta/names")
                set_path(nstate, "density", neutral_density, report, ids_name, "/eqsys/n_i Z0=0")

        # ------------------------------------------------------------
        # DREAM Z0=1..Z: charged ion densities
        # ------------------------------------------------------------
        charged_block = block[1:, :]

        if has_ion:
            ion = p.ion[iion]

            set_path(ion, "label", label, report, ids_name, "/settings/eqsys/n_i/names")
            set_element_atomic_properties(
                ion,
                label=label,
                z=z,
                A=A,
                report=report,
                ids_name=ids_name,
                source="/settings/eqsys/n_i/Z,/settings/eqsys/n_i/isotopes",
            )

            # ion.density should exclude neutrals.
            if charged_block.size:
                set_path(
                    ion,
                    "density",
                    np.sum(charged_block, axis=0),
                    report,
                    ids_name,
                    "/eqsys/n_i summed over charged states Z0=1..Z",
                )

            # One IMAS ion.state entry per charged state.
            # DREAM block[1] -> Z0=1, block[2] -> Z0=2, ...
            if hasattr(ion, "state") and resize_aos(ion.state, charged_block.shape[0]):
                for z0 in range(1, block.shape[0]):
                    state = ion.state[z0 - 1]

                    # Do NOT set z_min/z_max here, per your request.
                    # The charge state is encoded by ordering and label.
                    set_path(
                        state,
                        "label",
                        f"{label} Z0={z0}",
                        report,
                        ids_name,
                        "/eqsys/n_i charge-state index",
                    )
                    set_path(
                        state,
                        "density",
                        block[z0, :],
                        report,
                        ids_name,
                        "/eqsys/n_i",
                    )

def map_runaway_electrons(factory: Any, dream: DreamH5, grids: dict[str, Any], report: MappingReport):
    ids_name = "runaway_electrons"
    re_ids = make_ids(factory, ids_name, report)
    if re_ids is None:
        return None

    time = grids["time"]
    nt   = len(time)

    set_path(re_ids, "time", time, report, ids_name, "/grid/t")
    fill_vacuum_toroidal_field(re_ids, ids_name, grids, report)

    if not resize_child_aos(re_ids, "profiles_1d", nt):
        report.skip(ids_name, "profiles_1d", "could not resize AoS", "")
        return re_ids

    if not resize_child_aos(re_ids, "global_quantities", nt):
        report.skip(ids_name, "global_quantities", "could not resize AoS", "")
        return re_ids

    quantities = {
        "density": ("/eqsys/n_re", dream.arr("/eqsys/n_re", report)),
        "current_density": ("/eqsys/j_re", dream.arr("/eqsys/j_re", report)),
        #"energy_density_kinetic": ("/other/fluid/W_re", dream.arr("/other/fluid/W_re")),
        "ddensity_dt_total": ("/other/fluid/runawayRate", dream.arr("/other/fluid/runawayRate")),
        "ddensity_dt_dreicer": ("/other/fluid/gammaDreicer", dream.arr("/other/fluid/gammaDreicer")),
        "ddensity_dt_hot_tail": ("/other/fluid/gammaHottail", dream.arr("/other/fluid/gammaHottail")), 
        "ddensity_dt_tritium": ("/other/fluid/gammaTritium", dream.arr("/other/fluid/gammaTritium")),  
        "ddensity_dt_compton": ("/other/fluid/gammaCompton", dream.arr("/other/fluid/gammaCompton")),  
        #"ddensity_dt_avalanche": ("/other/fluid/GammaAva", dream.arr("/other/fluid/GammaAva")),  #not available in dictionary
        "momentum_critical_avalanche": ("/other/fluid/pCrit", dream.arr("/other/fluid/pCrit")*p_norm),
        "momentum_critical_hot_tail": ("/other/fluid/pCritHottail", dream.arr("/other/fluid/pCritHottail")*p_norm),
        "e_field_dreicer": ("/other/fluid/EDreic", dream.arr("/other/fluid/EDreic")),
        "e_field_critical": ("/other/fluid/Ectot", dream.arr("/other/fluid/Ectot")),
    }

    # Fill up profiles_1d for each time step
    for it in range(nt):
        p = re_ids.profiles_1d[it]
        set_path(p, "time", time[it], report, ids_name, "/grid/t")

        fill_1d_grid(p, grids, dream, nt, it, report, ids_name)

        for target, (source, data) in quantities.items():
            aligned = time_aligned(data, nt)
            if aligned is not None:
                set_path(p, target, aligned[it], report, ids_name, source)
            elif source:
                report.missing(source)

    # Fill global_quantities for each time step, if available.
    weight_int_area = grids.get("weight_int_area")
    j_RE = dream.arr("/eqsys/j_re")
    I_RE = np.sum(j_RE * weight_int_area[None,:], axis=1)

    set_path(re_ids, "global_quantities/current", I_RE, report, ids_name, "derived from j_re")

    return re_ids


import numpy as np

def cell_center_to_edges(xc):
    """
    Construct cell-edge values from cell-center values.

    Interior edges are midpoint averages.
    Boundary edges are linear extrapolations.
    """
    xc = np.asarray(xc)

    xf = np.empty(xc.size + 1, dtype=xc.dtype)

    xf[1:-1] = 0.5 * (xc[:-1] + xc[1:])
    xf[0] = xc[0] - 0.5 * (xc[1] - xc[0])
    xf[-1] = xc[-1] + 0.5 * (xc[-1] - xc[-2])

    return xf


def dVdpsi_from_dVdr(dVdr, dr, psi):
    """
    Convert dV/dr to dV/dpsi while preserving the cell-integrated volume.

    dVdr, dr, psi are all cell-centered arrays with the same length.
    """
    dVdr = np.asarray(dVdr)
    dr = np.asarray(dr)
    psi = np.asarray(psi)

    dV_cell = dVdr * dr

    psi_f = cell_center_to_edges(psi)
    dpsi = np.diff(psi_f)

    dVdpsi = dV_cell / dpsi

    return dVdpsi, dpsi, psi_f


def map_equilibrium(factory: Any, dream: DreamH5, grids: dict[str, Any], report: MappingReport):
    ids_name = "equilibrium"
    eq = make_ids(factory, ids_name, report)
    if eq is None:
        return None

    time = grids["time"]

    nt = len(time)
    R0 = grids.get("R0") or 0.0
    Z0 = dream.scalar("/grid/eq/Z0", 0.0) or 0.0

    set_path(eq, "time", time, report, ids_name, "/grid/t")
    fill_vacuum_toroidal_field(eq, ids_name, grids, report)

    if not resize_child_aos(eq, "time_slice", nt):
        report.skip(ids_name, "time_slice", "could not resize AoS", "")
        return eq

    psi_p = time_aligned(dream.arr("/eqsys/psi_p", report), nt)*psi_cocos
    j_tot = time_aligned(dream.arr("/eqsys/j_tot", report), nt)
    ip = flatten_1d(dream.arr("/eqsys/I_p", report))
    
    r_boundary, z_boundary = boundary_outline(dream, R0, Z0)

    # Some radial grid profiles
    rho_tor = grids["rho_tor"]
    rho_tor_norm = grids["rho_tor_norm"]
    phi_tor = grids["phi_tor"]
    r = grids["r"]
    dr = grids["dr"]

    R_outboard = grids["R_outboard"] 
    dVdr = grids.get("dVdr")
    
    for it in range(nt):
        ts = eq.time_slice[it]
        set_path(ts, "time", time[it], report, ids_name, "/grid/t")
        if ip is not None and it < len(ip):
            set_path(ts, "global_quantities/ip", float(ip[it]), report, ids_name, "/eqsys/I_p")
        
        set_path(ts, "profilesls_1d/rho_tor", rho_tor, report, ids_name, "derived from rho_tor /grid/geometry/toroidalFlux")
        set_path(ts, "profiles_1d/rho_tor_norm", rho_tor_norm, report, ids_name, "derived from rho_tor")
        set_path(ts, "profiles_1d/phi_tor", phi_tor, report, ids_name, "derived from phi_tor")
        set_path(ts, "profiles_1d/r_outboard", R_outboard, report, ids_name, "derived from R0 and r")

        if j_tot is not None:
            set_path(ts, "profiles_1d/j_phi", j_tot[it], report, ids_name, "/eqsys/j_tot")

        if psi_p is not None:
            psi_arr  = np.asarray(psi_p[it], dtype=float)
            psi_bnd  = float(psi_arr[-1]) if psi_arr.size > 0 else 0.0
            psi_axis = float(psi_arr[0]) if psi_arr.size > 0 else 0.0
            psi_norm = normalized_radius(psi_arr) if psi_arr is not None else None

            set_path(ts, "profiles_1d/psi", psi_arr, report, ids_name, "/eqsys/psi_p")
            set_path(ts, "profiles_1d/rho_pol_norm", psi_norm, report, ids_name, "derived")

            set_path(ts, "global_quantities/psi_boundary", psi_bnd, report, ids_name, "/eqsys/psi_p[:,-1]")
            set_path(ts, "profiles_1d/psi_magnetic_axis", psi_axis, report, ids_name, "/eqsys/psi_p[:,0]")

            dVdpsi, dpsi, psi_f = dVdpsi_from_dVdr(dVdr, dr, psi_arr)

            set_path(ts, "profiles_1d/dvolume_dpsi", dVdpsi, report, ids_name, "derived")

        if r_boundary is not None and z_boundary is not None:
            set_path(ts, "boundary/outline/r", r_boundary, report, ids_name, "/grid/eq/RMinusR0_f[:,-1] + R0")
            set_path(ts, "boundary/outline/z", z_boundary, report, ids_name, "/grid/eq/ZMinusZ0_f[:,-1] + Z0")
            set_path(ts, "boundary/geometric_axis/r", float(np.mean([np.min(r_boundary), np.max(r_boundary)])), report, ids_name, "/grid/eq/RMinusR0_f")
            set_path(ts, "boundary/geometric_axis/z", float(np.mean([np.min(z_boundary), np.max(z_boundary)])), report, ids_name, "/grid/eq/ZMinusZ0_f")

    

    return eq


def boundary_outline(dream: DreamH5, R0: float, Z0: float) -> tuple[Optional[np.ndarray], Optional[np.ndarray]]:
    Rm = dream.arr("/grid/eq/RMinusR0_f")
    Zm = dream.arr("/grid/eq/ZMinusZ0_f")
    if Rm is None or Zm is None:
        Rm = dream.arr("/grid/eq/RMinusR0")
        Zm = dream.arr("/grid/eq/ZMinusZ0")
    if Rm is None or Zm is None:
        return None, None
    Rm = np.asarray(Rm, dtype=float)
    Zm = np.asarray(Zm, dtype=float)
    if Rm.ndim != 2 or Zm.ndim != 2:
        return None, None
    r_outline = R0 + Rm[:, -1]
    z_outline = Z0 + Zm[:, -1]
    # Close contour if not already closed.
    if r_outline.size > 1 and (r_outline[0] != r_outline[-1] or z_outline[0] != z_outline[-1]):
        r_outline = np.r_[r_outline, r_outline[0]]
        z_outline = np.r_[z_outline, z_outline[0]]
    return r_outline, z_outline


def map_summary(factory: Any, dream: DreamH5, grids: dict[str, Any], report: MappingReport):
    ids_name = "summary"
    summary = make_ids(factory, ids_name, report)
    if summary is None:
        return None
    time = grids["time"]
    set_path(summary, "time", time, report, ids_name, "/grid/t")
    fields = {
        "global_quantities/ip/value": ("/eqsys/I_p", dream.arr("/eqsys/I_p")),
        "global_quantities/v_loop/value": ("/eqsys/V_loop_w", dream.arr("/eqsys/V_loop_w")),
        "global_quantities/energy_diamagnetic/value": ("/eqsys/W_cold", dream.arr("/eqsys/W_cold")),
        "global_quantities/psi_boundary/value": ("/eqsys/psi_edge", dream.arr("/eqsys/psi_edge")),
    }
    for target, (source, data) in fields.items():
        if data is None:
            report.missing(source)
            continue
        arr = np.asarray(data)
        if arr.ndim > 1:
            if arr.shape[1] == 1:
                arr = arr[:, 0]
            else:
                # For W_cold, use a simple radial sum as a scalar time trace only
                # if the summary path exists. Prefer users verify against volume metric.
                arr = np.sum(arr, axis=1)
                report.warn(f"summary/{target}: radial profile {source} was summed to make a scalar trace.")
        set_path(summary, target, time_aligned(arr, len(time)), report, ids_name, source)
    return summary


# ----------------------------- writing ------------------------------------

def write_ids(ids_list: list[Any], uri: str, report: MappingReport) -> None:
    ensure_imas()
    db = imas.DBEntry(uri, "w")
    try:
        for ids in ids_list:
            if ids is not None:
                db.put(ids)
        report.written_uri = uri
    finally:
        try:
            db.close()
        except Exception:
            pass


def build_ids(dream_file: str, dd_version: str | None, selected: Iterable[str]) -> tuple[list[Any], MappingReport]:
    factory = make_factory(dd_version)
    report = MappingReport(source_file=dream_file, dd_version_requested=dd_version)
    dream = DreamH5(dream_file)
    try:
        grids = common_grids(dream, report)
        ids_list = []
        selected_set = set(selected)
        if "plasma_profiles" in selected_set:
            ids_list.append(map_plasma_profiles(factory, dream, grids, report))
        if "runaway_electrons" in selected_set:
            ids_list.append(map_runaway_electrons(factory, dream, grids, report))
        if "equilibrium" in selected_set:
            ids_list.append(map_equilibrium(factory, dream, grids, report))
        if "summary" in selected_set:
            ids_list.append(map_summary(factory, dream, grids, report))
        return [ids for ids in ids_list if ids is not None], report
    finally:
        dream.close()


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Convert selected DREAM HDF5 output quantities to IMAS IDSs."
    )
    parser.add_argument("dream_h5", help="Input DREAM HDF5 output file")
    parser.add_argument(
        "--uri",
        default="dream_imas.nc",
        help=(
            "Output URI passed to imas.DBEntry. Examples: 'dream_imas.nc' for netCDF, "
            "or 'imas:hdf5?path=./imas_dream' for IMAS-Core HDF5."
        ),
    )
    parser.add_argument(
        "--ids",
        nargs="+",
        default=["plasma_profiles", "runaway_electrons", "equilibrium", "summary"],
        choices=["plasma_profiles", "runaway_electrons", "equilibrium", "summary"],
        help="IDSs to create/write.",
    )
    parser.add_argument(
        "--dd-version",
        default=None,
        help="Optional IMAS Data Dictionary version, e.g. 4.1.0. Default: environment/latest.",
    )
    parser.add_argument(
        "--report",
        default=None,
        help="Path to write JSON mapping report. Default: <uri or input stem>.mapping_report.json",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Build IDS objects and write only the mapping report; do not call DBEntry.put().",
    )
    return parser.parse_args(argv)


def default_report_path(args: argparse.Namespace) -> Path:
    if args.report:
        return Path(args.report)
    uri_name = args.uri.replace(":", "_").replace("?", "_").replace("/", "_").replace(";", "_")
    if uri_name:
        return Path(f"{uri_name}.mapping_report.json")
    return Path(args.dream_h5).with_suffix(".mapping_report.json")


def main(argv: list[str] | None = None) -> int:
    args = parse_args(sys.argv[1:] if argv is None else argv)
    if not Path(args.dream_h5).exists():
        print(f"ERROR: input file does not exist: {args.dream_h5}", file=sys.stderr)
        return 2
    ids_list, report = build_ids(args.dream_h5, args.dd_version, args.ids)
    if not args.dry_run:
        write_ids(ids_list, args.uri, report)
    report_path = default_report_path(args)
    report.write(report_path)
    print(f"Built IDSs: {', '.join(report.created_ids) if report.created_ids else '(none)'}")
    if args.dry_run:
        print("Dry run: no IMAS DBEntry was written.")
    else:
        print(f"Wrote IMAS data entry: {args.uri}")
    print(f"Mapping report: {report_path}")
    print(f"Set nodes: {len(report.set_nodes)} | Skipped nodes: {len(report.skipped_nodes)} | Missing sources: {len(report.missing_sources)}")
    if report.warnings:
        print("Warnings:")
        for w in report.warnings[:10]:
            print(f"  - {w}")
        if len(report.warnings) > 10:
            print(f"  ... {len(report.warnings) - 10} more warnings in the report")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
