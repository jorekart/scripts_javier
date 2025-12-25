import argparse
import sys
import imas
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

sys.path.append('/home/ITER/artolaj/scripts_javier/IMAS')
from imas_custom_utils2 import read_uris, get_entries, make_entry_label


def dd_major(dd_version) -> int:
    try:
        return int(str(dd_version).split(".")[0])
    except Exception:
        return 4

def get_dd_version(entry) -> str:
    try:
        ids = entry.get("summary",autoconvert=False)
        return imas.util.get_data_dictionary_version(ids)
    except Exception:
        print("Failed loading summary to check dd_version...")
        return entry.dd_version
        

def choose_profile_ids(dd_version) -> str:
    return "core_profiles" if dd_major(dd_version) == 3 else "plasma_profiles"


def get_quantity_spec(qtty: str):
    qtty = qtty.strip()
    if qtty == "Te":
        return r"$T_e$ [eV]", ("electrons", "temperature")
    if qtty == "Jphi":
        return r"$J_{\phi}$ [A/m$^2$]", ("j_total",)
    if qtty == "ne":
        return r"$n_e$ [m$^{-3}$]", ("electrons", "density")
    if qtty == "psi":
        return r"$\psi$ [Wb]", ("grid", "psi")
    raise ValueError("Unknown quantity. Use: Te, Jphi, ne, psi")


def slice_time_value(ids_slice) -> float | None:
    if hasattr(ids_slice, "time"):
        t = np.asarray(ids_slice.time).ravel()
        if t.size:
            return float(t[0])
    return None


def get_profiles_1d_node(ids_slice):
    p1d = getattr(ids_slice, "profiles_1d", None)
    if p1d is None:
        return None
    try:
        return p1d[0] if len(p1d) > 0 else None
    except Exception:
        return p1d


def get_grid_attr(grid, attr: str):
    # small helper to safely fetch grid.<attr>
    if hasattr(grid, attr):
        v = getattr(grid, attr)
        if v is not None:
            return v
    return None


def choose_xcoord_from_slice(ids_slice, xcoord: str, warn_prefix: str = ""):
    """
    Returns (x_array, x_label, used_attr_name).
    If requested xcoord isn't available, falls back to alternatives.
    """
    p1d = get_profiles_1d_node(ids_slice)
    if p1d is None:
        raise RuntimeError("Slice has no profiles_1d")

    grid = p1d.grid
    xcoord = xcoord.strip().lower()

    # Candidate lists: first is preferred, later are fallbacks.
    if xcoord == "rho_tor_norm":
        candidates = [
            ("rho_tor_norm", r"$\rho_{tor,norm}$"),
            ("rho_pol_norm", r"$\rho_{pol,norm}$"),
            ("rho_tor",      r"$\rho_{tor}$"),
            ("rho_pol",      r"$\rho_{pol}$"),
            ("psi",          r"$\psi$ [Wb]"),
        ]
    elif xcoord == "rho_pol_norm":
        candidates = [
            ("rho_pol_norm", r"$\rho_{pol,norm}$"),
            ("rho_tor_norm", r"$\rho_{tor,norm}$"),
            ("rho_pol",      r"$\rho_{pol}$"),
            ("rho_tor",      r"$\rho_{tor}$"),
            ("psi",          r"$\psi$ [Wb]"),
        ]
    elif xcoord == "rho_tor":
        candidates = [
            ("rho_tor",      r"$\rho_{tor}$"),
            ("rho_tor_norm", r"$\rho_{tor,norm}$"),
            ("rho_pol",      r"$\rho_{pol}$"),
            ("rho_pol_norm", r"$\rho_{pol,norm}$"),
            ("psi",          r"$\psi$ [Wb]"),
        ]
    elif xcoord == "rho_pol":
        candidates = [
            ("rho_pol",      r"$\rho_{pol}$"),
            ("rho_pol_norm", r"$\rho_{pol,norm}$"),
            ("rho_tor",      r"$\rho_{tor}$"),
            ("rho_tor_norm", r"$\rho_{tor,norm}$"),
            ("psi",          r"$\psi$ [Wb]"),
        ]
    elif xcoord == "psi":
        candidates = [
            ("psi",          r"$\psi$ [Wb]"),
            ("rho_tor_norm", r"$\rho_{tor,norm}$"),
            ("rho_pol_norm", r"$\rho_{pol,norm}$"),
            ("rho_tor",      r"$\rho_{tor}$"),
            ("rho_pol",      r"$\rho_{pol}$"),
        ]
    else:
        raise ValueError("Unknown --xcoord. Use: rho_tor_norm, rho_pol_norm, rho_tor, rho_pol, psi")

    for attr, lab in candidates:
        v = get_grid_attr(grid, attr)
        if v is not None:
            x = np.asarray(v)
            # Optional sanity: accept 1D vectors only
            if x.ndim == 1:
                if attr != xcoord and warn_prefix:
                    print(f"[WARN] {warn_prefix}: requested xcoord='{xcoord}' not available, using '{attr}'")
                return x, lab, attr

    raise RuntimeError(f"None of the xcoord candidates are available for requested '{xcoord}'")


def get_y_from_slice(ids_slice, qtty: str):
    p1d = get_profiles_1d_node(ids_slice)
    if p1d is None:
        raise RuntimeError("Slice has no profiles_1d")

    ylab, qpath = get_quantity_spec(qtty)
    obj = p1d
    for name in qpath:
        obj = getattr(obj, name)
    y = np.asarray(obj)
    if y.ndim != 1:
        raise RuntimeError(f"Expected 1D profile for y, got y.ndim={y.ndim}")
    return y, ylab


def get_timebase(entry, profile_ids: str) -> np.ndarray:
    # Lightweight time read to build common time grid if --times not given
    ids = entry.get(profile_ids, lazy=True, autoconvert=False)
    #ids = imas.convert_ids(ids, entry.dd_version)
    if hasattr(ids, "time"):
        t = np.asarray(ids.time).ravel()
        if t.size:
            return t
    raise RuntimeError(f"{profile_ids} has no usable time array")


def compute_common_times(time_arrays: list[np.ndarray], ncurves: int) -> np.ndarray:
    tmins = [float(np.min(t)) for t in time_arrays]
    tmaxs = [float(np.max(t)) for t in time_arrays]
    tmin_common = max(tmins)
    tmax_common = min(tmaxs)

    if tmax_common <= tmin_common:
        tmin_common = min(tmins)
        tmax_common = max(tmaxs)
        print("[WARN] Time ranges do not overlap; using full union time window.")

    ncurves = max(1, int(ncurves))
    if ncurves == 1:
        return np.array([(tmin_common + tmax_common) / 2.0], dtype=float)
    return np.linspace(tmin_common, tmax_common, ncurves)


def main():
    parser = argparse.ArgumentParser(
        description="Plot IMAS 1D profiles: colors=requested time (rainbow), linestyles=entry; uses get_slice(CLOSEST)."
    )
    src = parser.add_mutually_exclusive_group(required=True)
    src.add_argument("--uri", "-uri", type=str, help="Single IMAS URI")
    src.add_argument("--uri_file", "-uf", type=str, help="Text file with one IMAS URI per line")

    parser.add_argument("-q", "--quantity", default="Te", help="Te / Jphi / ne / psi")

    parser.add_argument("--xcoord", default="rho_tor_norm",
                        help="x-axis coordinate: rho_tor_norm / rho_pol_norm / rho_tor / rho_pol / psi")

    parser.add_argument("--times", nargs="+", type=float, default=None,
                        help="Requested times [s]. If omitted, uses --ncurves over common time window.")
    parser.add_argument("--ncurves", type=int, default=10,
                        help="Number of requested times if --times is not provided.")
    parser.add_argument("--warn_dt", type=float, default=5e-3,
                        help="Warn if |t_selected - t_requested| > warn_dt [s]. Default: 5e-3")
    parser.add_argument("--dd_version", default="None",
                        help="Pass-through to get_entries (e.g. 'None' or '4.0.0').")

    args = parser.parse_args()

    uris = read_uris(uri=args.uri, uri_file=args.uri_file)
    entries = get_entries(uris, mode="r", dd_version=args.dd_version)

    datasets = []
    time_arrays = []
    linestyles = ["-", "--", ":", "-."]

    # Track whether we've already warned about xcoord fallback for an entry (avoid spam)
    warned_xcoord = set()

    for i, (entry, uri) in enumerate(zip(entries, uris), start=1):
        dd_version = get_dd_version(entry)
        print(f" Getting uri={uri} with version = {dd_version}")
        profile_ids = choose_profile_ids(dd_version)
        print(f' Choosen profile IDS is {profile_ids}')

        try:
            t = get_timebase(entry, profile_ids)
            time_arrays.append(t)
        except Exception as e:
            print(f"[WARN] Skipping {make_entry_label(uri, i)} (cannot read timebase): {e}")
            continue

        datasets.append({
            "i": i,
            "uri": uri,
            "label": make_entry_label(uri, i),
            "entry": entry,
            "profile_ids": profile_ids,
            "ls": linestyles[(i - 1) % len(linestyles)],
        })

    if not datasets:
        raise RuntimeError("No usable datasets loaded.")

    # Requested times list
    if args.times and len(args.times) > 0:
        t_req = np.array([float(x) for x in args.times], dtype=float)
    else:
        t_req = compute_common_times(time_arrays, args.ncurves)

    # Rainbow colors for time
    cmap = plt.cm.rainbow
    colors = cmap(np.linspace(0.0, 1.0, len(t_req)))

    # Plot
    fig, ax = plt.subplots()
    fig.subplots_adjust(right=0.70)  # space for legends outside
    ax.grid(True)

    xlab_used = None
    ylab_used = None

    for d in datasets:
        entry = d["entry"]
        profile_ids = d["profile_ids"]
        ls = d["ls"]

        for tj, col in zip(t_req, colors):
            # Closest time slice via get_slice
            ids_slice = entry.get_slice(profile_ids, float(tj), imas.ids_defs.CLOSEST_INTERP, 0, autoconvert=False)
            ids_slice = imas.convert_ids(ids_slice, entry.dd_version)

            # Warn if selected time far from requested
            t_sel = slice_time_value(ids_slice)
            if t_sel is not None:
                dt = abs(t_sel - float(tj))
                if dt > float(args.warn_dt):
                    print(f"[WARN] {d['label']}: requested {tj:.6f}s, selected {t_sel:.6f}s (|dt|={dt:.3e}s)")

            # x coordinate with fallbacks
            warn_prefix = d["label"]
            if d["i"] in warned_xcoord:
                warn_prefix = ""  # warn only once per entry
            x, xlab, used_attr = choose_xcoord_from_slice(ids_slice, args.xcoord, warn_prefix=warn_prefix)
            if used_attr != args.xcoord:
                warned_xcoord.add(d["i"])

            # y quantity
            y, ylab = get_y_from_slice(ids_slice, args.quantity)

            xlab_used = xlab
            ylab_used = ylab

            ax.plot(x, y, linestyle=ls, color=col)

    ax.set_xlabel(xlab_used if xlab_used else args.xcoord)
    ax.set_ylabel(ylab_used if ylab_used else args.quantity)

    # Legends outside
    time_handles = [
        Line2D([0], [0], color=col, linestyle="-", linewidth=2, label=f"{tj*1e3:.1f} ms")
        for tj, col in zip(t_req, colors)
    ]
    leg_time = ax.legend(
        handles=time_handles,
        title="Requested time",
        fontsize=8,
        loc="upper left",
        bbox_to_anchor=(1.02, 1.00),
        borderaxespad=0.0,
    )
    ax.add_artist(leg_time)

    entry_handles = [
        Line2D([0], [0], color="black", linestyle=d["ls"], linewidth=2, label=d["label"])
        for d in datasets
    ]
    ax.legend(
        handles=entry_handles,
        title="Entry (pulse/run)",
        fontsize=8,
        loc="upper left",
        bbox_to_anchor=(1.02, 0.45),
        borderaxespad=0.0,
    )

    fig.tight_layout()
    plt.show()

    # Close entries
    for d in datasets:
        try:
            d["entry"].close()
        except Exception:
            pass


if __name__ == "__main__":
    main()
