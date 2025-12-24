import argparse
import sys
import imas
import matplotlib.pyplot as plt

sys.path.append('/home/ITER/artolaj/scripts_hub/scripts_javier/IMAS')
from imas_custom_utils2 import read_uris, get_entries, parse_imas_uri_kv, make_entry_label


def get_time_value_and_label(entry, qtty):
    if qtty == "I_RE":
        runaway = entry.get("runaway_electrons", autoconvert=False)
        runaway = imas.convert_ids(runaway, entry.dd_version)
        return runaway.time * 1e3, runaway.global_quantities.current_phi * 1e-6, "I_RE [MA]"

    summary = entry.get("summary", autoconvert=False)
    summary = imas.convert_ids(summary, entry.dd_version)
    time_ms = summary.time * 1e3

    if qtty == "Ip":
        return time_ms, summary.global_quantities.ip.value * 1e-6, "Ip [MA]"
    if qtty == "Zaxis":
        return time_ms, summary.local.magnetic_axis.position.z, "Zaxis [m]"
    if qtty == "li_3":
        return time_ms, summary.global_quantities.li_3.value, "li_3 [-]"
    if qtty == "q95":
        return time_ms, summary.global_quantities.q_95.value, "q95 [-]"
    if qtty == "Wth":
        return time_ms, summary.global_quantities.energy_thermal.value, "W_th [J]"
    if qtty == "Te_av":
        return time_ms, summary.volume_average.t_e.value, "Te_av [eV]"
    if qtty == "ne_av":
        return time_ms, summary.volume_average.n_e.value, "ne_av [eV]"

    raise ValueError(f"Unknown quantity: {qtty}")


def main():
    parser = argparse.ArgumentParser(description="Plot IMAS quantity for one or many URIs.")
    src = parser.add_mutually_exclusive_group(required=True)
    src.add_argument("--uri", "-uri", type=str, help='Single IMAS URI (wrap in quotes).')
    src.add_argument("--uri_file", "-uf", type=str, help="Text file with one IMAS URI per line.")
    parser.add_argument("-q", "--quantity", default="Ip",
                        help="Ip/Zaxis/li_3/q95/Wth/Te_av/ne_av/I_RE")

    args = parser.parse_args()

    uris = read_uris(uri=args.uri, uri_file=args.uri_file)
    entries = get_entries(uris, mode="r", dd_version="None") # Put dd_version to the version you want 
    #entries = get_entries(uris, mode="r", dd_version="4.0.0") # Works only between same major versions

    ylab_used = None
    for i, (entry, uri) in enumerate(zip(entries, uris), start=1):
        print(f' Getting uri={uri} with version = {entry.dd_version}')
        t_ms, val, ylab = get_time_value_and_label(entry, args.quantity)
        ylab_used = ylab
        plt.plot(t_ms, val, "*-", label=make_entry_label(uri, i))

    plt.xlabel("Time [ms]")
    plt.ylabel(ylab_used if ylab_used else args.quantity)
    plt.grid(True)
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
