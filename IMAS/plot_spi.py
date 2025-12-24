import argparse
import sys
import imas
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

sys.path.append('/home/ITER/artolaj/scripts_javier/IMAS')
from imas_custom_utils2 import read_uris, get_entries, parse_imas_uri_kv, make_entry_label


# Load vessel and first wall
vv = np.loadtxt('/home/ITER/artolaj/ITER_components/ITER_vessel_inner.txt')
fw = np.loadtxt('/home/ITER/artolaj/ITER_components/1st_wall.txt')

COLORS = ['red', 'blue', 'green', 'purple', 'orange']


def make_fig_title(uri: str, idx: int) -> str:
    """
    Similar spirit to make_entry_label but good for suptitle.
    """
    kv = parse_imas_uri_kv(uri)
    pulse = kv.get("pulse", "?")
    run = kv.get("run", "?")
    db = kv.get("database", kv.get("db_name", "?"))
    user = kv.get("user", kv.get("user_name", "?"))
    return f"{idx}: pulse={pulse} run={run} db={db} user={user}"


def data_to_display_sizes(radius, ax):
    # same as your original logic
    display_sizes = []
    for r in radius:
        size_in_data_coords = ax.transData.transform((r, 0)) - ax.transData.transform((0, 0))
        size_in_points = size_in_data_coords[0] * 3e1
        display_sizes.append(size_in_points)
    return display_sizes


def plot_spi_entry(entry, uri: str, idx: int):
    """
    Make one SPI slider figure for one entry.
    Returns (fig, slider) so objects stay referenced.
    """
    # Read SPI with autoconvert=False then convert_ids to entry.dd_version
    spi = entry.get('spi', autoconvert=False)
    spi = imas.convert_ids(spi, entry.dd_version)

    # Determine time base length (as in your original script)
    ntime = len(spi.injector[0].fragment[0].position.r)
    it0 = 0

    fig, ax = plt.subplots()
    plt.subplots_adjust(bottom=0.2)

    ax.plot(vv[:, 0] * 1e-3, vv[:, 1] * 1e-3, label='Vessel')
    ax.plot(fw[:, 0], fw[:, 1], label='First Wall')

    ax.set_aspect('equal', adjustable='box')
    ax.set_xlabel('R [m]')
    ax.set_ylabel('Z [m]')

    fragment_scatters = []
    for j, _injector in enumerate(spi.injector):
        sc = ax.scatter([], [], s=[],
                        label=f'SPI Injector {j + 1}',
                        color=COLORS[j % len(COLORS)])
        fragment_scatters.append(sc)

    ax.legend(fontsize=8)
    fig.suptitle(make_fig_title(uri, idx), fontsize=10)

    def update_plot(_val):
        it = int(slider.val)

        for j, injector in enumerate(spi.injector):
            frag_R = []
            frag_Z = []
            frag_radius = []

            for frag in injector.fragment:
                frag_R.append(frag.position.r[it])
                frag_Z.append(frag.position.z[it])
                frag_rad = (frag.volume[it] * 3 / (4 * np.pi)) ** (1 / 3.0)
                frag_radius.append(frag_rad)

            sizes = data_to_display_sizes(frag_radius, ax)
            fragment_scatters[j].set_offsets(np.c_[frag_R, frag_Z])
            fragment_scatters[j].set_sizes(sizes)

        ax.set_title(f"Fragments at Time = {spi.time[it] * 1e3:.1f} ms")
        fig.canvas.draw_idle()

    ax_slider = plt.axes([0.2, 0.05, 0.6, 0.03], facecolor='lightgoldenrodyellow')
    slider = Slider(ax_slider, 'Time Index', 0, ntime - 1, valinit=it0, valstep=1)
    slider.on_changed(update_plot)

    # Initial draw
    update_plot(it0)

    return fig, slider


def main():
    parser = argparse.ArgumentParser(description="SPI fragments viewer for one or many URIs.")
    src = parser.add_mutually_exclusive_group(required=True)
    src.add_argument("--uri", "-uri", type=str,
                     help="Single IMAS URI e.g. 'imas:hdf5?user=...;database=...;version=...;pulse=...;run=...'")
    src.add_argument("--uri_file", "-uf", type=str, help="Text file with one IMAS URI per line.")

    args = parser.parse_args()

    uris = read_uris(uri=args.uri, uri_file=args.uri_file)

    # mirror your other tool call style
    entries = get_entries(uris, mode="r", dd_version="None")  # or set dd_version="4.0.0" if desired

    figs_and_sliders = []
    for i, (entry, uri) in enumerate(zip(entries, uris), start=1):
        print(f" Getting uri={uri} with version = {entry.dd_version}")
        try:
            fig, slider = plot_spi_entry(entry, uri, i)
            figs_and_sliders.append((fig, slider))  # keep references alive
            plt.show()
        except Exception as e:
            print(f"[WARN] Skipping {make_entry_label(uri, i)} بسبب error:\n  {e}\n")

    if not figs_and_sliders:
        raise RuntimeError("No SPI figures created (all entries failed).")

    


if __name__ == "__main__":
    main()
