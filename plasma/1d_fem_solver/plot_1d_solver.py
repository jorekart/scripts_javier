import h5py
import numpy as np
import matplotlib.pyplot as plt


def plot_solution(
    filename,
    selected_times=None,
    n_slices=5,
    average_filename="x_average_values.dat",
):
    """
    Plot solution profiles from an HDF5 solution file.

    Parameters
    ----------
    filename : str
        HDF5 solution file.

    selected_times : list[float] or None, optional
        Physical times to plot. The closest stored
        timestep is used.

    n_slices : int, optional
        Number of equally spaced snapshots (in index
        space) to plot when selected_times is None.

    average_filename : str, optional
        Output text file for the x-average at every stored time.
    """

    with h5py.File(filename, "r") as f:

        xvals = f["xvals"][:]
        times = f["time"][:]
        varvals = f["varvals"][:]

        # Spatial average for every stored time. Trapezoidal integration also
        # handles a non-uniform x grid correctly.
        x_average = np.trapz(varvals, x=xvals, axis=1) / (xvals[-1] - xvals[0])

        np.savetxt(
            average_filename,
            np.column_stack((times, x_average)),
            header="time[s] x_average",
        )
        print(f"Wrote x-averaged values to {average_filename}")

        plt.figure(figsize=(8, 5))

        if selected_times is not None:

            # User provided physical times
            indices = [
                np.argmin(np.abs(times - t))
                for t in selected_times
            ]

        else:

            # Equally spaced indices
            n_steps = len(times)

            indices = np.linspace(
                0,
                n_steps - 1,
                min(n_slices, n_steps),
                dtype=int,
            )

        for it in indices:

            plt.plot(
                xvals,
                varvals[it, :],
                label=f"t = {times[it]*1e3} ms"
            )

        plt.xlabel("x")
        plt.ylabel("Solution")
        plt.title("Solution evolution")
        plt.grid(True)
        plt.legend()
        plt.tight_layout()
        plt.show()


plot_solution(
    "solution.h5",
    n_slices=6
)
