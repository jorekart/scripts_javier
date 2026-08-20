import h5py
import numpy as np
import matplotlib.pyplot as plt

import h5py
import numpy as np
import matplotlib.pyplot as plt


def plot_solution(
    filename,
    selected_times=None,
    n_slices=5,
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
    """

    with h5py.File(filename, "r") as f:

        xvals = f["xvals"][:]
        times = f["time"][:]
        varvals = f["varvals"]

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
                label=f"t = {times[it]} s"
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
