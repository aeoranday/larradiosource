import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt

from pathlib import Path


DATA_PATH: Path = Path("~/Code/Bi-207 Simulation/bi207_simulation/bi207spectra.txt").expanduser()


def main() -> int:
    data: npt.NDArray[float] = np.loadtxt(DATA_PATH, delimiter=" ")
    energies: npt.NDArray[float] = data[:, 0]
    spectra_out: npt.NDArray[float] = data[:, 1]
    spectra_in: npt.NDArray[float] = data[:, 4]

    edges = np.arange(0, 401, 1.0)
    counts = spectra_out

    plt.figure()
    plt.box(False)

    plt.stairs(counts, edges, fill=False, color='k')

    plt.xlabel("Energy (MeV)")
    plt.ylabel("Count")

    plt.show()
    plt.close()
    return 0


if __name__ == "__main__":
    main()
