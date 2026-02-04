#!/bin/python3

import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt

from pathlib import Path


DATA_PATH: Path = Path("~/Code/Bi-207 Simulation/bi207_simulation/bi207stream.txt").expanduser()


def main() -> int:
    data: npt.NDArray[float] = np.loadtxt(DATA_PATH)
    data = data[data[:, 1] == 0]  # Only get the long TPC
    print(data.shape)
    #data = data[data[:, -1] == 0]  # Only get branch A

    photons = data[data[:, 3] == 0]
    electrons = data[data[:, 3] == 1]

    photon_energy = photons[:, -2]
    electron_energy = electrons[:, -2]

    bins: npt.NDArray = np.arange(0.26, 2, 0.02)

    plt.figure(figsize=(6, 4), dpi=300, layout="constrained")

    plt.hist([photon_energy, electron_energy],
             bins=bins,
             edgecolor=["#EE442F", "#63ACBE"],
             fill=False,
             stacked=True,
             label=["Photon", "Electron"],
             )

    plt.title("Bi-207 Simulated Energy Emissions (Target)")
    plt.xlabel("Energy (MeV)")
    plt.ylabel("Count")

    plt.box(False)
    plt.legend(frameon=False)

    plt.savefig("energy_target.png")

    plt.show()
    plt.close()
    return 0


if __name__ == "__main__":
    main()
