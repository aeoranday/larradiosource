from bi207_simulation.radiation import Source, Emission

import numpy as np
import numpy.typing as npt


POSITION_COUNT: int = 10


def expected_distribution() -> npt.NDArray[float]:
    np.random.seed(42)
    positions: npt.NDArray[float] = np.zeros((POSITION_COUNT, 3))
    for idx in range(POSITION_COUNT):
        while True:
            x = 5 * np.random.rand() - 2.5
            y = 5 * np.random.rand() - 2.5
            if x**2 + y**2 < 6.25:
                break
        print(x,y)

        d = -140 * np.log(np.random.rand())
        print("d", d)
        cth = np.random.rand()
        print("cth", cth)
        phi = np.pi * 2 * np.random.rand()
        print("phi", phi)

        sth = np.sqrt(1 - cth**2)

        x1 = x + d * sth * np.cos(phi)
        y1 = y + d * sth * np.sin(phi)
        z1 = d * cth
        positions[idx] = [x1, y1, z1]

    return positions


def actual_distribution() -> npt.NDArray[float]:
    positions: npt.NDArray[float] = np.zeros((POSITION_COUNT, 3))
    source: Source = Source(
                geometry={'origin': [0,0,0], 'height': 0, 'type': 'cylinder', 'radius': 2.5},
                decay_branches=[{'probability': 1, "emission": {'type': "photon", "energy_mev": 1, "interaction_dist": 140}}],
                decay_rate=1,
    )
    np.random.seed(42)
    emission: Emission = Emission(type="photon", energy_mev=1, interaction_dist=140)
    for idx in range(POSITION_COUNT):
        position: npt.NDArray[float] = source.get_emission_position(emission)
        positions[idx] = position
    return positions


def main() -> int:
    expected_positions: npt.NDArray[float] = expected_distribution()
    actual_positions: npt.NDArray[float] = actual_distribution()

#    print(expected_positions)
#    print(actual_positions)

    print(np.all(np.isclose(expected_positions, actual_positions)))
    return  0


if __name__ == "__main__":
    main()
