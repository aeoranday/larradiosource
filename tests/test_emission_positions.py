from bi207_simulation.radiation import Source, Emission

import numpy as np
import numpy.typing as npt


POSITION_COUNT: int = 100_000
INTERACTION_DIST: float = 0


def expected_distribution() -> npt.NDArray[float]:
    np.random.seed(42)
    positions: npt.NDArray[float] = np.zeros((POSITION_COUNT, 3))
    for idx in range(POSITION_COUNT):
        while True:
            x = 5 * np.random.rand() - 2.5
            y = 5 * np.random.rand() - 2.5
            if x**2 + y**2 < 6.25:
                break

        d = -INTERACTION_DIST * np.log(np.random.rand())
        cth = np.random.rand()
        phi = np.pi * 2 * np.random.rand()

        sth = np.sqrt(1 - cth**2)

        x1 = x + d * sth * np.cos(phi)
        y1 = y + d * sth * np.sin(phi)
        z1 = d * cth
        positions[idx] = [x1, y1, z1]

    return positions


def actual_distribution() -> npt.NDArray[float]:
    positions: npt.NDArray[float] = np.zeros((POSITION_COUNT, 3))
    emission: Emission = Emission(type="photon", energy_mev=1, interaction_dist=INTERACTION_DIST)
    source: Source = Source(
                geometry={'origin': [0,0,0], 'height': 0, 'type': 'cylinder', 'radius': 2.5},
                decay_branches=[{'probability': 1, "emission": emission}],
                decay_rate=1,
    )
    np.random.seed(42)
    for idx in range(POSITION_COUNT):
        position: npt.NDArray[float] = source.get_emission_position(emission)
        positions[idx] = position
    return positions


def main() -> int:
    expected_positions: npt.NDArray[float] = expected_distribution()
    actual_positions: npt.NDArray[float] = actual_distribution()


    expected_position_mean: npt.NDArray[float] = np.mean(expected_positions, axis=0)
    actual_position_mean: npt.NDArray[float] = np.mean(actual_positions, axis=0)

    expected_position_std: npt.NDArray[float] = np.std(expected_positions, axis=0)
    actual_position_std: npt.NDArray[float] = np.std(actual_positions, axis=0)

    print(f"Expected position: [{expected_position_mean[0]} ± {expected_position_std[0]}, {expected_position_mean[1]} ± {expected_position_std[1]}, {expected_position_mean[2]} ± {expected_position_std[2]}]")
    print(f"  Actual position: [{actual_position_mean[0]} ± {actual_position_std[0]}, {actual_position_mean[1]} ± {actual_position_std[1]}, {actual_position_mean[2]} ± {actual_position_std[2]}]")
    return  0


if __name__ == "__main__":
    main()
