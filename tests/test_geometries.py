#!/usr/bin/env python3

from bi207_simulation.geometry.cylinder import Cylinder
from bi207_simulation.geometry.rectangular_prism import RectangularPrism

import numpy as np
import numpy.typing as npt


def test_containment():
    test_in_position: npt.NDArray[float] = np.array([0.5, 0.5, 0.5])
    test_out_position: npt.NDArray[float] = np.array([-1, -1, -1])

    origin: tuple[float, float, float] = (0, 0, 0)

    radius: float = 2
    height: float = 5
    cylinder: Cylinder = Cylinder(origin=origin,
                                  radius=radius,
                                  height=height)
    length_x: float = 5
    width_y: float = 4
    height_z: float = 3
    rectangle: RectangularPrism = RectangularPrism(origin=origin,
                                                   length_x=length_x,
                                                   width_y=width_y,
                                                   height_z=height_z,
                                                   )

    if not cylinder.is_position_contained(test_in_position):
        raise RuntimeError("Failed cylinder in-position containment.")
    if not rectangle.is_position_contained(test_in_position):
        raise RuntimeError("Failed rectangular prism position containment.")

    if cylinder.is_position_contained(test_out_position):
        raise RuntimeError("Failed cylinder out-position containment.")
    if rectangle.is_position_contained(test_out_position):
        raise RuntimeError("Failed rectangular prism out-position containment.")
    return


def main():
    print("Testing geometry containments.")
    test_containment()
    print("All geometries passed containment checks.")
    return


if __name__ == "__main__":
    main()
