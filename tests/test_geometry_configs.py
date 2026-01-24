#!/usr/bin/env python3

from bi207_simulation.config.geometry.registry import get_geometry_config
from bi207_simulation.geometry.cylinder import Cylinder

from typing import Any


def test_configs():
    cyl_config_dict: dict[str, Any] = {
                "geometry": {
                    "type": "cylinder",
                    "origin": [0, 0, 0],
                    "radius": 5,
                    "height": 1,
                }
            }
    cyl_config: CylinderConfig = get_geometry_config("cylinder").model_validate(cyl_config_dict)
    cylinder: Cylinder = cyl_config.build()
    return


def main():
    test_configs()
    return


if __name__ == "__main__":
    main()
