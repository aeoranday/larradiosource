from .random import BaseRNGModel

import numpy as np
import numpy.typing as npt
import pydantic_numpy.typing as pnpt
from pydantic import BaseModel, field_validator

from abc import ABC, abstractmethod
from typing import Any, Literal

type Vector2 = tuple[float, float]
type Vector3 = tuple[float, float, float]


class BaseGeometry(ABC, BaseRNGModel):
    """
    The abstract base geometry class.
    """

    origin: pnpt.Np1DArray
    height: float
    type: str

    def is_inside(self, position: Vector3) -> bool:
        """Check if the given position is inside the geometry."""
        if not self._check_xy_position(position[:2]):
            return False

        if not self._check_z_position(position[-1]):
            return False
        return True

    @abstractmethod
    def _check_xy_position(self, position: Vector2) -> bool:
        """Check that the given position is inside the geometry's x-y plane."""
        pass

    @abstractmethod
    def _check_z_position(self, z: float) -> bool:
        """Check that the given z position is inside the geometry's z range."""
        pass

    @abstractmethod
    def get_random_face_position(self) -> pnpt.Np1DArray:
        """
        Get a random face position at the base of the geometry.
        """
        pass


class Cylinder(BaseGeometry):
    """
    A cylindrical geometry.
    """

    type: Literal["cylinder"]
    radius: float
    _radius2: float = 0

    def model_post_init(self, context: Any) -> None:
        self._radius2 = self.radius**2
        return

    def _check_xy_position(self, position: Vector2) -> bool:
        pos: npt.NDArray[float] = np.asarray(position) - self.origin[:2]
        return np.sum(pos**2) < self._radius2

    def _check_z_position(self, z: float) -> bool:
        offset_z: float = z - self.origin[-1]
        return 0 <= offset_z < self.origin[-1] + self.height

    def get_random_face_position(self) -> pnpt.Np1DArray:
        """Get a random face position at the base of the cylinder."""
        rng: np.random.Generator = self.get_random_generator()
        theta: float = rng.random() * 2 * np.pi
        radius: float = np.sqrt(rng.random()) * self.radius
        return np.asarray([radius * np.cos(theta), radius * np.sin(theta), 0]) + self.origin


class RectangularPrism(BaseGeometry):
    """
    A rectangular prism geometry.
    """

    type: Literal["rectangular_prism"]
    face_size: pnpt.Np1DArray

    @field_validator("face_size", mode="before")
    @classmethod
    def validate_face_size(cls, vec: pnpt.Np1DArray) -> pnpt.Np1DArray:
        size: npt.NDArray[float] = np.asarray(vec)
        if np.any(size < 0):
            raise ValueError(f"RectangularPrism face_size must be non-negative; received {size}")
        if (len_size := len(size)) != 2:
            raise ValueError(f"RectangularPrism face_size is only 2D; received {len_size}D")
        return size

    def _check_xy_position(self, position: pnpt.Np1DArray) -> bool:
        pos: npt.NDArray[float] = np.asarray(position) - self.origin[:2]

        return not (np.any(pos < 0) or np.any(self.face_size < pos))

    def _check_z_position(self, z: float) -> bool:
        offset_z: float = z - self.origin[-1]
        return 0 <= offset_z < self.origin[-1] + self.height

    def get_random_face_position(self) -> pnpt.Np1DArray:
        """Get a random face position at the base of the rectangular prism."""
        rng: np.random.Generator = self.get_random_generator()
        rands: npt.NDArray[float] = rng.random(2) * self.face_size
        return np.asarray([rands[0], rands[1], self.origin[-1]])
