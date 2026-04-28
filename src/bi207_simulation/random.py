import numpy as np
from pydantic import BaseModel

from typing import Self


class BaseRNGModel(BaseModel):
    """A BaseModel with additional RNG features."""

    _rng: np.random.Generator | None = None

    def set_random_generator(self, rng: np.random.Generator) -> Self:
        """Set a NumPy random number Generator."""
        self._rng = rng
        return self

    def get_random_generator(self) -> np.random.Generator:
        """Get this class's random generator."""
        if self._rng is None:
            # Set the generator if it wasn't already set.
            self.set_random_generator(np.random.default_rng())
        return self._rng
