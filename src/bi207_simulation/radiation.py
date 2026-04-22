from .geometry import Cylinder, RectangularPrism
from .random import BaseRNGModel

import numpy as np
import numpy.typing as npt
from pydantic import BaseModel, Field, field_validator

from typing import Literal


class Emission(BaseRNGModel):
    type: Literal["photon", "electron"]
    energy_mev: float
    interaction_dist: float = 0

    def _get_scatter_probability(self) -> float:
        """
        Get the Compton scattering probability using the Klein-Nishina formula.
        """
        rng: np.random.Generator = self.get_random_generator()
        ratio: float = self.energy_mev / 0.511  # m_e * c**2: MeV
        cos_theta: float = 2 * rng.random() - 1
        eps: float = 1 / (1 + ratio * (1 - cos_theta))
        return (eps**2 * (eps + 1 / eps - (1 - cos_theta**2)) / 2, eps)

    def get_compton_energy(self) -> float:
        """
        Get the electron energy from Compton scattering using the Klein-Nishina formula.
        """
        rng: np.random.Generator = self.get_random_generator()
        # It would only really make sense to make this call if the type is a photon.
        if self.type == "electron":
            return self.energy_mev

        rand: float = rng.random()
        scatter_prob, eps = self._get_scatter_probability()

        while rand > scatter_prob:
            rand = np.random.rand()
            scatter_prob, eps = self._get_scatter_probability()

        return self.energy_mev * (1 - eps)


class DecayBranch(BaseModel):
    probability: float = Field(gt=0, le=1)
    emission: Emission


class Source(BaseRNGModel):
    decay_branches: tuple[DecayBranch, ...]
    decay_rate: float  # Bq
    geometry: Cylinder | RectangularPrism = Field(discriminator="type")

    @field_validator("decay_branches")
    @classmethod
    def validate_decay_branches(cls, decay_branches: tuple[DecayBranch, ...]) -> tuple[DecayBranch, ...]:
        prob_total: float = 0
        for branch in decay_branches:
            prob_total += branch.probability

        return decay_branches

    def _get_random_decay_branch(self) -> DecayBranch:
        """ Get a random decay branch. """
        rng: np.random.Generator = self.get_random_generator()

        rand: float = rng.random()
        prob_total: float = 0
        for branch in self.decay_branches:
            prob_total += branch.probability
            if rand < prob_total:
                return branch
        return branch

        # Assuming a rounding error, default to returning the last branch.
        return branch.emission

    def get_emission_position(self, emission: Emission) -> npt.NDArray[float]:
        """
        Get the emission position on the x-y plane and depth for photons.
        """
        rng: np.random.Generator = self.get_random_generator()
        face_pos: npt.NDArray[float] = self.geometry.get_random_face_position()
        if emission.interaction_dist == 0:
            return face_pos

        r: float = -np.log(rng.random()) * emission.interaction_dist
        cos_theta: float = 2*rng.random() - 1
        phi = rng.random() * 2 * np.pi
        cos_phi: float = np.cos(phi)
        sin_phi: float = np.sin(phi)

        sin_theta: float = np.sqrt(1 - cos_theta**2)

        vec: npt.NDArray[float] = np.array([r * cos_phi * sin_theta,
                                            r * sin_phi * sin_theta,
                                            r * cos_theta],
                                           dtype=float)

        return face_pos + vec
