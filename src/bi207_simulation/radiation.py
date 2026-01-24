from .geometry import Cylinder, RectangularPrism

import numpy as np
import numpy.typing as npt
from pydantic import BaseModel, Field, field_validator

from typing import Literal


class Emission(BaseModel):
    type: Literal["photon", "electron"]
    energy_mev: float
    interaction_dist: float = 0

    def _get_scatter_probability(self) -> float:
        """
        Get the Compton scattering probability using the Klein-Nishina formula.
        """
        ratio: float = self.energy_mev / 0.511  # m_e * c**2: MeV
        cos_theta: float = 2 * np.random.rand() - 1
        eps: float = 1 / (1 + ratio * (1 - cos_theta))
        return (eps**2 * (eps + 1 / eps - (1 - cos_theta**2)) / 2, eps)

    def get_compton_energy(self) -> float:
        """
        Get the electron energy from Compton scattering using the Klein-Nishina formula.
        """
        # It would only really make sense to make this call if the type is a photon.
        if self.type == "electron":
            return self.energy_mev

        rand: float = np.random.rand()
        scatter_prob, eps = self._get_scatter_probability()

        while rand <= scatter_prob:
            rand = np.random.rand()
            scatter_prob, eps = self._get_scatter_probability()

        return self.energy_mev * (1 - eps)


class DecayBranch(BaseModel):
    probability: float = Field(gt=0, le=1)
    emission: Emission


class Source(BaseModel):
    decay_branches: tuple[DecayBranch, ...]
    decay_rate: float  # Bq
    geometry: Cylinder | RectangularPrism = Field(discriminator="type")

    @field_validator("decay_branches")
    @classmethod
    def validate_decay_branches(cls, decay_branches: tuple[DecayBranch, ...]) -> tuple[DecayBranch, ...]:
        prob_total: float = 0
        for branch in decay_branches:
            prob_total += branch.probability
#            print(f"Total P.: {prob_total:.4f}. Int. Dist.: {branch.emission.interaction_dist:>5.1f}. Prob.: {branch.probability}. Energy: {branch.emission.energy_mev:.3f}.")

#        print()
#        for branch in decay_branches:
#            branch.probability /= prob_total

        return decay_branches

    def get_emission(self) -> Emission:
        """
        Get a random decay branch to produce an emission with.
        """
        random: float = np.random.rand()
        prob: float = 0
        for branch in self.decay_branches:
            prob += branch.probability
            if random < prob:
                return branch.emission

        # Assuming a rounding error, default to returning the last branch.
        return branch.emission

    def get_emission_position(self, emission: Emission) -> npt.NDArray[float]:
        """
        Get the emission position on the x-y plane and depth for photons.
        """
        face_pos: npt.NDArray[float] = self.geometry.get_random_face_position()
        if emission.interaction_dist == 0:
            return face_pos

#        random: npt.NDArray[float] = np.random.rand(3)

        r: float = -np.log(np.random.rand()) * emission.interaction_dist
#        print("r", r)
        cos_theta: float = np.random.rand()#2*random[1] - 1
#        print("cos_theta", cos_theta)
        phi = np.random.rand() * 2 * np.pi
#        print("phi", phi)
        cos_phi: float = np.cos(phi)#2*random[0] - 1
        sin_phi: float = np.sin(phi)#np.sqrt(1 - cos_phi**2)

        sin_theta: float = np.sqrt(1 - cos_theta**2)

        vec: npt.NDArray[float] = np.array([r * cos_phi * sin_theta,
                                            r * sin_phi * sin_theta,
                                            r * cos_theta],
                                           dtype=float)

#        print("vec", vec)
        return face_pos + vec

