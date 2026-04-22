from .geometry import Cylinder, RectangularPrism
from .random import BaseRNGModel

import numpy as np
import numpy.typing as npt
from pydantic import BaseModel, Field, field_validator

from typing import Literal, Self


class Emission(BaseRNGModel):
    type: Literal["photon", "electron"]
    energy_mev: float
    interaction_dist: float = 0
    internal_conversion_coefficients: dict[str, float]
    is_corrected: bool = False

    @field_validator("internal_conversion_coefficients")
    @classmethod
    def validate_internal_conversion_coefficients(cls, internal_conversion_coefficients) -> dict[str, float]:
        """ Check the keys are reasonable shells names. """
        for shell_str in internal_conversion_coefficients:
            if shell_str.upper() not in ("K", "L", "M", "N"):
                raise KeyError(f"Internal conversion shell string {shell_str} is not implemented. Use one of [K, L, M, N].")
        return internal_conversion_coefficients

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

    def get_corrected_emission(self, shell_to_binding_energy: dict[str, float]) -> Emission:
        """
        Get the energy-corrected emission from the known binding energies.

        Parameter:
            shell_to_binding_energy (dict[str, float]):
                A map of the shells and binding energies to use for potential corrections.

        Returns a copy of this emission with corrected energy.
        """
        if self.is_corrected:
            return self.model_copy()

        rng: np.random.Generator = self.get_random_generator()

        rand: float = rng.random()
        total_prob: float = 0.0
        for shell_str, ic_prob in self.internal_conversion_coefficients.items():
            total_prob += ic_prob
            if rand < total_prob:
                try:
                    corrected_energy: float = self.energy_mev - shell_to_binding_energy[shell_str]
                except KeyError:
                    raise KeyError("Mismatch on Emission's IC coefficients and shell binding energies.\n"
                                   f"Missing {shell_str} in the shell binding energies.")
                return Emission.model_validate(dict(type="electron",
                                energy_mev=corrected_energy,
                                interaction_dist=0,
                                internal_conversion_coefficients={},
                                is_corrected=True,
                                ))

        # If we did not convert to an electron, then it stayed a photon and should not be corrected again.
        return self.model_copy(update=dict(is_corrected=True))


class DecayBranch(BaseModel):
    probability: float = Field(gt=0, le=1)
    emissions: tuple[Emission, ...]

    def get_corrected_emissions(self, shell_to_binding_energy: dict[str, float]) -> list[Emission, ...]:
        """
        Get the corrected emissions for this decay branch.

        Parameter:
            shell_to_binding_energy (dict[str, float]):
                The shell to binding energy map to perform potential energy corrections with.

        Returns list of Emissions.
        """
        corrected_emissions: list[Emission, ...] = [emission.get_corrected_emission(shell_to_binding_energy) for emission in self.emissions]
        return corrected_emissions


class Source(BaseRNGModel):
    decay_branches: tuple[DecayBranch, ...]
    decay_rate: float  # Bq
    daughter_shell_to_binding_energy: dict[str, float]
    geometry: Cylinder | RectangularPrism = Field(discriminator="type")

    @field_validator("decay_branches")
    @classmethod
    def validate_decay_branches(cls, decay_branches: tuple[DecayBranch, ...]) -> tuple[DecayBranch, ...]:
        prob_total: float = 0
        for branch in decay_branches:
            prob_total += branch.probability

        if np.abs(prob_total - 1) > 0.01:
            raise ValueError(f"Total decay branch probabilities is {prob_total}. This is outside a tolerance of 1 +- 0.01.")

        return decay_branches

    @field_validator("daughter_shell_to_binding_energy")
    @classmethod
    def validate_daughter_shell_to_binding_energy(cls, daughter_shell_to_binding_energy: dict[str, float]) -> dict[str, float]:
        """ Check that the keys are sensible to atomic shells. """
        for shell_str in daughter_shell_to_binding_energy.keys():
            if shell_str.upper() not in ("K", "L", "M", "N"):
                raise KeyError(f"Daughter shell string {shell_str} is not implemented. Use one of [K, L, M, N].")
        return daughter_shell_to_binding_energy

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

    def get_random_decay_emissions(self) -> list[Emission, ...]:
        """
        Get a random decay branch and return the corrected emissions.
        """
        branch: DecayBranch = self._get_random_decay_branch()
        emissions: list[Emission, ...] = branch.get_corrected_emissions(self.daughter_shell_to_binding_energy)
        return emissions

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
