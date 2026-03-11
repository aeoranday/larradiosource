from .geometry import Cylinder, RectangularPrism
from .radiation import Emission

import numpy as np
import numpy.typing as npt
from pydantic import BaseModel, Field, field_validator

from typing import Any, Literal


class Event(BaseModel):
    main_channel_index: int = 0
    channels: tuple[int, ...] = []
    amplitudes: tuple[float, ...] = []
    type: Literal["electron", "photon", None] = None

    def model_post_init(self, ctx: Any) -> None:
        if (num_chan := len(self.channels)) != (num_amps := len(self.amplitudes)):
            raise ValueError(f"The number of channels and amplitudes in this Event are not equal. len(channels) == {num_chan} and len(amplitudes) == {num_amps}")
        return

    def __len__(self) -> int:
        return len(self.channels)


class Detector(BaseModel):
    geometry: Cylinder | RectangularPrism = Field(discriminator="type")
    num_channels: int
    channel_widths: tuple[float, ...]
    channel_spacing: float
    channel_start_pos: float
    channel_increment_dir: Literal["x", "y", "r"]
    threshold: float = 0.26
    drift_velocity: float = 1.6  # mm / us
    electron_lifetime_ms: float = 3  # ms
    channel_resolution_mev: float = 0.08
    signal_width: int = 12
    fake_scale: float = 0.8
    _Rc: float = 0.6980  # https://lar.bnl.gov/properties/ : R_c (Birks Model)
    _fC_per_MeV: float = 6.79  # fC / MeV

    def model_post_init(self, ctx: Any) -> None:
        self._attenuation_dist_mm: float = self.electron_lifetime_ms * 1000 * self.drift_velocity

        num_channel_widths: int = len(self.channel_widths)
        if num_channel_widths> self.num_channels:
            raise ValueError(f"There are more channel widths {num_channel_widths} than there are channels {self.num_channels}.")
        if num_channel_widths < self.num_channels and num_channel_widths != 1:
            raise ValueError(f"There are less channel widths {num_channel_widths} than there are channels {self.num_channels}.")
        channel_widths: npt.NDArray[float] = np.asarray(self.channel_widths)
        if len(self.channel_widths) == 1:
            channel_widths = np.asarray((self.channel_widths[0],) * self.num_channels)

        self._channel_center_positions: npt.NDArray[float] = np.arange(self.channel_start_pos,
                                                                       (self.num_channels + 1) * self.channel_spacing,
                                                                       self.channel_spacing)
        self._channel_coverage: npt.NDArray[float] = np.zeros((self.num_channels, 2))
        for channel in range(self.num_channels):
            center: float = self._channel_center_positions[channel]
            self._channel_coverage[channel] = (center - channel_widths[channel],
                                               center + channel_widths[channel])

        self._channel_resolution_fC: float = self.channel_resolution_mev * self._fC_per_MeV
        return

    def _reduce_electron_energy_to_charge(self, energy_mev: float, dist: float) -> float:
        """
        Apply recombination and attachment reduction factors for electron transport.

        Makes use of the starting energy and travel distance and converts from energy to charge.

        Parameters:
            energy_mev (float):
                Initial energy in MeV of the electron.
            dist (float):
                Distance to travel in mm.

        Returns:
            (float):
                Resultant charge after transportation.
        """
        return self._Rc * np.exp(dist / self._attenuation_dist_mm) * energy_mev * self._fC_per_MeV * self.fake_scale

    def process_emission(self, emission: Emission, position: npt.NDArray[float]) -> Event:
        """
        Process the given emission to produce an event.
        """
        if not self.geometry.is_inside(position):
            return Event(channels=[], amplitudes=[])

        r: float = 0
        match self.channel_increment_dir:
            case 'x':
                r = position[0]
            case 'y':
                r = position[1]
            case 'r':
                r = np.linalg.norm(position[:2])
            case _:
                raise ValueError(f"Unsupported channel direction: {self.channel_increment_dir}.")


        energy: float = emission.energy_mev
        if emission.type == "photon":
            energy = emission.get_compton_energy()

        dist: float = position[2] - (self.geometry.origin[2] + self.geometry.height)
        charge: float = self._reduce_electron_energy_to_charge(energy, dist)

        main_channel_index: int = 0
        idx: int = 0
        min_channel_dist: float = np.inf
        channels: list[int] = []
        amplitudes: list[float] = []
        for channel in range(self.num_channels):
            coverage: npt.NDArray[float] = self._channel_coverage[channel]
            center: float = self._channel_center_positions[channel]
            if r >= coverage[0] and r < coverage[1]:
                rand_noise: float = self._channel_resolution_fC * (np.sum(np.random.rand(self.signal_width)) - self.signal_width / 2)
                response: float = charge + rand_noise
                if response < self.threshold:
                    continue
                if (new_min := np.abs(r - center)) < min_channel_dist:
                    min_channel_dist = new_min
                    main_channel_index = idx
                idx += 1
                channels.append(channel)
                amplitudes.append(response)

        return Event(main_channel_index=main_channel_index, channels=channels, amplitudes=amplitudes, type=emission.type)
