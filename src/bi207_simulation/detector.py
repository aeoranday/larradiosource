from .geometry import Cylinder, RectangularPrism
from .radiation import Emission

import numpy as np
import numpy.typing as npt
from pydantic import BaseModel, Field, field_validator

from typing import Any, Literal


class Event(BaseModel):
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
                                                                       self.num_channels * self.channel_spacing,
                                                                       self.channel_spacing)
        self._channel_coverage: npt.NDArray[float] = np.zeros((self.num_channels, 2))
        for channel in range(self.num_channels):
            center: float = self._channel_center_positions[channel]
            self._channel_coverage[channel] = (center - channel_widths[channel],
                                               center + channel_widths[channel])
        return


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
        energy *= np.exp(dist / self._attenuation_dist_mm)  # Electron survival

        channels: list[int] = []
        amplitudes: list[float] = []
        for channel in range(self.num_channels):
            coverage: npt.NDArray[float] = self._channel_coverage[channel]
            if r >= coverage[0] and r < coverage[1]:
                rand_noise: float = self.channel_resolution_mev * (np.sum(np.random.rand(self.signal_width)) - self.signal_width / 2)
                response: float = energy + rand_noise
                if response < self.threshold:
                    continue
                channels.append(channel)
                amplitudes.append(energy + rand_noise)

        return Event(channels=channels, amplitudes=amplitudes, type=emission.type)
