from bi207_simulation.detector import Event

import numpy as np
import numpy.typing as npt

from abc import ABC, abstractmethod
from dataclasses import dataclass


@dataclass
class Waveform:
    waveform: npt.NDArray[int]
    peak: int
    integral: int
    channel: int


class SignalProcessor(ABC):
    window_width: int

    def __init__(self, window_width: int = 0) -> None:
        self.window_width = window_width
        return

    def process_event(self, event: Event) -> list[Waveform]:
        """
        Process the given Event and return a list of Waveforms.
        """
        waveforms: list[Waveform] = []
        for channel, amplitude in zip(event.channels, event.amplitudes):
            waveform: npt.NDArray[int] = self._process_channel_amplitude(amplitude)
            wf: Waveform = Waveform(waveform, waveform.max(), waveform.sum(), channel)
            waveforms.append(wf)
        return waveforms

    @abstractmethod
    def _process_channel_amplitude(self, amplitude: float) -> npt.NDArray[int]:
        """Signal shaping that is dependant on the electronic response function."""
        pass


class PurityMonitorSignalProcessor(SignalProcessor):
    convolution_function: npt.NDArray[float]

    def __init__(self, window_width: int) -> None:
        super().__init__(window_width)
        x: npt.NDArray[float] = np.arange(1, window_width + 1) * 0.1
        self.convolution_function = np.exp(-(x - 10) / 5) / (0.56988 * (1 + np.exp(-(x - 10) / 1.25)))
        return

    def _process_channel_amplitude(self, amplitude: float) -> npt.NDArray[int]:
        """Signal shaping that models the electronic response function of the NP02 Purity Monitor."""
        return self.convolution_function * amplitude
