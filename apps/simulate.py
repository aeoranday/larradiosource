from bi207_simulation.detector import Detector, Event
from bi207_simulation.radiation import Source, DecayBranch, Emission
from bi207_simulation.signals import PurityMonitorSignalProcessor, Waveform

import click
import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt

from collections import defaultdict
from pathlib import Path
from pprint import pprint
from tomllib import load
from typing import Any


def plot_waveforms(waveforms: list[Waveform]) -> None:
    peaks: npt.NDArray[float] = np.asarray([waveform.peak for waveform in waveforms]) * 200
    plt.figure()
    plt.box(False)

    plt.hist(peaks, bins=50, fill=False, edgecolor='k')

    plt.xlabel("Peak ADC")
    plt.ylabel("Count")

    plt.show()
    plt.close()
    return


def plot_amplitude_histogram(events: list[Event]) -> None:
    amplitudes: list[list[float, ...], list[float, ...]] = [[], []]

    for event in events:
        if event.type == "electron":
            amplitudes[1].append(event.amplitudes[0])
        else:
            amplitudes[0].append(event.amplitudes[0])

    bins: npt.NDArray[float] = np.arange(0.26, 2, 0.02)
    plt.figure(figsize=(6, 4), dpi=300, layout="constrained")
    plt.hist(amplitudes,
             bins=bins,
             edgecolor=["#EE442F", "#63ACBE"],
             stacked=True,
             fill=False,
             label=["Photon", "Electron"],
             )
    plt.legend(frameon=False)

    plt.title("Bi-207 Simulated Energy Emissions")
    plt.xlabel("Energy (MeV)")
    plt.ylabel("Count")
    plt.box(False)

    plt.show()
    plt.close()
    return


@click.command()
@click.argument("config_path", type=click.Path(readable=True, resolve_path=True, path_type=Path))
def main(config_path: Path) -> int:
    np.random.seed(42)
    with open(config_path, 'rb') as f:
        config_dict: dict[str, Any] = load(f)

    num_events: int = config_dict.get("num_events", 2)

    if "source_A" not in config_dict:
        raise KeyError(f"Configuration is missing the `source` config section.")
    if "source_B" not in config_dict:
        raise KeyError(f"Configuration is missing the `source` config section.")
    sources: list[Source] = []
    sources.append(Source.model_validate(config_dict['source_A']))
    sources.append(Source.model_validate(config_dict['source_B']))

    if "detector_short" not in config_dict:
        raise KeyError(f"Configuration is missing the `detector_short` config section.")
    detector: Detector = Detector.model_validate(config_dict['detector_short'])

    events: list[Event] = []

    event_count: int = 0
    while event_count < num_events:
        had_event: bool = False
        for source in sources:
            emission: Emission = source.get_emission()
            position: npt.NDArray[float] = source.get_emission_position(emission)
            if not detector.geometry.is_inside(position):
                continue

            event: Event = detector.process_emission(emission, position)
            if len(event) > 0:
                events.append(event)
                had_event = True
        if had_event:
            event_count += 1

    plot_amplitude_histogram(events)

    signal_processor: PurityMonitorSignalProcessor = PurityMonitorSignalProcessor(500)
    waveforms: list[Waveform] = []
    for event in events:
        waveforms += signal_processor.process_event(event)

    plot_waveforms(waveforms)
    return 0


if __name__ == "__main__":
    main()
