from bi207_simulation.detector import (Detector, Event)
from bi207_simulation.radiation import (Source, DecayBranch, Emission)

import click
import h5py
import numpy as np
import numpy.typing as npt
from tqdm import tqdm

from collections import defaultdict
from collections.abc import Callable
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from tomllib import load
from typing import Any


BIN_STEP: int = 20
ADC_TO_FC: npt.NDArray[float] = np.array([2.209, 0.063]) * 1e-3  # fC/ADC
BINS: ArrayLike = np.arange(0, 4000+BIN_STEP, BIN_STEP) * ADC_TO_FC[0]


def get_primary_channel(events: list[Event]) -> int:
    channels, counts = np.unique([event.channels[event.main_channel_index] for event in events],
                                 return_counts=True)

    return channels[np.argmax(counts)]


def calculate_charge_distribution(events: list[Event], max_count_channel: int) -> npt.NDArray[int]:
    """
    Calculate the output charge distribution.
    """
    amplitudes = {"photon": defaultdict(list), "electron": defaultdict(list)}

    # Process all events and find the amplitudes for the channel and event type
    for event in events:
        amplitude: float = event.amplitudes[event.main_channel_index]
        main_channel: int = event.channels[event.main_channel_index]
        amplitudes[event.type][main_channel].append(amplitude)

    # Calculate the histograms for each of the channels according to the
    # emission type.
    source_electron_counts = None
    source_photon_counts = None
    gamma_electron_counts = None
    gamma_photon_counts = None
    for channel in range(max_count_channel - 11, max_count_channel + 12):
        if channel >= max_count_channel - 6 and channel < max_count_channel + 7:
            if source_electron_counts is None:
                source_electron_counts, _ = np.histogram(amplitudes['electron'][channel], BINS)
                source_photon_counts, _ = np.histogram(amplitudes['photon'][channel], BINS)
                continue
            new_se_counts, _ = np.histogram(amplitudes['electron'][channel], BINS)
            new_sp_counts, _ = np.histogram(amplitudes['photon'][channel], BINS)
            source_electron_counts += new_se_counts
            source_photon_counts += new_sp_counts
            continue
        if gamma_electron_counts is None:
            gamma_electron_counts, _ = np.histogram(amplitudes['electron'][channel], BINS)
            gamma_photon_counts, _ = np.histogram(amplitudes['photon'][channel], BINS)
            continue
        new_se_counts, _ = np.histogram(amplitudes['electron'][channel], BINS)
        new_sp_counts, _ = np.histogram(amplitudes['photon'][channel], BINS)
        gamma_electron_counts += new_se_counts
        gamma_photon_counts += new_sp_counts


    source_counts: npt.NDArray[int] = source_electron_counts + source_photon_counts
    gamma_counts: npt.NDArray[int] = gamma_electron_counts + gamma_photon_counts

    return (source_electron_counts, source_photon_counts), (gamma_electron_counts, gamma_photon_counts)


def simulate_events(num_events: int, worker_id: int, detector: Detector, sources: list[Source]) -> list[Event]:
    """
    Simulate for n events with the given detector and sources.

    Parameters:
        num_events (int):
            The number of events to simulate.
        detector (Detector):
            The detector to read out from.
        sources (list[Source]):
            The sources that will be decaying.

    Returns a list of Events that has a length equal to `num_events`.
    """
    events: list[Event] = []

    event_count: int = 0
    pbar = tqdm(total=num_events, position=worker_id)

    # This condition can overshoot: a decay branch can have multiple emissions and multiple events.
    while event_count < num_events:
        new_event_count: int = 0
        for source in sources:
            emissions: list[Emission, ...] = source.get_random_decay_emissions()
            for emission in emissions:
                position: npt.NDArray[float] = source.get_emission_position(emission)
                if not detector.geometry.is_inside(position):
                    continue

                event: Event = detector.process_emission(emission, position)
                # Emissions may be below threshold and result in an empty Event.
                if len(event.amplitudes) == 0:
                    continue
                events.append(event)
                new_event_count += 1
        event_count += new_event_count
        pbar.update(new_event_count)
    return events


def categorize_events(events: list[Event]) -> dict[str, dict[str, list[float]]]:
    """
    For the list of events, categorize the event type, energy, and resultant charge.
    """
    categories: dict[str, dict[str, list[float]]] = {"photon": defaultdict(list), "electron": defaultdict(list)}
    for event in events:
        emission: Emission = event.init_emission
        subcategory: dict[str, list[float]] = categories[emission.type]
        # Prefer to index as keV and integers.
        subcategory[(int(emission.energy_mev * 1e3), event.channels[event.main_channel_index])].append(event.amplitudes[event.main_channel_index])

    return categories


def run_simulation(processor_count: int,
                   num_events: int,
                   detector: Detector,
                   sources: list[Source],
                   sim_save_path: Path,
                   ) -> None:
    """
    Run the simulation and save the results to an HDF5 file.

    Arguments:
        processor_count (int):
            The number of processors to use.
        num_events (int):
            The number of events to produce.
        detector (Detector):
            The detector that will process emissions.
        sources (list[Source]):
            The list of Sources that will produce emissions.
    """
    num_events_per_processor: list[int] = [num_events // processor_count,] * processor_count
    if (remainder_event_count := num_events % processor_count) != 0:
        num_events_per_processor[-1] = num_events_per_processor[-1] + remainder_event_count

    events: list[Event] = []

    with ProcessPoolExecutor(max_workers=processor_count) as executor:
        futures = [executor.submit(simulate_events, n, idx, detector, sources) for idx, n in enumerate(num_events_per_processor)]
        for future in futures:
            events += future.result()

    max_count_channel: int = get_primary_channel(events)
    source_distrs, gamma_distrs = calculate_charge_distribution(events, max_count_channel)

    categories: dict[str, dict[tuple[int, int], list[float]]] = categorize_events(events)

    sim_file: h5py.File = h5py.File(sim_save_path, 'w')

    sim_file.attrs['Simulation'] = True
    sim_file.attrs['bins'] = BINS
    sim_file.create_dataset("Source", data=source_distrs[0] + source_distrs[1])
    sim_file.create_dataset("Source Electrons", data=source_distrs[0])
    sim_file.create_dataset("Source Photons", data=source_distrs[1])
    sim_file.create_dataset("Gamma", data=gamma_distrs[0] + gamma_distrs[1])
    sim_file.create_dataset("Gamma Electrons", data=gamma_distrs[0])
    sim_file.create_dataset("Gamma Photons", data=gamma_distrs[1])

    for category in categories:
        subcategories: dict[str, list[float]] = categories[category]
        event_type_group: h5py.Group = sim_file.create_group(category)
        for subcategory in subcategories:
            energy_group: h5py.Group = event_type_group.require_group(str(subcategory[0]))
            energy_group.create_dataset(str(subcategory[1]), data=subcategories[subcategory])

    sim_file.close()
    return


@click.command()
@click.argument("config_path", type=click.Path(readable=True, resolve_path=True, path_type=Path))
@click.option("--save-path", '-s', type=click.Path(writable=True, resolve_path=True, path_type=Path))
def main(config_path: Path, save_path: Path) -> int:
    with open(config_path, 'rb') as f:
        config_dict: dict[str, Any] = load(f)

    if "simulation" not in config_dict:
        raise KeyError(f"Configuration is missing the `simulation` config section.")

    simulation_config: dict[str, Any] = config_dict["simulation"]
    processor_count: int = simulation_config.get("num_workers", 1)
    num_events: int = simulation_config.get("num_events", 2)
    source_list: list[str] = simulation_config.get("sources", ["source"])
    detector_name: str = simulation_config.get("detector_name", "detector")
    config_save_path: Path | None = simulation_config.get("save_path", None)

    # Check that there is a save path to use.
    if config_save_path is None and save_path is None:
        raise ValueError("No path to save simulation. Add `save_path` to config's `simulation` section or use `--save-path`.")
    # Prefer the optional save path. If not set, then use the config save path.
    if save_path is None:
        save_path = config_save_path

    if save_path.suffix != ".hdf5":
        print("Missing '.hdf5' in file name. Appending.")
        save_path = save_path.with_name(f"{save_path.name}.hdf5")
    print(f"Saving to: {save_path}")

    sources: list[Source] = []
    for source_name in source_list:
        if source_name not in config_dict:
            raise KeyError(f"Configuration is missing the `{source_name}` config section.")
        sources.append(Source.model_validate(config_dict[source_name]))

    if detector_name not in config_dict:
        raise KeyError(f"Configuration is missing the `{detector_name}` config section.")
    detector: Detector = Detector.model_validate(config_dict[detector_name])
    run_simulation(processor_count, num_events, detector, sources, save_path)
    return 0


if __name__ == "__main__":
    main()
