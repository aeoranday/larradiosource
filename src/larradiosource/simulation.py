from larradiosource.detector import Detector, Event
from larradiosource.radiation import Source, DecayBranch, Emission

import h5py
import numpy as np
import numpy.typing as npt
from tqdm import tqdm

from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path


def simulate_events(
    num_events: int, worker_id: int, detector: Detector, sources: list[Source], use_tqdm: bool
) -> list[Event]:
    """
    Simulate for n events with the given detector and sources.

    Parameters:
        num_events (int):
            The number of events to simulate.
        detector (Detector):
            The detector to read out from.
        sources (list[Source]):
            The sources that will be decaying.
        use_tqdm (bool):
            Flag on whether to use tqdm progress tracking or not.

    Returns a list of Events that has a length equal to `num_events`.
    """
    events: list[Event] = []

    event_count: int = 0
    pbar = tqdm(total=num_events, position=worker_id, disable=not use_tqdm)

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


def run_simulation(
    processor_count: int,
    num_events: int,
    detector: Detector,
    sources: list[Source],
    use_tqdm: bool,
    sim_save_path: Path,
) -> None:
    """
    Run the simulation and save the results to an HDF5 file.

    Parameters:
        processor_count (int):
            The number of processors to use.
        num_events (int):
            The number of events to produce.
        detector (Detector):
            The detector that will process emissions.
        sources (list[Source]):
            The list of Sources that will produce emissions.
        use_tqdm (bool):
            Flag on whether to use tqdm progress tracking or not.
    """
    num_events_per_processor: list[int] = [
        num_events // processor_count,
    ] * processor_count
    if (remainder_event_count := num_events % processor_count) != 0:
        num_events_per_processor[-1] = num_events_per_processor[-1] + remainder_event_count

    events: list[Event] = []

    with ProcessPoolExecutor(max_workers=processor_count) as executor:
        futures = [
            executor.submit(simulate_events, n, idx, detector, sources, use_tqdm)
            for idx, n in enumerate(num_events_per_processor)
        ]
        for future in futures:
            events += future.result()

    categories: dict[str, dict[tuple[float, int], list[float]]] = categorize_events(events)

    sim_file: h5py.File = h5py.File(sim_save_path, "w")
    sim_file.attrs["Simulation"] = True

    for category in categories:
        subcategories: dict[str, list[float]] = categories[category]
        event_type_group: h5py.Group = sim_file.create_group(category)
        for subcategory in subcategories:
            energy_group: h5py.Group = event_type_group.require_group(str(subcategory[0]))
            energy_group.create_dataset(str(subcategory[1]), data=subcategories[subcategory])

    sim_file.close()
    return


def categorize_events(events: list[Event]) -> dict[str, dict[tuple[float, int], list[float]]]:
    """
    For the list of Events, categorize the event type, energy, and resultant charge.

    Parameters:
        events (list[Event]):
            A list of already processed Events to categorize.

    Returns:
        dict[str, dict[tuple[float, int], list[float]]]:
            A dictionary with 'photon' and 'electron' as base keys and the values as dictionaries
            with 2-tuples for (energy, readout channel) and a list of energies (floats) as the values.
    """
    categories: dict[str, dict[str, list[float]]] = {"photon": defaultdict(list), "electron": defaultdict(list)}
    for event in events:
        emission: Emission = event.init_emission
        subcategory: dict[tuple[float, int], list[float]] = categories[emission.type]
        # Prefer to index as keV and integers.
        subcategory[(emission.energy_mev, event.channels[event.main_channel_index])].append(
            event.amplitudes[event.main_channel_index]
        )

    return categories
