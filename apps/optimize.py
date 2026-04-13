from bi207_simulation.detector import Detector, Event
from bi207_simulation.radiation import Source, DecayBranch, Emission
from bi207_simulation.signals import PurityMonitorSignalProcessor, Waveform

import click
import h5py
import numpy as np
import numpy.typing as npt
from tqdm import tqdm

from collections import defaultdict
from collections.abc import Callable
import csv
from pathlib import Path
from pprint import pprint
from tomllib import load
from typing import Any

from concurrent.futures import ProcessPoolExecutor

BIN_STEP: int = 20
ADC_TO_FC: npt.NDArray[float] = np.array([2.209, 0.063]) * 1e-3  # fC/ADC
BINS: ArrayLike = np.arange(0, 4000+BIN_STEP, BIN_STEP) * ADC_TO_FC[0]
FC_STEP: float = BIN_STEP * ADC_TO_FC[0]

DATA_PATH: Path = Path("~/Documents/Code/Bi-207 Simulation/Validation/data/anode_tp_triplets_charge_integral_spectrum_20241122T101809.hdf5").expanduser()
SIM_SAVE_PATH: Path = Path("~/Documents/Code/Bi-207 Simulation/Validation/data/simulated_charge_spectrum.hdf5").expanduser()
CSV_BASE_PATH: Path = Path("~/Documents/Code/Bi-207 Simulation/").expanduser()

MASK: npt.NDArray[int] = np.arange(19, 110)
DATA: h5py.File = h5py.File(DATA_PATH, 'r')
DATA_DISTRIBUTION: npt.NDArray[int] = DATA["Source"][MASK]# - DATA["Mirror"][MASK]


def calculate_charge_distribution(events: list[Event], max_count_channel: int) -> npt.NDArray[int]:
    """
    Calculate the output charge distribution.
    """
    amplitudes = {"photon": defaultdict(list), "electron": defaultdict(list)}

    for event in events:
        amplitude: float = event.amplitudes[event.main_channel_index]
        main_channel: int = event.channels[event.main_channel_index]
        amplitudes[event.type][main_channel].append(amplitude)

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


def calculate_chi2(sim_distribution: npt.NDArray[int]) -> float:
    chi2: float = np.sum( (sim_distribution[MASK] - DATA_DISTRIBUTION)**2 / DATA_DISTRIBUTION )
    return chi2


def calculate_emd(sim_distribution: npt.NDArray[int]) -> float:
    sim_distribution = sim_distribution[MASK] / np.sum(sim_distribution[MASK])
    delta: npt.NDArray[float] = np.abs(np.cumsum(sim_distribution) - np.cumsum(DATA_DISTRIBUTION))
    return np.sum(delta[:-1] * FC_STEP)


def get_primary_channel(events: list[Event]) -> int:
    channels, counts = np.unique([event.channels[event.main_channel_index] for event in events],
                                 return_counts=True)

    return channels[np.argmax(counts)]


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
    pbar = tqdm(total=num_events, position=worker_id, disable=True)
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
            pbar.update(1)
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


def collect_stat_distances(processor_count: int,
                           num_events: int,
                           stat_divergence: str,
                           detector: Detector,
                           sources: list[Source],
                           fake_scale: float = 0,
                           elifetime: float = 0,
                           signal_width: float = 0,
                           channel_resolution_mev: float = 0,
                           save: bool = False,
                           ) -> None:
    np.random.seed(42)
    num_events_per_processor: list[int] = [num_events // processor_count,] * processor_count
    if (remainder_event_count := num_events % processor_count) != 0:
        num_events_per_processor[-1] = num_events_per_processor[-1] + remainder_event_count

    match stat_divergence:
        case "chi2":
            distance_fn: Callable[[npt.NDArray[int]], float] = calculate_chi2
        case "emd":
            distance_fn: Callable[[npt.NDArray[int]], float] = calculate_emd
        case _:
            raise NotImplemented(f"Requested statistical divergence '{stat_divergence}' is not implemented.")

    events: list[Event] = []
    if fake_scale > 0:
        detector.fake_scale = fake_scale
    if elifetime > 0:
        detector.electron_lifetime_ms = elifetime
    if signal_width > 0:
        detector.signal_width = signal_width
    if channel_resolution_mev > 0:
        detector.channel_resolution_mev = channel_resolution_mev

    with ProcessPoolExecutor(max_workers=processor_count) as executor:
        futures = [executor.submit(simulate_events, n, idx, detector, sources) for idx, n in enumerate(num_events_per_processor)]
        for future in futures:
            events += future.result()

    max_count_channel: int = get_primary_channel(events)
    print("Primary Channel Number:", max_count_channel)
    source_distrs, gamma_distrs = calculate_charge_distribution(events, max_count_channel)
    distance: float = distance_fn(source_distrs[0] + source_distrs[1])

    if save:
        categories: dict[str, dict[tuple[int, int], list[float]]] = categorize_events(events)
        print(f"{stat_divergence} = {distance:.2f}")
        sim: h5py.File = h5py.File(SIM_SAVE_PATH, 'w')
        sim.attrs['bins'] = BINS
        sim.attrs[stat_divergence] = distance
        sim.create_dataset("Source", data=source_distrs[0] + source_distrs[1])
        sim.create_dataset("Source Electrons", data=source_distrs[0])
        sim.create_dataset("Source Photons", data=source_distrs[1])
        sim.create_dataset("Gamma", data=gamma_distrs[0] + gamma_distrs[1])
        sim.create_dataset("Gamma Electrons", data=gamma_distrs[0])
        sim.create_dataset("Gamma Photons", data=gamma_distrs[1])

        for category in categories:
            subcategories: dict[str, list[float]] = categories[category]
            event_type_group: h5py.Group = sim.create_group(category)
            for subcategory in subcategories:
                energy_group: h5py.Group = event_type_group.require_group(str(subcategory[0]))
                energy_group.create_dataset(str(subcategory[1]), data=subcategories[subcategory])


        sim.attrs['Simulation'] = True
        sim.close()


    csv_filename: str = ""
    if fake_scale > 0:
        csv_filename = f"locations_{stat_divergence}_w{detector.signal_width}-r{detector.channel_resolution_mev}.csv"
    elif signal_width > 0:
        csv_filename = f"widths_{stat_divergence}_f{detector.fake_scale}-e{detector.electron_lifetime_ms}.csv"

    if len(csv_filename) > 0:
        csv_file_path: Path = CSV_BASE_PATH.joinpath(csv_filename)
        with open(csv_file_path, 'a') as f:
            writer = csv.writer(f)
            writer.writerow([detector.fake_scale,
                             detector.electron_lifetime_ms,
                             detector.signal_width,
                             detector.channel_resolution_mev,
                             distance])
    return


@click.command()
@click.argument("config_path", type=click.Path(readable=True, resolve_path=True, path_type=Path))
@click.option("--save", '-s', is_flag=True, default=False, help="Saves the simulation event histograms to an HDF5 file.")
def main(config_path: Pathm, save: bool) -> int:
    with open(config_path, 'rb') as f:
        config_dict: dict[str, Any] = load(f)

    if "simulation" not in config_dict:
        raise KeyError(f"Configuration is missing the `simulation` config section.")

    simulation_config: dict[str, Any] = config_dict["simulation"]
    processor_count: int = simulation_config.get("num_workers", 1)
    num_events: int = simulation_config.get("num_events", 2)
    source_list: list[str] = simulation_config.get("sources", ["source"])
    detector_name: str = simulation_config.get("detector_name", "detector")
    stat_divergence: str = simulation_config.get("stat_divergence", "chi2")

    sources: list[Source] = []
    for source_name in source_list:
        if source_name not in config_dict:
            raise KeyError(f"Configuration is missing the `{source_name}` config section.")
        sources.append(Source.model_validate(config_dict[source_name]))

    if detector_name not in config_dict:
        raise KeyError(f"Configuration is missing the `{detector_name}` config section.")
    detector: Detector = Detector.model_validate(config_dict[detector_name])

    fake_scales: npt.NDArray[float] = np.arange(0.80, 0.83, 0.0025)
    elifetimes: npt.NDArray[float] = np.arange(3.0, 10, 0.5)
    f, e = np.meshgrid(fake_scales, elifetimes)
#    signal_widths: npt.NDArray[int] = np.arange(2, 11)
#    channel_resolutions: npt.NDArray[float] = np.arange(0.005, 0.1, 0.005)
#    w, r = np.meshgrid(signal_widths, channel_resolutions)

    if save:
        print(f"Saving simulation to {SIM_SAVE_PATH}")
        collect_stat_distances(processor_count, num_events, stat_divergence, detector, sources, save=save)
        return 0

    pbar = tqdm(total=len(elifetimes) * len(fake_scales))
    for row in range(len(elifetimes)):
        for col in range(len(fake_scales)):
            collect_stat_distances(processor_count, num_events, stat_divergence, detector, sources,
                                   fake_scale=f[row, col],
                                   elifetime=e[row, col],
#                                   signal_width=w[row, col],
#                                   channel_resolution_mev=r[row, col],
                                   )
            pbar.update(1)
    return 0


if __name__ == "__main__":
    main()
