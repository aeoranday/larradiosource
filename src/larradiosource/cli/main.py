from larradiosource.detector import (Detector, Event)
from larradiosource.radiation import (Source, DecayBranch, Emission)
from larradiosource.simulation import run_simulation

import click
import numpy as np
import numpy.typing as npt

from pathlib import Path
from tomllib import load
from typing import Any


@click.command()
@click.argument("config_path", type=click.Path(readable=True, resolve_path=True, path_type=Path))
@click.option("--save-path", '-s', type=click.Path(writable=True, resolve_path=True, path_type=Path))
@click.option("--use-tqdm", '-t', is_flag=True, default=False)
def cli(config_path: Path, save_path: Path, use_tqdm: bool) -> int:
    """ Read configuration, run the LArRadioSource simulation, and save to HDF5.  """
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
    random_seed: int | None = simulation_config.get("random_seed")
    config_use_tqdm: bool | None = simulation_config.get("use_tqdm", False)

    rng: np.random.Generator = np.random.default_rng(seed=random_seed)

    # If the flag is false and config is true, use the config value.
    if not use_tqdm and config_use_tqdm:
        use_tqdm = config_use_tqdm

    # Check that there is a save path to use.
    if config_save_path is None and save_path is None:
        raise ValueError(
                "No path to save simulation. Add `save_path` to config's `simulation` section or use `--save-path`."
        )
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
        source: Source = Source.model_validate(config_dict[source_name])
        source.set_random_generator(rng.spawn(1)[0])
        source.geometry.set_random_generator(rng.spawn(1)[0])
        sources.append(source)

    if detector_name not in config_dict:
        raise KeyError(f"Configuration is missing the `{detector_name}` config section.")
    detector: Detector = Detector.model_validate(config_dict[detector_name])
    detector.set_random_generator(rng.spawn(1)[0])
    detector.geometry.set_random_generator(rng.spawn(1)[0])
    run_simulation(processor_count, num_events, detector, sources, use_tqdm, save_path)
    return 0


if __name__ == "__main__":
    cli()
