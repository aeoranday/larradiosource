from larradiosource.detector import (Detector, Event)
from larradiosource.radiation import (Source, DecayBranch, Emission)
from larradiosource.simulation import run_simulation

import click
import numpy as np
import numpy.typing as npt

from pathlib import Path
from tomllib import load
from typing import Any


@click.group()
def cli():
    """LArRadioSource -- Simulate radioactive sources deployed in LArTPCs."""
    pass


@cli.command("simulate")
@click.argument("config_path", type=click.Path(readable=True, resolve_path=True, path_type=Path))
@click.option("--output", '-o', type=click.Path(writable=True, resolve_path=True, path_type=Path))
@click.option("--use-tqdm", '-t', is_flag=True, default=False)
def simulate(config_path: Path, output: Path, use_tqdm: bool) -> int:
    """ Read configuration, simulate, and save to HDF5.  """
    with open(config_path, 'rb') as f:
        config_dict: dict[str, Any] = load(f)

    if "simulation" not in config_dict:
        raise KeyError(f"Configuration is missing the `simulation` config section.")

    simulation_config: dict[str, Any] = config_dict["simulation"]
    processor_count: int = simulation_config.get("num_workers", 1)
    num_events: int = simulation_config.get("num_events", 2)
    source_list: list[str] = simulation_config.get("sources", ["source"])
    detector_name: str = simulation_config.get("detector_name", "detector")
    config_output: Path | None = simulation_config.get("output", None)
    random_seed: int | None = simulation_config.get("random_seed")
    config_use_tqdm: bool | None = simulation_config.get("use_tqdm", False)

    rng: np.random.Generator = np.random.default_rng(seed=random_seed)

    # If the flag is false and config is true, use the config value.
    if not use_tqdm and config_use_tqdm:
        use_tqdm = config_use_tqdm

    # Check that there is a save path to use.
    if config_output is None and output is None:
        raise ValueError(
                "No output path to save simulation. Add `output` to config's `simulation` section or use `--output`/`-o`."
        )
    # Prefer the optional save path. If not set, then use the config save path.
    if output is None:
        output = config_output

    if output.suffix != ".hdf5":
        print("Missing '.hdf5' in file name. Appending.")
        output = output.with_name(f"{output.name}.hdf5")
    print(f"Saving to: {output}")

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
