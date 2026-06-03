# LArTPC Radioactive Source Simulation `larradiosource`
This is a simulation package for radioactive sources in liquid argon time projection chambers (LArTPC).
The target of this simulation is to produce the charge response from a radioactive source.
The initial use-case of this package was for Bi-207 internal conversion decays, so the development into other decay schemes is left open.

![LArRadioSource UML Diagram](docs/uml/larradiosource.svg "LArRadioSource UML Diagram")

## Detector Geometries Configuration
The focus of this simulation limits detector geometries to prism-like TPCs, such as cylinders and rectangular prisms.
It is expected that the two example geometries will be sufficient for most use cases for 1D channel coverage through radial extension or Cartesian extension.
Details on the `Detector` and `Geometry` configurations are given in [docs/configuration/detector.md](docs/configuration/detector.md) and [docs/configuration/geometry.md](docs/configuration/geometry.md).

## Decay Configuration
Decays branches are written according to their emissions and probability.
Each emission should have a type, e.g. `photon` or `electron`, an energy (MeV), and conversion coefficients in the case of IC.
Since the configuration is fully user-controlled, it is on the user to know these values and any calculations required to getting these values.
References and using Nuclear Data Sheets and similar would be useful for these configurations.
Details on the radiation-related configurations are in [docs/configuration/radiation.md](docs/configuration/radiation.md).

## Simulation Configuration
There are a few fields that are used by the top-level simulator, such as a save path, number of workers to spawn, number of emissions to create, and the names of the detector and source(s) to be used.
For scripting purposes, a few of these configurable parameters are also CLI options and can be seen using `larradiosource simulate --help`.

## Examples
There are various example configurations kept in the `config` directory; the most extensive is the `np04_simulation.toml`, modeled after the [ProtoDUNE 2 Horizontal Drift detector](https://arxiv.org/pdf/1706.07081).
Running a simulation with the `np04_simulation.toml` config and saving to the current working directory would use the command
```bash
larradiosource simulate config/np04_simulation.toml -o test_simulation.hdf5
```

## File Structure
The simulated contents are written to HDF5 with the specified output (either CLI option or written in configuration file).
The structure of the HDF5 file contains several groups:
1. Top level `photon` or `electron` group for the original emission type.
2. Energy level in MeV as a float (there will be awkward precision keys due to this).
3. Channel number (if any activity on that channel).

The contents of group 3 will only be channels that overlap with the energy deposition; channels that experience no energy throughout the simulation will not be written to file.
Similarly, contents of group 2 will only be for energies that occurred randomly, e.g., a low probability emission that happened to never occur in the simulation will not be written to file.
However, group 1 will always have the `photon` and `electron` split.
It may be that the contents of either group is empty, but the group will still exist.

For example, an HDF5 path could appear as `/photon/0.569698/100` and the underlying dataset would contain the charge deposits on channel 100 from a 0.569698 MeV photon emission.
A path with awkward precision may be `/electron/1.0598049999999999/203` for channel 203 from (ideally) a 1.059805 MeV electron emission.
It is more intuitive to use HDF5's keys to traverse the contents of these files.
