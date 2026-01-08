# Bismuth 207 Simulation
This is a simulation package for Bi-207 radioactive decay in liquid argon time projection chambers (LArTPC).
The target of this simulation is to decay the charge response from the various decays, with a primary focus on the 976 keV and 1048 keV electron emissions through internal conversion.
The gamma decays are included as they will also produce a charge response.

## Detector Geometries Configuration
The focus of this simulation limits detector geometries to prism-like, such as cylindrical and rectangular TPCs.
It is expected that the two example geometries will be sufficient for most use cases for 1D channel coverage through radial extension or Cartesian extension.
The abstract base class of detectors will bring forward whether the emission is contained in the detector and which channels it may be detected on.
The inheriting classes define the calculations required to evaluate containment and activated channels.

## Decay Configuration
This package only focuses on Bi-207 decays in LAr.
However, these decays are parsed from a configuration file that contain the various decay branches.
In theory, one could exchange these decay branches with another radioactive isotope, but the test cases aim to compare against data that was collected from Bi-207.
