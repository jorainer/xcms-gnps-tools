# Utility functions for use of xcms with GNPS

Authors: Johannes Rainer, Mar Garcia-Aloy and Michael Witting.

This repository contains some utility functions to integrate `xcms`-based data
processing into the GNPS [Feature-Based Molecular
Networking](https://ccms-ucsd.github.io/GNPSDocumentation/featurebasedmolecularnetworking/)
workflow (FBMN). See [this
page](https://ccms-ucsd.github.io/GNPSDocumentation/featurebasedmolecularnetworking-with-xcms3)
for details on xcms integration into FBMN, and the following [GitHub
repo](https://github.com/DorresteinLaboratory/XCMS3_FeatureBasedMN) for example
scripts as Jupyter notebook and RCommander script.
 
The main provided functions are:

- `formatSpectraForGNPS`: format spectra for MGF export in the format expected
  by GNPS.
- `getEdgelist`: extract a list of *edges* between co-eluting ions potentially
  representing adducts or isotopes of the same metabolite as defined by
  `CAMERA`.
- `getFeatureAnnotations`: extract adduct annotations for features from a
  `CAMERA` result in the format required for Ion Identity Networking (IIN) in
  FBMN/GNPS.
  

