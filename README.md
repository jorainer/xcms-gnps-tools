# Utility functions for use of xcms with GNPS

Authors: Johannes Rainer, Mar Garcia-Aloy and Michael Witting.

This repository contains some utility functions to integrate `xcms`-based data
processing into the GNPS [Feature-Based Molecular
Networking](https://ccms-ucsd.github.io/GNPSDocumentation/featurebasedmolecularnetworking/)
workflow (FBMN). See [this
page](https://ccms-ucsd.github.io/GNPSDocumentation/featurebasedmolecularnetworkingwith-xcms3)
for details on xcms integration into FBMN.
 
Two main functions are defined:

- `formatSpectraForGNPS`: format spectra for MGF export in the format expected
  by GNPS.
- `getEdgelist`: extract a list of *edges* between co-eluting ions potentially
  representing adducts or isotopes of the same metabolite as defined by
  `CAMERA`.

