# Utility functions for use of xcms with GNPS

Authors: Johannes Rainer, Mar Garcia-Aloy, Philippine Louail and Michael 
         Witting.

This repository contains some utility functions to integrate `xcms`-based data
processing into the GNPS [Feature-Based Molecular
Networking](https://ccms-ucsd.github.io/GNPSDocumentation/featurebasedmolecularnetworking/)
workflow (FBMN). See [this
page](https://ccms-ucsd.github.io/GNPSDocumentation/featurebasedmolecularnetworking-with-xcms3)
for details on xcms integration into FBMN, and the following [GitHub
repo](https://github.com/DorresteinLaboratory/XCMS3_FeatureBasedMN) for example
scripts Rmarkdown documents as Jupyter notebooks.
 
The main provided functions are:

- `formatSpectraForGNPS()`: format spectra for MGF export in the format expected
  by GNPS.
- `maxTicPeaksData()`: helper function to be used with
  `Spectra::combineSpectra()` to select the peak matrix from the fragment
  spectra with the largest TIC (sum of all fragment intensities) for a group of
  fragment spectra.
- `getEdgelist()`: extract a list of *edges* between co-eluting ions potentially
  representing adducts or isotopes of the same metabolite as defined by
  *CAMERA*.
- `getFeatureAnnotations()`: extract adduct annotations for features from a
  *CAMERA* result in the format required for Ion Identity Networking (IIN) in
  FBMN/GNPS.
  

