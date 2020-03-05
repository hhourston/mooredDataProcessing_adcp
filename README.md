# ADCP
`ADCP` is an R package that allows the user to open raw ADCP files, perform QA/QC procedures on them, create preliminary and QA/QC plots from them, and export ADCP data in netCDF format. `ADCP` acts as an intermediary between the user and the R packages `oce`, written by Dan Kelley and Clark Richards, and `ncdf4`, written by David Pierce. `ADCP` was written by Emily Chisholm for processing of moored ADCP data at the Bedford Institute of Ocean Sciences, NS. This version is customized for the Institute of Ocean Sciences, BC. All modifications to the package were made in the adcpToolbox_P01.R file; the rest of the files are unchanged.

Currently, there is only one low level of processing supported by this package (level 1). This level comprises setting leading and trailing ensembles of a raw ADCP dataset that were taken during deployment and recovery to NA's and applying magnetic declination to the dataset.

## Installation
  To install dependencies in R:
  
    install.packages("devtools")
    install.packages("testthat")
    install.packages("gsw")
    install_github("dankelley/oce", ref="develop")
    install_github("dankelley/ocedata", ref="master")
    install.packages("ncdf4")
    install.packages("maps")
  
  The `oce` "develop" version is installed because it uses the 13th generation of the International Geomagnetic Reference Field (IGRF-13) for magnetic field calculations while the "master" version still uses IGRF-12. See https://github.com/dankelley/oce/issues/1473 for more.
  
  To install ADCP in R:
  
    library(devtools)
    install_github("hhourston/mooredDataProcessing_adcp", ref="master")
  
## Documentation of dependencies
  `oce`: https://cran.r-project.org/web/packages/oce/oce.pdf \
  `ncdf4`: https://cran.r-project.org/web/packages/ncdf4/ncdf4.pdf

## Usage
  Please see the repository https://github.com/hhourston/ADCP_processing_visualization for usage instructions. A sample script called *ADCP_lvl1_process.R* that makes use of this package can be found in this repository along with a csv metadata template, a sample raw file, filled-out csv metadata file and output netCDF file, and R and Python scripts for viewing and plotting netCDF file data.
  
  A version of this script that is callable from Python can be found at https://github.com/hhourston/ADCP_processing_visualization/tree/master/callR_fromPython. A sample Python script that runs this R script can also be found in that folder. Usage instructions can be found in https://github.com/hhourston/ADCP_processing_visualization.
  
## Credits
  `ADCP` was created by Emily Chisholm (https://github.com/Echisholm21).
  
