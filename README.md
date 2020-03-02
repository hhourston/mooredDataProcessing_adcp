# ADCP
`ADCP` is an R package that allows you to open raw ADCP files, perform QA/QC procedures on them, create preliminary and QA/QC plots from them, and export ADCP data in netCDF format. `ADCP` acts as an intermediary between the user and the R packages `oce`, written by Dan Kelley and Clark Richards, and `ncdf4`, written by David Pierce. `ADCP` was written by Emily Chisholm for the Bedford Institute of Ocean Sciences. This version is customized for the Institute of Ocean Sciences.

Currently, there is only one low level of processing supported by this package (level 1). This level comprises setting leading and trailing ensembles of a raw ADCP dataset that were taken during deployment and recovery to NA's and applying magnetic declination to the dataset.

This repository is under construction.

## Installation
  To install dependencies in R:
  
    install.packages("devtools")
    install_github("dankelley/oce", ref="develop")
    install_github("dankelley/ocedata", ref="master")
    install.packages("ncdf4")
    install.packages("maps")
  
  To install ADCP in R:
  
    library(devtools)
    install_github("hhourston/mooredDataProcessing_adcp", ref="master")
  
## Documentation of dependencies
  `oce`: https://cran.r-project.org/web/packages/oce/oce.pdf \
  `ncdf4`: https://cran.r-project.org/web/packages/ncdf4/ncdf4.pdf

## Usage
  A sample script called *ADCP_lvl1_process.R* that makes use of this package can be found at https://github.com/hhourston/ADCP_processing_visualization.
  
  A version of this script that is callable from Python can be found at https://github.com/hhourston/ADCP_processing_visualization/tree/master/callR_fromPython. A sample Python script that runs this R script can also be found in that folder.
  
## Credits
  `ADCP` was created by Emily Chisholm (https://github.com/Echisholm21).
  
