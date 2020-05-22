# ADCP
`ADCP` is an R package that allows the user to open raw ADCP files, perform QA/QC procedures on them, create preliminary and QA/QC plots from them, and export ADCP data in netCDF format. `ADCP` acts as an intermediary between the user and the R packages `oce`, written by Dan Kelley and Clark Richards, and `ncdf4`, written by David Pierce. `ADCP` was written by Emily Chisholm for processing of moored ADCP data at the Bedford Institute of Ocean Sciences, NS. This version is customized for the Institute of Ocean Sciences, BC. 

All modifications to the package were made in the adcpToolbox_P01.R file; the rest of the files are unchanged.

Currently, there is only one low level of processing supported by this package (level 1). This level comprises setting leading and trailing ensembles of a raw ADCP dataset that were taken during deployment and recovery to NA's and applying magnetic declination to the dataset.

## Raw ADCP file types supported
File types supported:
- ".000", ".664" (corresponding to Teledyne RDI Narrowband, Broadband, and Workhorse instruments)
- ".pd0", ".PD0" (corresponding to Teledyne RDI Sentinel V instruments)

File types not supported:
- ".001" (same file format as ".000" and ".664")
- ".arg" (corresponding to Sontek Argonaut XR Multi-Cell Doppler Current Profiler instruments): See `oce` support issue at https://github.com/dankelley/oce/issues/1637

File types for which support is unknown:
- ".666" (same file format as ".000" and ".664")

## Installation
1. Open R and install `ADCP` dependencies:
    install.packages(c("devtools", "testthat", "gsw", "ncdf4", "maps"))
    library(devtools)
    install_github("dankelley/oce", ref="develop")
    install_github("dankelley/ocedata", ref="master")
  
  Note: The `oce` "develop" version is installed because it uses the 13th generation of the International Geomagnetic Reference Field (IGRF-13) for magnetic field calculations while the "master" version still uses IGRF-12. See https://github.com/dankelley/oce/issues/1473 for more.
  
2. Install `ADCP` from R
    install_github("hhourston/mooredDataProcessing_adcp", ref="master")
  
  Note: To update `ADCP` (and other R packages) in R with `devtools`: 
    update_packages()

## Usage
  Please see the the *sample_scripts* folder for processing scripts, and the *sample_data* folder for a sample raw ADCP file and filled-out metadata file.
  
## Credits
  `ADCP` was created by Emily Chisholm (https://github.com/Echisholm21).

## Helpful links
Documentation: \
  `oce` documentation: https://cran.r-project.org/web/packages/oce/oce.pdf \
  `ncdf4` documentation: https://cran.r-project.org/web/packages/ncdf4/ncdf4.pdf \
  netCDF documentation: https://www.unidata.ucar.edu/software/netcdf/docs/index.html \

Conventions: \
  BODC SeaDataNet quality flags: https://www.bodc.ac.uk/data/documents/series/37006/#QCflags \
  BODC SeaDataNet P01 vocabulary search: http://seadatanet.maris2.nl/v_bodc_vocab_v2/search.asp?lib=p01&screen=0 \
  GF3 codes (no longer maintained): https://www.nodc.noaa.gov/woce/woce_v3/wocedata_1/sss/documents/liste_param.htm \
  CF Conventions: http://cfconventions.org/standard-names.html \


