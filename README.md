# ADCP
`ADCP` is an R package that allows the user to open raw ADCP files, perform QA/QC procedures on them, create preliminary and QA/QC plots from them, and export ADCP data in netCDF format. `ADCP` was written by Emily Chisholm for processing of moored ADCP data at the Bedford Institute of Ocean Sciences, NS. This version is customized for the Institute of Ocean Sciences, BC. 

Currently, there is only one low level of processing supported by this package (level 1). This level comprises:
* Corrections for magnetic declination
* Calculation of sea surface height from pressure values and latitude
* Rotation into enu coordinates if this is not already the coordinate system of the dataset
* Flagging leading and trailing ensembles from before and after deployment and setting them to nan's
* Flagging negative pressure values

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
1. Open R and install `ADCP` dependencies: \
    install.packages(c("devtools", "testthat", "gsw", "ncdf4", "maps"))  
    library(devtools)  
    install_github("dankelley/oce", ref="develop")  
    install_github("dankelley/ocedata", ref="master")  
  
  Note: The `oce` "develop" version is installed because it uses the 13th generation of the International Geomagnetic Reference Field (IGRF-13) for magnetic field calculations while the "master" version still uses IGRF-12. See https://github.com/dankelley/oce/issues/1473 for more.
  
2. Install `ADCP` from R \
    install_github("hhourston/mooredDataProcessing_adcp", ref="master")  
  
  Note: To update `ADCP` (and other R packages) in R with `devtools`: \
    update_packages()  

Alternative method for installing `ADCP` using command line:  
    git clone https://github.com/hhourston/ADCP_processing_visualization.git  

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


