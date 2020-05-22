# Read in then close an nc file

library(ncdf4)

# Open a netCDF file
nc_file = './sample_nc/a1_20050503_20050504_0221m.adcp.L1.nc'
ncin <- nc_open(nc_file, write=FALSE, readunlim = FALSE, verbose=FALSE)

View(ncin)

ncin

nc_close(ncin)
