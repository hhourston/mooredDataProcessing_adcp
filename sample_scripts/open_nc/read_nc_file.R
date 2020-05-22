# Read in then close an nc file

library(ncdf4)

#The first file to do
ncin1 <- nc_open('D:/ADCP_processing/ADCP/a1_20160713_20170513_0480m/a1_20160713_20170513_0480m_.adcp.L1.nc', write=FALSE, readunlim = FALSE, verbose=FALSE)

View(ncin1)

ncin1

nc_close(ncin1)
