import xarray as xr

f = './sample_nc/a1_20050503_20050504_0221m.adcp.L1.nc'

ncdata = xr.open_dataset(f)

# Can access variables and metadata in the netCDF using ncdata.YOUR_VARIABLE.ATTRIBUTE
# Use tab key to view options for YOUR_VARIABLE and ATTRIBUTE

# Close ncdata
xr.Dataset.close(ncdata)
