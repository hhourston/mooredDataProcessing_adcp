import xarray as xr

f = 'your_file_name.nc'

d = xr.open_dataset(f)

# Can access variables and metadata in the netCDF using d.YOUR_VARIABLE.COMPONENT
# Use tab key to view options for YOUR_VARIABLE and COMPONENT

# Close d
xr.Dataset.close(d)
