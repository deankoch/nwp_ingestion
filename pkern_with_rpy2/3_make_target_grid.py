# Define the output grid to interpolate onto 
# Dean Koch
# January 11, 2022
#
# OVERVIEW
#
# This script creates an xarray object specifying the desired grid
# dimensions of the kriging output from pkern. It writes the file
# target_grid_simple.nc to testdata/pkern_examples
#
# For this example I simply take the original temperature data grid
# (pkern_examples/temperature_series.nc) and increase its resolution
# by 4x.
#
# However the target grid can be configured however you like. In
# another example script I'll show how to use an elevation raster
# as a target grid.

import numpy as np
import xarray as xr


# %%
'''------- define working directory and filepaths --------'''

# path to data subdirectory for output
nwp_dir = '/mnt/e/coding-projects/nwp_ingestion'
pkern_data_dir = nwp_dir + '/testdata/pkern_examples'

# define path to temperature file and desired output file
temperature_path = pkern_data_dir + '/temperature_series.nc'
grid_out_path = pkern_data_dir + '/target_grid_simple.nc'

# define a downscaling factor
ds_fact = 4


# %%
'''------- open the temperature grid --------'''

# open the netcdf file as xarray.Dataset
temperature_data = xr.open_dataset(temperature_path)

# extract projection information
crs = temperature_data['spatial_ref'].crs_wkt

# extract grid dimensional info for x
nx = temperature_data.dims['x']
minx = temperature_data['x'].min()
maxx = temperature_data['x'].max()

# and for y
ny = temperature_data.dims['y']
miny = temperature_data['y'].min()
maxy = temperature_data['y'].max()


# %%
'''------- define a higher-res grid --------'''

# define new grid lines with increased resolution on each axis
glx = np.linspace(minx, maxx, ds_fact * nx)
gly = np.linspace(miny, maxy, ds_fact * ny)

# construct the DataArray and append projection info
target_grid = xr.DataArray(coords={'x': glx, 'y': gly}, dims=['x', 'y'])
target_grid = target_grid.rio.write_crs(crs)

# save to netcdf
target_grid.to_netcdf(path=grid_out_path)