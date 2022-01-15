# Get example data for pkern demo
# Dean Koch
# January 11, 2022
#
# OVERVIEW
#
# This script uses Herbie to download some temperature grids to
# use as examples in the pkern kriging demo scripts.
#
# It creates the directory nwp_ingestion/testdata/pkern_examples
# and writes the file temperature_series.nc, containing 4 time
# slices of the NBM 2m temperature variable.

import os
import geopandas

# import my extensions of Herbie package for downloading GRIBs
from herbie_helpers import get_nbm


# %%
'''------- define working directory and filepaths --------'''

# path to nwp_ingestion project and new data subdirectory
nwp_dir = '/mnt/e/coding-projects/nwp_ingestion'
data_dir = nwp_dir + '/testdata'
pkern_data_dir = data_dir + '/pkern_examples'
if not os.path.isdir(pkern_data_dir):
    print('creating data directory at: ' + pkern_data_dir)
    os.mkdir(pkern_data_dir)

# define the polygon file to use as a boundary for AOI
aoi_poly_path = data_dir + '/input_example/UYRW_boundary_25kmpadding.geojson'

# define path to output file
temperature_out_path = pkern_data_dir + '/temperature_series.nc'


# %%
'''------- specify an example item from NBM collection --------'''

# to list options, use get_nbm_times and get_nbm_variable_names
# from herbie_helpers eg: get_nbm_times(date=eg_date, fhr=0)

# set an example date-time series to download
eg_date = '2022-01-11'                          # release date
eg_fhr = [0, 6, 12, 18]                         # release hour(s)
eg_fxx = 1                                      # forecast hour
eg_vname = 'TMP_2_m_above_ground_1_hour_fcst'   # variable name


# %%
'''------- download the data and crop to study area --------'''

# get the data and open as xarray
xdata = get_nbm(eg_date, {eg_vname}, fxx=eg_fxx, fhr=eg_fhr)
crs = xdata.rio.crs

# open  polygon, transform (as needed) to match xdata projection
aoi_polygon = geopandas.read_file(aoi_poly_path).to_crs(crs)
aoi_bbox = aoi_polygon.bounds
minx = aoi_bbox['minx'][0]
miny = aoi_bbox['miny'][0]
maxx = aoi_bbox['maxx'][0]
maxy = aoi_bbox['maxy'][0]

# clip to AOI and save a copy as netcdf
temperature_data = xdata.rio.clip_box(minx, miny, maxx, maxy)
temperature_data.to_netcdf(path=temperature_out_path)

