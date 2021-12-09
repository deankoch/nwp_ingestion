# Test download from NOAA's National Blended Model (NBM)
# Based on: readme tutorial at Herbie (AKA "hrrrb") package github repository
# https://blaylockbk.github.io/Herbie/_build/html/ (visited 30/11/2021)
#
# DESCRIPTION
#
# The script downloads forecast data for the specified times/dates/variables
# then imports it into python as an xarray dataset, where it can be more easily
# converted into netCDF, ZARR, geoJSON, or ASCII tables (as needed by SWAT).
# 
# INSTALLATION
#
# See herbie_install.txt for instructions on getting started in Windows 10
#
# WHAT'S NEXT 
#
# We may want to load a polygon (path supplied by the user) and crop the
# raster layers to its extent. Wote that we can't download geographical subsets
# of the GRIBs, so ultimately we are stuck downloading for the full continental
# USA, and we may wish to keep these files archived in case we want to look at
# other study areas in the future.
#
# For the next step in the SWAT+ workflow I will want to save the cropped data
# as GeoJSON (containing the grid point locations), and a JSON containing the
# data in name-value pairs (by date?).
#
# TODO: Input/output variables used in this script
#
# name                  description
# ---------------------------------------------------------------------------
# aoi_polygon_path      local path to input polygon (GeoJSON) file
# ---------------------------------------------------------------------------
#

# %%
#
# There are two essential pieces of information we need from the GRIB files:
# the appropriately referenced data layers (ie with CRS info), and the layer
# attributes (ie what is band 1? band 2?).
#
# I have not been able to find a single python package that can do both, so
# we use "rioxarray" to load the data layers as xarray objects, then assign
# attributes loaded from the GRIB using "pygrib". These modules are loaded
# in "helpers.py"

# handle dataframes, geometries
import pandas
import geopandas
import xarray
import rioxarray

# Herbie package for downloading GRIBs
from herbie.archive import Herbie

# my own helper functions
from helpers import H2x, get_variable_names, get_nbm

# %%
# browse the NBM collection
#
# I'm using my own function instead of Herbie.xarray to open the GRIB files
# because there seems to be a bug in cfgrib (the module that Herbie uses)
# that causes problems with NBM. I think this has to do with non-uniqueness of
# variable names in the NBM GRIBs.

# fetch a table of variables for the NBM
vnames_all = get_nbm('2021-01-01')
# note that some dates/times will have only a subset of these

# define a set of variables of interest
# based on: https://vlab.noaa.gov/web/mdl/nbm-wx-elements
short_names = {'DSWRF', 'RH', 'TMP', 'APCP', 'WIND'}

# Note that in general, "short_name" does not uniquely identify a layer in
# a file. There can be many variables with the same short name, but
# differing in some other attribute (often a vertical or temporal "level").

# print the layers available for these short names
vnames_short = vnames_all[vnames_all['variable'].isin(short_names)]
print(vnames_short['name'])

# %%
# define layers, date range and time interval to fetch

# note that the APCP variable is an aggregated precip value for the previous
# hour - this means we really need hourly data for it to be useful going
# forward. eg we would take the sum over the 24 hours of the day to get total
# precip for a particular date.
# variable names from on vnames_all
vnames_toget = {
    'TMP_2_m_above_ground_1_hour_fcst',
    'APCP_surface_01_hour_acc_fcst',
    'RH_2_m_above_ground_1_hour_fcst',
    'DSWRF_surface_1_hour_fcst',
    'WIND_10_m_above_ground_1_hour_fcst'}

# for now I use a few recent dates as examples. Ultimately it would be good to
# run this on the complete available time series (2021-01-01 to present)
# dates to compile
dates = pandas.date_range('2021-11-06', '2021-11-08', freq='D')

# We will want to download all hours with fhr = [0, 1, ..., 23].
# forecast release time of day: integer or list with elements in 0,..23
fhr = list(range(0, 24))

# %%
# get the data and open as xarray
xdata = get_nbm(dates, vnames_toget, fhr=fhr)
# this function call does a lot of things in sequence - it uses Herbie to
# download the GRIBs (if they aren't in the local storage path), loads them
# and assigns missing attributes/projection info, then combines them into a
# single xarray dataset with a time axis

# print the resulting data object
xdata

# %%
# save to netCDF
# note: setting encoding fields to specify compression using gzip causes
# problems with the CRS fields in the netCDF file (files could be loaded
# but CRS field not present in both R and QGIS)
# enc_template = dict(zlib=True, complevel=5)
# enc_dict = {var: enc_template for var in xarray_output.data_vars}
xdata.to_netcdf('test_xarray.nc')

# check that the file can be reloaded
x = xarray.open_dataset('test_xarray.nc', decode_coords=True)
print(x)

# %%
# plot an example (temperature on the first time in the series)
x['TMP_2_m_above_ground_1_hour_fcst'][0,:,:].plot()


# %%
'''----------------- Area of interest ----------------'''

# define the polygon file to use as a boundary for area of interest
aoi_polygon_dir = 'D:/UYRW_data/side_projects/python_NWP/data/input_example/'
aoi_polygon_filename = 'UYRW_boundary_25kmpadding.geojson'
aoi_polygon_path = aoi_polygon_dir + aoi_polygon_filename

# open the polygon and transform (as needed) to match projection of xdata 
aoi_polygon = geopandas.read_file(aoi_polygon_path).to_crs(xdata.rio.crs)
aoi_bbox = aoi_polygon.bounds
minx = aoi_bbox['minx']
miny = aoi_bbox['miny']
maxx = aoi_bbox['maxx']
maxy = aoi_bbox['maxy']

# clip to AOI and save a copy as netcdf
xdata_clipped = xdata.rio.clip_box(minx, miny, maxx, maxy)
xdata.to_netcdf('test_xarray2.nc')

# %%
# plot an example ()
xdata_clipped['APCP_surface_01_hour_acc_fcst'][fhr,:,:].plot(col='time', col_wrap=4)