# Test download from NOAA's High Resolution Rapid Refresh  
# Based on: readme tutorial at Herbie (AKA "hrrrb") package github repository
# https://blaylockbk.github.io/Herbie/_build/html/ (visited 30/11/2021)
#
# DESCRIPTION
#
# The script loads a polygon (path supplied by the user) and determines
# the corresponding GRIB files to download. These rasters are opened and
# a mosaic is created, clipped, and masked, for a selection of model variables.
# The result is then saved in two files: a GeoJSON containing the grid point
# locations, and a JSON containing the data in name-value pairs.
#
# TODO: Input/output variables used in this script
#
# name                  description
# ---------------------------------------------------------------------------
# aoi_polygon_path      local path to input polygon (GeoJSON) file
# ---------------------------------------------------------------------------
#
# INSTALLATION
#
# See herbie_install.txt for instructions on getting started in Windows 10
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

# print the layers available for these short names
vnames_short = vnames_all[vnames_all['variable'].isin(short_names)]
print(vnames_short['name'])

# %%
# define layers, date range and time interval to fetch

# variable names from on vnames_all
vnames_toget = {
    'TMP_surface_1_hour_fcst',
    'APCP_surface_01_hour_acc_fcst',
    'RH_2_m_above_ground_1_hour_fcst',
    'DSWRF_surface_1_hour_fcst',
    'WIND_10_m_above_ground_1_hour_fcst'}

# dates to compile
dates = pandas.date_range('2021-11-06', '2021-11-08', freq='D')

# forecast release time of day: integer or list with elements in 0,..23
fhr = [0, 12]

# %%
# get the data and open as xarray
xdata = get_nbm(dates, vnames_toget, fhr=fhr)
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

# It's not clear how to name the layers in a GRIB file. Layers typically
# have dozens to hundreds of attributes, many of which are identical. The
# obvious solution would be to concatenate a selection of keyword attributes
# like "short_name" and "level", but I haven't found a combination that
# yields simple and unique names accross the whole file. This makes it
# challenging to write code that will work for any GRIB file
# 
# Why not just refer to a layer by its name? Because GRIB layers don't have
# unique names! They are numbered, but the order and content varies from
# file to file. Instead, in each GRIB layer you will find dozens to hundreds
# of attributes. We can't rely on any of these to be unique across layers.
#
#
# Note that in general, "short_name" does not uniquely identify a layer in
# a file. There can be many variables with the same short name, but
# differing in some other attribute (often a vertical or temporal "level").
#
# In this case we are looking for specific variable short names, paired with
# specific height (m) levels, and in the case of precip, a specific time
# period (0-1 hours). This should yield a unique set of 5 variables found


