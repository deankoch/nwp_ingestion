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
import datetime
import numpy
import pandas as pd
import geopandas
import xarray

# Herbie package for downloading GRIBs and my own extensions
from herbie.archive import Herbie
from helpers import get_nbm_times, get_nbm_variable_names, get_nbm

# %%
'''----------------- browse the NBM collection ----------------'''

# The get_nbm_times function looks in the AWS bucket for directories and files
# corresponding to GRIB2 forecast data. Call the function without arguments to
# get a list of dates for which there is a data directory
dates = get_nbm_times()

# call with a date to get a list of forecast release hours (fhr)
earliest_date = pd.to_datetime('2020-09-30')
print(get_nbm_times(earliest_date))

# the above calls only look for directories. The next one actually identifies
# files that can be downloaded:
# call with a date and fhr to get a list of available forecast hours (fxx)
fxx_avail = get_nbm_times(earliest_date, 0)
print(fxx_avail)

# Not all of the dates returned by get_nbm_times() are accessible. This is
# because we need an index file (*grib2.idx) to read the GRIB2 properly.
# As of December 2021, these are available from 2020-09-30 to present.

# You can check yourself by calling get_nbm_times in a loop. The code below
# looks for index files for the 0-hour release times in a date range around
# the cutoff:
start_date = earliest_date - datetime.timedelta(5)
end_date = earliest_date + datetime.timedelta(5)
for date in pd.date_range(start_date, end_date, freq='D'):
    print(str(date) + ': ' + str(len(get_nbm_times(date, 0))) + ' files')

# The forecast hours are also slightly irregular. This code prints the fxx
# value along with the step length (time since previous forecast)
fxx_step = numpy.concatenate([numpy.array([1]), numpy.diff(fxx_avail)])
print('\n'.join('{} {}'.format(x, y) for x, y in zip(fxx_avail, fxx_step)))
# notice that up to fxx=36 we have hourly forecasts. Beyond that point we
# have 3 and 6 hour intervals (up to hours 192 and 264, respectively)

# %%
# Once you've found a file (date/fhr/fxx) of interest you can check the
# available layers for that timestamp using get_nbm_variable_names
time_dict = {'date': earliest_date}
vnames_all = get_nbm_variable_names(time_dict)
print(vnames_all['name'])
# note that get_nbm_variable_names sets defaults fhr=0, fxx=1, product='co'

# get_nbm_variable_names builds variable names based on 3 fields in the
# GRIB index file, including "variable" (which is not itself unique!) and
# level. This creates a unique naming scheme for the NBM GRIBs. I don't know
# enough about GRIBs yet to know if this will work for other datasets

# %%
'''----------------- Identify layers/times of interest ----------------'''

# We are mostly interested in the subset of variables that can be passed as
# input to SWAT+. We can focus on those by specifying the short "variable"
# names listed below (based on: https://vlab.noaa.gov/web/mdl/nbm-wx-elements)
short_names = {'DSWRF', 'RH', 'TMP', 'WIND', 'APCP'}
vnames_short = vnames_all[vnames_all['variable'].isin(short_names)]
print(vnames_short['name'])

# This shrinks the table to an understandable length - we can see that
# there are often several variables with the same short name, but
# differing in some other attribute (often a vertical "level").
# for example TMP (temperature) has two vertical levels (surface and 2m)
# and a third layer with uncertainty data
vnames_short[vnames_short['variable']=='TMP']
print(vnames_short['name'])
# we will use the 2m TMP as temperature inputs for SWAT ("surface"
# would be better but after downloading it I found it was all NaNs??).

# Notice the variable names depend on the forecast hour (fxx). For example
# the following are valid variable names for fxx=1 hour forecasts:
vnames_toget = {
    'TMP_2_m_above_ground_1_hour_fcst',
    'APCP_surface_01_hour_acc_fcst',
    'RH_2_m_above_ground_1_hour_fcst',
    'DSWRF_surface_1_hour_fcst',
    'WIND_10_m_above_ground_1_hour_fcst'}

# %%
'''----------------- Example: build an fxx=1 archive  ----------------'''
# Given the variable names defined above, let's download and look at some data
# using the get_nbm function.
#
# In SWAT+ we need a fairly long time series to prime the model (building up
# soil water and snowpack, etc) before feeding it future forecast data. For
# this purpose we will build an archival time series of 1-hour-ahead forecasts
# (fxx=1) using the variable names defined above (vnames_toget).
#
# We will ultimately want to store all of the hourly forecasts,
# fhr = [0, 1, ... 23], from the earliest available date to present. But for
# illustration purposes we show how this is done with the earliest 2 dates:
example_enddate = earliest_date + datetime.timedelta(2)
dates = pd.date_range(earliest_date, example_enddate, freq='D')
fhr = list(range(0, 23+1))

# get the data and open as xarray
xdata = get_nbm(dates, vnames_toget, fhr=fhr)
# this function call does a lot of things in sequence - it uses Herbie to
# download the GRIBs (if they aren't in the local storage path), loads them
# and assigns missing attributes/projection info, then combines them into a
# single xarray dataset with a time axis

# print the resulting data object
print(xdata)

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

# plot an example (temperature on the first time in the series)
x['TMP_2_m_above_ground_1_hour_fcst'][0,:,:].plot()

# %%
'''----------------- Area of interest ----------------'''
# It's not possible to get a geographic subset from the GRIB2 files, so we
# are stuck downloading the entire continental USA. However we can clip
# that data very easily once it's opened as xarray:

# define the polygon file to use as a boundary for area of interest
aoi_polygon_dir = 'D:/UYRW_data/side_projects/NWP_ingestion/testdata/input_example/'
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
# plot an example (first 24 hours in series for precip)
example_tseries = xdata_clipped['APCP_surface_01_hour_acc_fcst'][fhr, :, :]
example_tseries.plot(col='time', col_wrap=4)

# %%
# plot an example (first 24 hours in series for precip)
vname = list(vnames_toget)[4]
example_tseries = xdata_clipped[vname][fhr,:,:]
example_tseries.plot(col='time', col_wrap=4)


# %%
'''----------------- Example: APCP and fxx>1 ----------------'''

# We can specify the desired levels (2m and 10m) easily enough for all of the
# variables except APCP (accumulated precipitation).

# APCP has a temporal level that depends on the forecast hour. Up to fxx=36,
# there is only one level, the accumulation over the preceeding hour...
time_dict.update({'fxx': 35})
vnames_check = get_nbm_variable_names(time_dict)
vnames_check[vnames_check['variable']=='APCP']
# ...This is good because the forecast interval is hourly up to fxx=36...

# ...But in later forecast hours the one-hour accumulation is less useful,
# because it only describes what happens in one of the three preceeding hours
# between forecasts. What we really need is the 3-hour accumulated precip, but
# this is not found in all subsequent files!
#
# Instead what we find is that every six hours, the forecast includes a 6-hour
# precip accumulation variable. eg looking at fxx=42, we have the uncertainty
# layer, the 1-hour accumulation, and the 6-hour accumulation:
time_dict.update({'fxx': 42})
vnames_check = get_nbm_variable_names(time_dict)
vnames_toprint = vnames_check[vnames_check['variable']=='APCP']
print(vnames_toprint[['variable', 'forecast_time']])

# So if we fetch the hourly forecasts up to fxx=36, then the six-hourly
# forecasts beyond that (up to fxx=264), we will end up with a time series
# without gaps that can be used to stitch together a daily accumulation.
# Unfortunately this means discarding half of the 3-hourly forecasts.