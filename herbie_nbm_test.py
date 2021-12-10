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

# %%
'''----------------- Identify layers/times of interest ----------------'''

# We are mostly interested in the subset of variables that can be passed as
# input to SWAT+. We can focus on those by specifying the short names listed
# below (based on: https://vlab.noaa.gov/web/mdl/nbm-wx-elements)
short_names = {'DSWRF', 'RH', 'TMP', 'WIND', 'APCP'}
vnames_short = vnames_all[vnames_all['variable'].isin(short_names)]
print(vnames_short['name'])

# In general, "short_name" does not uniquely identify a layer in
# a file. There can be many variables with the same short name, but
# differing in some other attribute (often a vertical "level").
# for example TMP (temperature) has two vertical levels (surface and 2m)
# and a third layer with uncertainty data
vnames_short[vnames_short['variable']=='TMP']
print(vnames_short['name'])

# get_nbm_variable_names builds variable names based on 3 fields in the
# GRIB index file, including name (which is not itself unique!) and level.
# This creates a unique naming scheme for the NBM GRIBs. I don't know enough
# about GRIBs yet to know if this will work for other datasets

# We can specify the desired levels (2m and 10m) easily enough for all of the
# variables except APCP (accumulated precipitation).
# 
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

# %%
'''----------------- Example: building an fxx=1 archive ----------------'''

# In SWAT+ we need a fairly long time series to prime the model (building up
# soil water and snowpack, etc) before feeding it future forecast data. For
# this purpose we will build an uninterrupted archival time series of
# 1-hour-ahead forecasts (fxx=1).
# 
# Later on, we can append the latest forecasts (on current date, fxx > 1) to
# extend this time series 10 days into the future.
#
# We saw earlier that variable names, Since the fxx is constant, every forecast request will be for the same
# five variable names:

# variable names from on vnames_all
vnames_toget = {
    'TMP_2_m_above_ground_1_hour_fcst',
    'APCP_surface_01_hour_acc_fcst',
    'RH_2_m_above_ground_1_hour_fcst',
    'DSWRF_surface_1_hour_fcst',
    'WIND_10_m_above_ground_1_hour_fcst'}



# %%
H = Herbie('2021-10-01 00:00', fhr=23, fxx=48, model='nbm', product='co')
vnames_all = get_variable_names(H)
short_names = {'DSWRF', 'RH', 'TMP', 'APCP', 'WIND'}
vnames_short = vnames_all[vnames_all['variable'].isin(short_names)]
print(vnames_short['name'])

date = '2021-10-01'

# identify the latest available forecast release time
fhr_avail = get_nbm_times(date)
fhr = fhr_avail[-1]

# identify the available forecast hours
fxx_avail = get_nbm_times(date, fhr)

# identify the timestep separating a given forecast hour with the previous
fxx_step = np.concatenate([np.array([1]), np.diff(np.array(fxx_avail))])

np.hstack([fxx_avail, fxx_step])
zip(fxx_avail, fxx_step)
res="\n".join("{} {}".format(x, y) for x, y in zip(fxx_avail, fxx_step))
print(res)

# %%
'''----------------- browse the NBM collection ----------------'''

# I'm using my own function instead of Herbie.xarray to open the GRIB files
# because there seems to be a bug in cfgrib (the module that Herbie uses)
# that causes problems with NBM. I think this has to do with non-uniqueness of
# variable names in the NBM GRIBs.

# define a herbie object, then copy the table of variable names
H = Herbie('2021-01-01', fhr=0, fxx=1, model='nbm', product='co')
vnames_all = get_variable_names(H)
print(vnames_all)

# We have a shorthand for this in the get_nbm function. When it is passed
# only a date, it calls get_variable_names on the Herbie object with
# fhr=0, fxx=1 and returns the table  
get_nbm('2021-01-01')

# note that some dates/times may only include subset of these 72 variables

# define a set of variables of interest
# based on: https://vlab.noaa.gov/web/mdl/nbm-wx-elements
short_names = {'DSWRF', 'RH', 'TMP', 'APCP', 'WIND'}



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
dates = pd.date_range('2021-11-06', '2021-11-08', freq='D')

# We will want to download all hours with fhr = [0, 1, ..., 23].
# forecast release time of day: integer or list with elements in 0,..23
fhr = list(range(0, 24))


# %%
'''----------------- Download and import data ----------------'''

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
# 

H = Herbie('2021-10-01 00:00', fhr=0, fxx=2, model='nbm', product='co')




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
# plot an example (first 24 hours in series for precip)
example_tseries = xdata_clipped['APCP_surface_01_hour_acc_fcst'][fhr, :, :]
example_tseries.plot(col='time', col_wrap=4)

# %%
# plot an example (first 24 hours in series for precip)
vname = list(vnames_toget)[4]
example_tseries = xdata_clipped[vname][fhr,:,:]
example_tseries.plot(col='time', col_wrap=4)