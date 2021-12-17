# Scan AWS bucket for NOAA's National Blended Model (NBM) for forecast gribs
#
# MOTIVATION
#
# I initially assumed the NBM forecast files would follow a straightforward
# and regular pattern forecast hours (fxx) and temporal levels for the
# variables. This doesn't seem to be the case, so before moving ahead with code
# to download the full archive it will be a good idea to make a list of all
# available date/time/fxx values in the collection.
#
# DESCRIPTION
#
# The script searches subdirectories of the bucket "noaa-nbm-grib2-pds" for
# GRIB files having an associated index file. These are used to make a list of
# available forecast release date-times and forecast hours, which is saved as
# CSV for later reference.
#
# TODO: modify the script so we can run it any time to update the CSV

import datetime
import numpy
import pandas as pd
import matplotlib
import seaborn

# Herbie package for downloading GRIBs and my own extensions
from helpers import get_nbm_times, get_nbm_variable_names, get_nbm

# %%
'''----------------- browse the NBM collection ----------------'''

# The get_nbm_times function looks in the AWS bucket for directories and files
# corresponding to grib2 forecast data. Call the function without arguments to
# get a list of dates for which there is a data directory
dates = get_nbm_times()
# Not all of these are accessible! More on this in next section...

# convert to Timestamps
dates = pd.to_datetime(dates).sort_values()

# We find almost a complete sequence of daily folders going back to May 18,
# 2020, except for three missing dates in October of 2020
start_date = min(dates)
end_date = max(dates)
dates_expected = pd.date_range(start_date, end_date, freq='D')
print(dates_expected.difference(dates))

# %%
# Within each of the daily folders, there are subfolders for each release time.
# Pass a date to get_nbm_times to get an array of release hours (fhr) for which
# the subfolder could be found
example_date = pd.to_datetime('2020-09-30')
print(get_nbm_times(example_date))

# %%
# The above calls only check for directories, not the data files themselves.
# To import the data, we will need both the grib2 file and its index file
# (*grib2.idx). Check for these by passing both a date and a time to
# get_nbm_times
print(get_nbm_times(date=example_date, fhr=0))

# The function returns an array of available forecast prediction hours (fxx)
# at supplied date/time. Each of these corresponds to exactly one GRIB file.
# Notice that date and hour (fhr) are passed separately. Any embedded time
# in argument "date" (eg Timestamp type) is ignored.

# %%
# Once you've found a file (date/fhr/fxx) of interest you can check the
# available layers in the file using get_nbm_variable_names
print(get_nbm_variable_names({'date': example_date, 'fhr': 0, 'fxx': 1}))

# This takes a moment to complete because the function uses Herbie to
# download the index data (a small ASCII file). Notice that the variable
# names depend on temporal level, among other things, meaning that you will
# have to specify different variable names to access a given variable at
# different fxx hours.

# %%
'''----------------- Compile available dates ----------------'''

# %%
# Not all of the dates returned by get_nbm_times() are accessible to us, as
# we need the grib2 index files, which are unavailable for early dates
# in the collection. Another hiccup is that the available set of forecast
# prediction hours (fxx) varies depending on the release hour (fhr).
#
# For some release times (usually, fhr = 1, 7, 13, 19) we see hourly forecasts
# up to 11 days (fxx in 1, 2, ..., 264). But for most others we see a schedule
# of hourly forecasts for the first 36 hours, then 3-hourly until day 8, then
# 6-hourly for the remainder. There are also (rarely) unexpected gaps within
# this framework, such as entire missing days.
#
# To get a better view of what's going on and plan a strategy for filling gaps,
# we start by building a complete list of available files, ie which fxx values
# are posted on which release timestamps. The code below calls get_nbm_times
# in a nested loop (over dates, then fhr) and reshapes the result as a
# Timestamp-indexed dataframe.

# build a list of info about forecast accessibility
# This may take a few hours to complete:
fxx_expected = range(1, 264+1)
fxx_names = ['fxx_' + str(fxx) for fxx in fxx_expected]
fxx_loop_result = []
for date in dates:

    # scan for subfolders representing forecast release times
    fhr_avail = get_nbm_times(date)

    # loop over forecast release times, scanning contents of subfolder
    for fhr in fhr_avail:

        # scan for available grib2 index files and print a message
        fxx_avail = get_nbm_times(date, fhr)
        n_files = len(fxx_avail)
        datetime_avail = date.replace(hour=fhr)
        print(str(datetime_avail) + ': ' + str(n_files) + ' files found')

        # compare to expected (hourly) fxx to get dataframe of booleans
        has_fxx = numpy.isin(fxx_expected, fxx_avail)
        df_out = pd.DataFrame(has_fxx, fxx_names, columns=[datetime_avail])

        # transpose and add to storage list
        fxx_loop_result.append(df_out.transpose())

# define an hourly sequence of release times from start to end
start_datetime = start_date.replace(hour=min(get_nbm_times(start_date)))
end_datetime = end_date.replace(hour=max(get_nbm_times(end_date)))
datetimes_expected = pd.date_range(start_datetime, end_datetime, freq='h')

# compile into a single dataframe, filling NaNs on missing Timestamps
fxx_all = pd.concat(fxx_loop_result, axis=0).reindex(datetimes_expected)

# replace NaN gaps with 'False', invoke copy to defragment
fxx_all = fxx_all.fillna(False).copy()

# add a timestamp indicating when this information was last checked
date_today = datetime.datetime.utcnow().date()
fxx_all.insert(loc=0, column='UTC_last_checked', value=date_today)

# %%
'''----------------- Save to disk as CSV ----------------'''

# save as CSV
csv_parent = pathlib.Path('D:/UYRW_data/side_projects/nwp_ingestion/data')
csv_path = pathlib.Path(csv_parent, 'available_fxx.csv')
fxx_all.to_csv(csv_path, index_label='UTC_release_date')

# reload CSV
fxx_all = pd.read_csv(csv_path)

# %%
'''----------------- Visualization ----------------'''

# select a week of dates to display
start_datetime_test = '2021-07-01 01:00'
end_datetime_test = '2021-07-03 00:00'
dates_test = pd.date_range(start_datetime_test, end_datetime_test, freq='h')

# build the plotting dataframe, simplifying axis tick labels
fxx_toplot = fxx_all.loc[dates_test].drop(['UTC_last_checked'], axis=1)
fxx_toplot = fxx_toplot.rename(columns=lambda x: int(x.replace('fxx_', '')))
simple_dates = fxx_toplot.index.map(lambda x: x.strftime('%Y-%m-%d %H'))
fxx_toplot.set_index(simple_dates, inplace=True)

# define a large plot window so we can draw all the dates/times
matplotlib.pyplot.subplots(figsize=(20, 10))

# plot the data as a heatmap
fig = seaborn.heatmap(
    fxx_toplot,
    linewidths=0.5,
    cmap='GnBu',
    cbar=False)

# set axis labels and title
title_string = 'NBM forecast files on AWS (checked ' + str(date_today) + ')'
xstring = 'UTC prediction time (in hours ahead of release time)'
ystring = 'UTC forecast release time as <date hour>'
fig.set(xlabel=xstring, ylabel='', title=tstring)

# save as PNG
fig.figure.savefig('fxx_sample.png', format='png', dpi=300, facecolor='white')

# %%
