# for browsing AWS files
import fsspec

# for copying attributes dictionary and simplifying its values
import copy
import pathlib

# for opening GRIBs and converting to xarray
import rioxarray
import pygrib

# Herbie methods for downloading GRIBs
from herbie.archive import Herbie
from herbie.tools import bulk_download

# miscellaneous
import datetime
import numpy
import pandas as pd
import xarray


def get_variable_names(H, join=True):
    """
    Returns a list of variable names and regexp to use as searchString

    Herbie requests can be for a subset of variables in the file, but we
    need to know the variable names (and the appropriate regexp) ahead of
    time. This function produces a dataframe with this information.

    H can be a list of Herbie class objects, in which case the function
    returns the corresponding list of dataframes; or if join=True, a single
    table containing all unique variable names discovered.

    Arguments
    ----------
    H = a Herbie class object or list of them
    join = boolean, indicates to concatenate list output into single dataframe

    Return
    ----------
    A pandas dataframe or list of them

    """
    # handle list input
    if isinstance(H, list):
        out_df = [get_variable_names(h) for h in H]
        if join:
            out_df = pd.concat(out_df).drop_duplicates()
        return(out_df)

    # define relevant column names
    naming_keys = ['variable', 'level', 'forecast_time', 'prob']

    # fetch the index file data
    if hasattr(H, 'idx_df'):
        H_idx = copy.deepcopy(H.idx_df)
    else:
        H_idx = H.read_idx()

    # name the unnamed column
    H_idx = H_idx.rename({'?': 'prob'}, axis=1)

    # build layer names and regexp
    name_strings = H_idx[naming_keys].agg(' '.join, axis=1)
    regexp_strings = H_idx[naming_keys].agg(':'.join, axis=1)
    regexp_strings = regexp_strings.str.replace(':$', '', regex=True)

    # strip the names of special characters to avoid netCDF issues
    name_strings = name_strings.str.replace('[^A-Za-z0-9_ ]+', '', regex=True)
    name_strings = name_strings.str.strip()
    name_strings = name_strings.str.replace(' ', '_')

    # concantate into dataframe and return
    new_dict = {'name': name_strings, 'searchString': regexp_strings}
    out_df = pd.DataFrame(new_dict)
    out_df = pd.concat([out_df, H_idx[naming_keys]], axis=1)
    return(out_df)


def simplify_dict(D):
    """
    Simplify a dictionary by coercing various things to string.
    This is for preparing attributes for a netCDF file.

    Arguments
    ----------
    D = a dictionary (or list of them) to simplify

    Return
    ----------
    A dictionary of netCDF friendly attributes (or list of them)
    """
    for key in D:
        key_value = D[key]
        # coerce simple types: bool, datetime, ndarry, Path and Timestamp
        stype = (
            bool,
            datetime.datetime,
            numpy.ndarray,
            pathlib.WindowsPath,
            pd.Timestamp)
        if isinstance(key_value, stype):
            D[key] = str(key_value)
        # coerce NoneTypes
        if key_value is None:
            D[key] = 'empty'
        # coerce lists using delimiter ', '
        if isinstance(key_value, list):
            key_value = ', '.join(str(key_value))
        # coerce nested dictionaries using delimiter '='
        if isinstance(key_value, dict):
            dict_as_list = []
            for name in key_value:
                dict_as_list.append(name + '=' + key_value.get(name))
            D[key] = ', '.join(dict_as_list)

    return(D)


def H_dict(H):
    """
    Copy attributes from the Herbie object to netCDF friendly dictionary

    This function flattens and stringifies attributes that otherwise
    couldn't be stored as netCDF attributes. We use this to copy the
    metadata about the download into the output dataset

    Note that the dictionary is deep copied, so nothing is changed in
    the input Herbie object H

    Dictionary values of type bool, WindowsPath, Timestamp, None, and
    dict cause problems with backwards compatibility of the netCDF file
    written by xarray's to_netCDF. This functions converts these values
    to a string representation.

    When H is a list of Herbie objects (eg as the return value of
    Herbie.tools.bulk_download), the function returns a list of
    dictionaries.

    Arguments
    ----------
    H = a Herbie class object or list of them

    Return
    ----------
    A dictionary of attributes (or list of them)


    """
    # recursive call for list input
    if isinstance(H, list):
        return([H_dict(h) for h in H])

    # copy the metadata associated with the download (omit indexing table)
    herbie_dict = copy.deepcopy(H.__dict__)
    herbie_dict.__delitem__('idx_df')

    # reshape dictionary returned by Herbie to write as netcdf attributes
    return(simplify_dict(herbie_dict))


def merge_attributes(in_list, simplify=False):
    """
    Merge like attributes from a list of objects

    This function checks the attributes of each element of in_list against
    each of the attributes in in_list[0]. In cases where the attribute is
    different in one or more elements of in_list, the function replaces the
    attribute in in_list[0] with a string representation of all unique
    values of that attribute.

    Only the attributes of in_list[0] are modified

    Arguments
    ----------
    in_list = a list of objects with like attribute names
    simplify = boolean indicated to omit list-type attributes
    """
    # get attributes from first list entry
    attr_global = in_list[0].attrs

    # loop over attribute keys, then over all over list entries after the first
    attr_to_drop = []
    for attr in attr_global:
        for x in in_list[1:]:

            # skip if this attribute value is identical to the first
            if hasattr(x, attr):
                if str(x.attrs[attr]) != str(attr_global[attr]):

                    # turn the global attribute into a list if it isn't already
                    if not isinstance(attr_global[attr], list):
                        attr_global[attr] = [str(attr_global[attr])]

                    # add the new attribute to the list
                    attr_global[attr].append(str(x.attrs[attr]))
            else:
                attr_to_drop.append(attr)

        # remove or flatten lists into strings, depending on "simplify"
        if isinstance(attr_global[attr], list):
            if simplify:
                attr_to_drop.append(attr)
            else:
                attr_global[attr] = ', '.join(attr_global[attr])

    # pop in different loop so dictionary doesn't change size during loop
    for attr in attr_to_drop:
        if hasattr(x, attr):
            attr_global.pop(attr)


def H2x(H, times=[], simplify=True):
    """
    Open a GRIB file (or list of them) from a Herbie object

    This function opens a GRIB file downloaded by Herbie and returns
    an xarray dataset with global attributes copied from H and layer
    attribues copied from the GRIB file (these are read separately
    with pygrib)

    When H is a list of Herbie objects, the function returns a list of
    datasets. If "times" is also supplied, the function concatenates the
    datasets by creating a new dimension (time) and assigning the values
    of "times" as coordinates. The output is a single xarray dataset

    Arguments
    ----------
    H = a Herbie class object or list of them
    times = optional list of datetime values for new time coordinate
    simplify = boolean indicated to omit list-type attributes

    Return
    ----------
    An xarray dataset (or list of them)

    """

    # define GRIB keys to omit from layer attributes
    keys_omit = {
        'latLonValues',
        'latitudes',
        'longitudes',
        'distinctLatitudes',
        'distinctLongitudes',
        'values',
        'codedValues'}

    # recursive call for list argument H
    if isinstance(H, list):
        # call the function on every element of H
        out_list = [H2x(h) for h in H]
        ntimes = len(times)

        # return in list if times not supplied
        if ntimes == 0:
            return(out_list)

        # proceed only if times is the expected length
        if ntimes == len(H):
            # merge global attributes
            merge_attributes(out_list, simplify=simplify)

            # merge layer attributes
            varnames = list(out_list[0].data_vars.keys())
            for var in varnames:
                xarray_byvar = [x.data_vars[var] for x in out_list]
                merge_attributes(xarray_byvar, simplify=simplify)

            # reshape xarray with new dimension/coordinate "time" and return
            xarray_output = xarray.concat(out_list, 'time')
            return(xarray_output.assign_coords({'time': times}))

        # should have hit a return statement by now
        raise Exception('The length of times must match the length of H')

    # translate Herbie dictionary
    herbie_dict = H_dict(H)

    # open the index file as a table
    vnames = get_variable_names(H)

    # open the GRIB2 separately with pygrib and rioxarray
    grib_pyg = pygrib.open(str(H.local_grib_subset))
    grib_rio = rioxarray.open_rasterio(str(H.local_grib_subset))

    # drop attributes from the rioxarray object (we will replace them below)
    grib_rio.attrs = {}

    # loop over GRIB "messages" (layers) saving grid data and metadata
    grib_data_dict = {}

    # note that grib messages are 1-indexed
    for i in range(1, 1 + grib_pyg.messages):

        # ````load metadata using pygrib````

        # load message summary and use it to initiate dictionary
        msg = grib_pyg[i]
        pygrib_summary = msg.__str__()
        msg_dict = {'pygrib_band': i, 'pygrib_summary': pygrib_summary}

        # make a list of available/desired attribute key names
        keys_avail = set(msg.keys())
        keys = list(keys_avail.difference(keys_omit))

        # make a dictionary of name-value pairs for these attributes
        for key in keys:
            msg_dict[key] = getattr(msg, key)
        simplify_dict(msg_dict)

        # ````load raster data using pygrib````

        # (with 0-indexing) extract the data, dropping "band" dimension
        msg_data = grib_rio[i-1].drop('band')

        # copy the name and metadata from index file
        idx_dict = vnames.iloc[i-1].to_dict()
        lyrname = idx_dict['name']

        # build new dictionary entries to append
        idx_dict.pop('name')
        for key in idx_dict:
            idx_dict['herbie_' + key] = idx_dict.pop(key)

        # append the dictionary from above as new attributes for the layer
        idx_dict.update(msg_dict)
        msg_data.attrs = idx_dict

        # append to dictionaries in storage
        grib_data_dict.update({lyrname: msg_data})

    # compile as xarray Dataset
    xarray_output = xarray.Dataset(grib_data_dict)

    # this attribute required to preserve CRS info when saving/reloading
    herbie_dict.update({'coordinates': 'spatial_ref'})

    # add global attributes
    xarray_output.attrs = herbie_dict
    return(xarray_output)


def get_nbm(
    dates,
    vnames_toget={},
    fhr=list(range(0, 23+1)),
    fxx=1,
    product='co',
    simplify=True
):
    """
    Download/open a GRIB file (or list of them) using Herbie

    This function uses Herbie's bulk_download to get GRIBS for a range of
    of forecast release hours on a given day, then uses H2x to open the
    data and reshape as a single, well-labelled xarray object.

    For a fixed forecast hour (fxx), the function flattens all requested
    forecast datetimes into a new dimension "time", and adds a bunch of
    global and layer-wise attributes

    When vnames_toget is empty (default), the function instead returns a
    pandas dataframe of variable name options for the first requested
    datetime.

    Arguments
    ----------
    date_string = the requested forecast release date, as a string
    fhr = a list of forecast release hours to download
    vnames_toget = a set of requested variable name strings
    fxx = (integer) forecast hour, one of 1, ... 23
    product = one of 'ak', 'co', etc
    simplify = boolean indicated to omit list-type attributes

    Return
    ----------
    An xarray dataset

    """
    # coerce forecast release hours (time of day) to list
    fhr = [fhr] if not isinstance(fhr, list) else fhr

    # handle lists of dates
    if isinstance(dates, (list, pd.DatetimeIndex)):
        # recursive calls looping over dates to build list of datasets
        xdata_list = []
        for date in dates:
            print(f"processing date: {date.date()} ({len(fhr)} files)...")
            xout = get_nbm(date, vnames_toget, fhr, fxx, product, simplify)
            xdata_list.append(xout)

        # merge global attributes
        print(f"merging data for {len(dates)} dates...")
        merge_attributes(xdata_list, simplify=simplify)

        # merge layer attributes
        varnames = list(xdata_list[0].data_vars.keys())
        for var in varnames:
            xarray_byvar = [x.data_vars[var] for x in xdata_list]
            merge_attributes(xarray_byvar, simplify=simplify)

        # concatenate datasets along time
        xout = xarray.concat(xdata_list, 'time')
        return(xout)
    # single date requests handled below...

    # express input date/times as datetime lists
    date = pd.to_datetime(dates)
    rdates = [date + datetime.timedelta(hours=hr) for hr in fhr]
    fdates = [date + datetime.timedelta(hours=hr + fxx) for hr in fhr]

    # fetch a table of variable names for the first forecast time
    date_string = rdates[0].__str__()
    H_probe = Herbie(date_string, fxx=fxx, model='nbm', product=product)
    vnames_all = get_variable_names(H_probe)

    # handle variable name requests
    if len(vnames_toget) == 0:
        return(vnames_all)

    # define a regexp for Herbie to request these layers only
    vnames = vnames_all[vnames_all['name'].isin(vnames_toget)]
    search_words = vnames['searchString']
    search_expr = '^' + '$|^'.join(search_words) + '$'

    # download the day's GRIB2 files using Herbie package
    H = bulk_download(
        rdates,
        search_expr,
        fxx=fxx,
        model='nbm',
        product='co')

    # open as xarray, creating time dimensions/coordinate
    return(H2x(H, times=fdates, simplify=simplify))