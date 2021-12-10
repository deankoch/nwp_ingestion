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


def get_nbm_times(date=None, fhr=None, product='co'):
    """
    Scan the AWS collection for filenames, returning the available dates,
    forecast release times, or forecast hours (depending on the arguments
    dates and fhr)

    If date is not supplied, the function returns all dates for
    which there is a subdirectory in the AWS bucket (as an array).

    If date but not fhr is supplied, the function returns an array of integers
    corresponding to the forecast time subdirectories found in /blend.{date}

    If both date and fhr are supplied, the function returns an array of
    integers corresponding to all forecast hours for which there is a
    *.grib2.idx (GRIB index) file. These files should be accesible using
    Herbie.

    Arguments
    ----------
    date = any date-like object understood by pd.to_datetime
    fhr = any object coercible to an integer by int(). Must be in [0,...23]
    product = one of 'co', 'ak', etc

    Return
    ----------
    A numpy array of date strings or integers representing hours (see above)
    """
    # define the bucket name and a directory prefix to search for
    bucket_prefix = 'noaa-nbm-grib2-pds' + '/blend.'

    # initialize the s3 protocol object
    fs = fsspec.filesystem('s3', anon=True)

    # find valid dates
    if date is None:
        dirs_list = fs.glob(bucket_prefix + '*')
        dstring_list = [d.replace(bucket_prefix, '') for d in dirs_list]
        dlist = [pd.to_datetime(d).strftime('%Y-%m-%d') for d in dstring_list]
        return(numpy.array(dlist))

    # find valid forecast release time
    dstring = pd.to_datetime(date).strftime('%Y%m%d')
    if fhr is None:
        search_root = bucket_prefix + dstring + '/'
        dirs_list = fs.glob(search_root + '*')
        fhr_list = [int(d.replace(search_root, '')) for d in dirs_list]
        return(numpy.array(fhr_list))

    # find valid forecast hours
    fhr = int(fhr)
    search_root = bucket_prefix + dstring + '/' + str(fhr).zfill(2) + '/core/'
    file_prefix = 'blend.t' + str(fhr).zfill(2) + 'z.core.f'
    search_suffix = '.' + product + '.grib2.idx'
    files_list = fs.glob(search_root + file_prefix + '*' + search_suffix)
    idx_list = [int(f.replace(search_suffix, '')[-3:]) for f in files_list]
    return(numpy.array(idx_list))


def get_nbm_variable_names(H, join=True):
    """
    Returns a list of variable names to use with a Herbie object

    Herbie requests can be for a subset of variables in the file, but we
    need to know the variable names (and the appropriate regexp) ahead of
    time. This function produces a dataframe with this information.

    Specify the GRIB2 file of interest by passing a Herbie class object or a
    list of them. In the case of list input, the function returns the
    corresponding list of dataframes - or, if join=True, a single table
    containing all unique variable names discovered.

    Alternatively, a dictionary (or list of them) containing key "date", and
    optionally "fhr", "fxx" (with default 1), "product" (with default 'co')
    can be passed as H, and the function will create the Herbie object
    internally. If "fhr" is not supplied, the function attempts to extract it
    from "date" or otherwise sets default fhr=0

    Arguments
    ----------
    H = a Herbie class object or dictionary, or list of them
    join = boolean, indicates to concatenate list output into single dataframe

    Return
    ----------
    A pandas dataframe or list of them

    """
    # handle list input
    if isinstance(H, list):
        out_df = [get_nbm_variable_names(h) for h in H]
        if join:
            out_df = pd.concat(out_df).drop_duplicates()
        return(out_df)

    # handle dictionary input
    if isinstance(H, dict):

        # check for invalid input and create datetime object from date
        if 'date' not in H: 
            raise Exception('dictionary H must contain key "date"')
        date = pd.to_datetime(H['date'])

        # set defaults and coerce hours to integers
        H['fxx'] = 1 if 'fxx' not in H else int(H['fxx'])
        H['product'] = 'co' if 'product' not in H else H['product']   
        H['fhr'] = date.time().hour if 'fhr' not in H else int(H['fhr'])

        # define date string for Herbie call and instantiate
        hour_string = str(H['fhr']).zfill(2) + ':00'
        date_string = date.strftime('%Y-%m-%d') + ' ' + hour_string
        H_out = Herbie(
            date_string,
            fxx=H['fxx'],
            model='nbm',
            product=H['product'])

        # raise an error if the file is not found, otherwise do recursive call
        if H_out.grib is None:
            t_info = date_string + ', fxx=' + str(H['fxx'])
            p_info = ' (' + H['product'] + ')'
            raise Exception('file not found for input ' + t_info + p_info)
        return(get_nbm_variable_names(H_out))

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
    vnames = get_nbm_variable_names(H)

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
    vnames_all = get_nbm_variable_names(H_probe)

    # handle variable name requests
    if len(vnames_toget) == 0:
        print('returning dataframe of variable names')
        return(vnames_all)

    # define a regexp for Herbie to request these layers only
    vnames = vnames_all[vnames_all['name'].isin(vnames_toget)]
    search_words = vnames['searchString']
    search_expr = '^' + '$|^'.join(search_words) + '$'

    # download the day's GRIB2 files using Herbie package
    print('accessing GRIB files using Herbie.tools.bulk_download')
    H = bulk_download(
        rdates,
        search_expr,
        fxx=fxx,
        model='nbm',
        product='co')

    # open as xarray, creating time dimensions/coordinate
    print('opening GRIB files and building xarray')
    return(H2x(H, times=fdates, simplify=simplify))
