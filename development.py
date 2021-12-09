import fsspec

def get_times(date=None, fhr=None, product='co'):
    fs = fsspec.filesystem('s3', anon=True)
    

    bucket_prefix = 'noaa-nbm-grib2-pds/blend.'
    #fs.glob(bucket_prefix + '*/00/core')

    # find valid dates
    if date is None:
        dirs_list = fs.glob(bucket_prefix + '*')
        dstring_list = [d.replace(bucket_prefix, '') for d in dirs_list]
        dlist = [pd.to_datetime(d).strftime('%Y-%m-%d') for d in dstring_list]
        return(dlist)

    # find valid forecast release time
    dstring = pandas.to_datetime(date).strftime('%y%m%d')
    if fhr is None:
        search_root = bucket_prefix + dstring + '/'
        dirs_list = fs.glob(search_root + '*')
        fhr_list = [int(d.replace(search_root, '')) for d in dirs_list]
        return(fhr_list.list)

    # find valid forecast hours
    fhr = int(fhr)
    search_root = bucket_prefix + dstring + '/' + str(fhr).zfill(2) + '/core/'
    file_prefix = 'blend.t' + str(fhr).zfill(2) + 'z.core.f'
    search_suffix = '.' + product + '.grib2.idx'
    files_list = fs.glob(search_root + file_prefix + '*' + search_suffix)
    idx_list = [int(f.replace(search_suffix, '')[-3:]) for f in files_list]
    return(idx_list)