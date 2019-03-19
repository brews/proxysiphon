"""Read proxies netCDF5 file and write proxies HDF5 file for LMR Data Assimilation workflow
"""

import os
import logging
import numpy as np
import pandas as pd
import xarray as xr
import netCDF4
import shapely.affinity
import shapely.geometry

import erebusfall as ef


log = logging.getLogger(__name__)


class DistanceThresholdError(Exception):
    """Raised when the distance between two points is further than a threshold

    Parameters
    ----------
    target_distance : int or float
        The distance between two target points (km).
    distance_threshold : int or float
        The distance threshold.

    """

    def __init__(self, target_distance, distance_threshold):
        self.target_distance = target_distance
        self.distance_threshold = distance_threshold


def chord_distance(latlon1, latlon2):
    """Chordal distance between two sequences of (lat, lon) points

    Parameters
    ----------
    latlon1 : sequence of tuples
        (latitude, longitude) for one set of points.
    latlon2 : sequence of tuples
        A sequence of (latitude, longitude) for another set of points.

    Returns
    -------
    dists : 2d array
        An mxn array of Earth chordal distances [1]_ (km) between points in
        latlon1 and latlon2.

    References
    ----------
    .. [1] https://en.wikipedia.org/wiki/Chord_(geometry)

    """
    earth_radius = 6378.137  # in km

    latlon1 = np.atleast_2d(latlon1)
    latlon2 = np.atleast_2d(latlon2)

    n = latlon1.shape[0]
    m = latlon2.shape[0]

    paired = np.hstack((np.kron(latlon1, np.ones((m, 1))),
                        np.kron(np.ones((n, 1)), latlon2)))

    latdif = np.deg2rad(paired[:, 0] - paired[:, 2])
    londif = np.deg2rad(paired[:, 1] - paired[:, 3])

    a = np.sin(latdif / 2) ** 2
    b = np.cos(np.deg2rad(paired[:, 0]))
    c = np.cos(np.deg2rad(paired[:, 2]))
    d = np.sin(np.abs(londif) / 2) ** 2

    half_angles = np.arcsin(np.sqrt(a + b * c * d))

    dists = 2 * earth_radius * np.sin(half_angles)

    return dists.reshape(m, n)


def get_nearest(latlon, dain, depth=None, lat_coord='Latitude', lon_coord='Longitude',
                depth_coord='depth', distance_threshold=1500):
    """Get nearest non-NaN to latlon from xarray.DataArray obj

    Finds the nearest not NaN to latlon, and optionally depth. It searches for a
    nearest depth first, if given, and then searches for nearest latlon. Note
    that this does not work with irregular grids, such as rotated polar, etc.

    Parameters
    ----------
    latlon : tuple
        Target latitude and longitude. Must be within -90 to 90 and -180 to 180.
    dain : xarray.DataArray
        Field with regular latlon coordinates.
    depth : float or int, optional
        Target depth to get nearest.
    lat_coord : str, optional
        Name of the latitude coordinate in ``da``.
    lon_coord : str, optional
        Name of the longitude coordinate in ``da``.
    depth_coord : str, optional
        Name of the depth coordinate in ``da``.
    distance_threshold : float or int, optional
        If the nearest distance is larger than this, raise

    Returns
    -------
    nearest : xarray.DataArray
        Nearest points.
    nearest_distance : float
        Chordal distance (km) from target to matched gridpoint.

    Raises
    ------
    DistanceThresholdError
    """
    da = dain.copy()

    assert latlon[0] <= 90 and latlon[0] >= -90
    assert latlon[1] <= 180 and latlon[1] >= -180

    assert lat_coord in da.coords
    assert lon_coord in da.coords

    assert (da[lat_coord].ndim == 1) and (da[lon_coord].ndim == 1)

    # First, find the nearest depth index, if given.
    if depth is not None:
        assert depth_coord in da.coords

        # Note use 'pad' because want next upper water column level value.
        da = da.sortby('depth')
        da = da.sel(**{depth_coord: depth}, method='nearest')

    # Now search for nearest latlon point.
    da_stack = da.stack(yx=[lat_coord, lon_coord]).dropna('yx')
    da_latlon_stack = np.vstack((da_stack[lat_coord], da_stack[lon_coord])).T

    # Any values above 180 become negative -- needed for 0-360 longitudes.
    highlon_msk = da_latlon_stack > 180
    da_latlon_stack[highlon_msk] = da_latlon_stack[highlon_msk] - 360

    distance = chord_distance(np.array([latlon]), da_latlon_stack)
    nearest = da_stack.isel(yx=np.argmin(distance))
    nearest_distance = np.min(distance)

    if nearest_distance > distance_threshold:
        raise DistanceThresholdError(nearest_distance, distance_threshold)

    return nearest, nearest_distance


def get_netcdf_resource(fl, **kwargs):
    """Read NetCDF files as package resource, output for xarray.Dataset

    Parameters
    ----------
    fl : str
        NetCDF resource name.
    **kwargs :
        Passed to ``xarray.open_dataset``.

    Returns
    -------
    data : xarray.Dataset
    """
    here = os.path.abspath(os.path.dirname(__file__))
    flpath = os.path.join(here, fl)
    data = xr.open_dataset(flpath, **kwargs)
    return data


def poly_dateline_wrap(p):
    """Split dateline crossing polygon into multipoly so wraps.

    Parameters
    ----------
    p : shapely.geometry.Polygon

    Returns
    -------
    out : shapely.geometry.Polygon or shapely.geometry.MultiPolygon
    """
    dateline = shapely.geometry.LineString([(180, 90), (180, -90)])

    if not dateline.crosses(p):
        return p

    right_ls = shapely.geometry.LineString([(360, 90), (360, -90)])
    left_ls = shapely.geometry.LineString([(-180, 90), (-180, -90)])
    left_clip = shapely.geometry.MultiLineString([left_ls, dateline]).convex_hull
    right_clip = shapely.geometry.MultiLineString([dateline, right_ls]).convex_hull

    northpacific_left = p.intersection(left_clip)
    dl_right_off = p.intersection(right_clip)

    # Shift right of dateline back so lon (x) is within [-180, 180]
    northpacific_lr = shapely.affinity.translate(shapely.geometry.LinearRing(dl_right_off.exterior.coords),
                                                 xoff=-360)
    northpacific_right = shapely.geometry.Polygon(northpacific_lr)
    out = shapely.geometry.MultiPolygon([northpacific_left, northpacific_right])
    return out


def find_seasonality(sitegrp, proxy_grp):
    """Return string list of ints giving site variable seasonality.

    Parameters
    ----------
    sitegrp : netCDF4.Group
        Proxy site netcdf group.
    proxy_grp: netCDF4.Group
        Proxy variable netcdf group

    Returns
    -------
    out : str
    """
    out = list(range(1, 13))

    proxy_type = str(proxy_grp.long_name)
    if proxy_type == "UK'37":
        log.debug('finding seasonality for UK37 proxy')
        # Jess Tierney polygons from BAYSPLINE
        mediterranean = shapely.geometry.Polygon([(-5.5, 36.25),
                                                  (3, 47.5),
                                                  (45, 47.5),
                                                  (45, 30),
                                                  (-5.5, 30)])
        northatlantic = shapely.geometry.Polygon([(-55, 48),
                                                  (-50, 70),
                                                  (20, 70),
                                                  (10, 62.5),
                                                  (-4.5, 58.2),
                                                  (-4.5, 48)])
        northpacific_raw = shapely.geometry.Polygon([(135, 45),
                                                     (135, 70),
                                                     (250, 70),
                                                     (232, 52.5),
                                                     (180, 45)])
        northpacific = poly_dateline_wrap(northpacific_raw)

        latlon = (float(sitegrp.latitude), float(sitegrp.longitude))
        assert -90 < latlon[0] < 90, 'site latitude must be -90 < lat < 90'
        assert -180 < latlon[1] < 180, 'site longitude must be -180 < lon < 180'

        site_location = shapely.geometry.Point(latlon[::-1])
        if mediterranean.contains(site_location):
            out = [1, 2, 3, 4, 5, 11, 12]
        if northatlantic.contains(site_location):
            out = [8, 9, 10]
        if northpacific.contains(site_location):
            out = [6, 7, 8]

    elif proxy_type == 'd18O':
        foram_type = str(proxy_grp.foraminifera_type)
        latlon = (float(sitegrp.latitude), float(sitegrp.longitude))
        log.debug('finding seasonality for d18O proxy ({}) @ {}'.format(foram_type, latlon))
        foram_seas = get_netcdf_resource('foram_seasons.nc', group='d18oc')
        foram_seas = foram_seas.where(foram_seas != '')  # Replace '' with nan.
        # Need to normalize full species/subspecies names to the variable names used
        # in the FORAM_SEASONS_NETCDF group.
        foraminifera_map = {'Globigerina bulloides': 'G. bulloides',
                            'Neogloboquadrina pachyderma sinistral': 'N. pachyderma',
                            'Neogloboquadrina incompta': 'N. incompta',
                            'Globigerinoides ruber white': 'G. ruber',
                            'Globigerinoides ruber pink': 'G. ruber',
                            'Trilobatus sacculifer': 'T. sacculifer',
                            }
        foram_spp_str = foraminifera_map.get(foram_type)
        if foram_spp_str is not None:
            nearest_seas, _ = get_nearest(latlon, foram_seas[foram_spp_str],
                                          lat_coord='lat', lon_coord='lon')
            out = nearest_seas.item()

        else:
            pass  # foram type not in our species mapping so use annual.

    elif proxy_type == 'Mg/Ca':
        foram_type = str(proxy_grp.foraminifera_type)
        latlon = (float(sitegrp.latitude), float(sitegrp.longitude))
        log.debug('finding seasonality for Mg/Ca proxy ({}) @ {}'.format(foram_type, latlon))
        foram_seas = get_netcdf_resource('foram_seasons.nc', group='mgca')
        foram_seas = foram_seas.where(foram_seas != '')  # Replace '' with nan.
        # Need to normalize full species/subspecies names to the variable names used
        # in the FORAM_SEASONS_NETCDF group.
        foraminifera_map = {'Globigerina bulloides': 'G. bulloides',
                            'Neogloboquadrina pachyderma sinistral': 'N. pachyderma sinistral',
                            'Neogloboquadrina incompta': 'N. pachyderma dextral',
                            'Globigerinoides ruber white': 'G. ruber white',
                            'Globigerinoides ruber pink': 'G. ruber pink',
                            'Trilobatus sacculifer': 'G. sacculifer',
                            }
        foram_spp_str = foraminifera_map.get(foram_type)
        if foram_spp_str is not None:
            nearest_seas, _ = get_nearest(latlon, foram_seas[foram_spp_str],
                                          lat_coord='lat', lon_coord='lon')
            out = nearest_seas.item()
        else:
            pass  # foram type not in our species mapping so use annual.
    log.debug('season months are: {}'.format(out))
    return str(out)


def icevol_correction(proxies, metadata):
    """Apply ice-volume correction to d18O proxies in LMR proxy dataframes.

    Returns modified copies of the original data.
    """
    p = proxies.copy()
    m = metadata.copy()

    matches = ['d18o' in c.lower() for c in p.columns]

    if not any(matches):
        return p, m

    m.set_index('Proxy ID', inplace=True)

    matches_idx = np.where(matches)
    matched_columns = p.columns[matches_idx]

    for c in matched_columns:
        proxy_raw = proxies[c].values
        age_raw = proxies[c].index.values
        age_yr = 1950 - age_raw  # CE/BC to Yr BP
        proxy_adjusted = ef.icevol_correction(age_yr, proxy_raw, proxytype='d18o',
                                              timeunit='ya')
        p.loc[:, c] = proxy_adjusted

        m.loc[c, 'Oldest (C.E.)'] = p[c].dropna().index.min()
        m.loc[c, 'Youngest (C.E.)'] = p[c].dropna().index.max()

    m.reset_index(inplace=True)
    return p, m


def lmr_da_dfs(sitegrp=None, agemodel_iter=None, find_modern_seasonality=True):
    """Return proxy data and metadata pandas df needed for LMR DA proxy input

    If no args are passed, return empty proxy and metadata (template) pandas
    dataframes.

    Parameters
    ----------
    sitegrp : netCDF4.Group or None, optional
        Proxy site netcdf group.
    agemodel_iter : int or None, optional
        Age-depth model iteration to use, if available. This is the index of the
        'agemodel_ensemble' column to use for the output data. If ``None`` (default),
        uses median age.
    find_modern_seasonality : bool
        Do you want to estimate sample seasonality using the modern record?
        Sets annual seasonality if False. Default is True.

    Returns
    -------
    data : pandas.DataFrame
    meta : pandas.DataFrame
    """
    data_template = pd.DataFrame()
    meta_template = pd.DataFrame({'Proxy ID': [], 'Site': [], 'Lat (N)': [], 'Lon (E)': [],
                                  'Archive type': [], 'Proxy measurement': [],
                                  'Resolution (yr)': [], 'Reference': [], 'Databases': [],
                                  'Seasonality': [], 'Elev': [],
                                  'Oldest (C.E.)': [], 'Youngest (C.E.)': []})
    meta_template = meta_template.astype(dtype={'Proxy ID': 'object', 'Site': 'object',
                                                'Lat (N)': 'float64', 'Lon (E)': 'float64',
                                                'Archive type': 'object',
                                                'Proxy measurement': 'object',
                                                'Resolution (yr)': 'float64',
                                                'Reference': 'object', 'Databases': 'object',
                                                'Seasonality': 'object', 'Elev': 'float64',
                                                'Oldest (C.E.)': 'float64',
                                                'Youngest (C.E.)': 'float64'})

    if sitegrp is None:
        return data_template.copy(), meta_template.copy()

    log.info('extracting LMR data from {}'.format(str(sitegrp.site_name)))

    variables_to_skip = ['depth', 'age_original', 'age_median', 'age_ensemble']
    proxy_variables = [(k, v) for k, v in sitegrp['data'].variables.items() if k not in variables_to_skip]

    latitude = float(sitegrp.latitude)
    longitude = float(sitegrp.longitude)
    elevation = float(sitegrp.elevation)
    data_df = data_template.copy()
    meta_df = meta_template.copy()

    for k, v in proxy_variables:
        log.debug('processing variable {}'.format(str(k)))

        # Convert years BP to CE/BP and find min max.
        if agemodel_iter is None:
            try:
                age_yrs_ce = 1950 - sitegrp['data'].variables['age_median'][:]
            except KeyError:
                age_yrs_ce = 1950 - sitegrp['data'].variables['age_original'][:]
        else:
            idx = int(agemodel_iter)
            age_yrs_ce = 1950 - sitegrp['data'].variables['age_ensemble'][:, idx]

        is_nan = np.isnan(v[:].filled(np.nan))
        youngest_ce = np.max(age_yrs_ce[~is_nan])
        oldest_ce = np.min(age_yrs_ce[~is_nan])

        # Put together proxy ID and proxy measurement strings.
        siteid = str(sitegrp.site_name).strip().lower()
        pmeasurement = str(k).strip().lower()
        # Append cleaning protocol info if available for Mg/Ca
        try:
            if v.long_name == 'Mg/Ca':
                log.debug('Mg/Ca variable, attempting to find cleaning protocol')
                # Should throw error if have Mg/Ca record w/o cleaning protocol info.
                cleaning_protocol = str(v.mgca_cleaning_protocol)
                cleaning_str = None
                if 'reductive' in cleaning_protocol.lower():
                    cleaning_str = ':reductive'
                elif 'barker' in cleaning_protocol.lower():
                    cleaning_str = ':barker'
                pmeasurement += cleaning_str
        except AttributeError:  # caught if v doesn't have long_name attrib...
            pass
        proxyid = siteid + ':' + pmeasurement

        this_meta = meta_template.copy()
        this_meta['Proxy ID'] = [proxyid]
        this_meta['Site'] = [siteid]
        this_meta['Lat (N)'] = latitude
        this_meta['Lon (E)'] = longitude
        this_meta['Archive type'] = ['Marine sediments']
        this_meta['Proxy measurement'] = [pmeasurement]
        this_meta['Resolution (yr)'] = [(youngest_ce - oldest_ce) / len(age_yrs_ce[~is_nan])]
        this_meta['Reference'] = [str(None)]
        this_meta['Databases'] = ['[DTDA]']
        if find_modern_seasonality:
            this_meta['Seasonality'] = find_seasonality(sitegrp, v)
        else:
            this_meta['Seasonality'] = str(list(range(1, 13)))
        this_meta['Elev'] = elevation
        this_meta['Oldest (C.E.)'] = oldest_ce
        this_meta['Youngest (C.E.)'] = youngest_ce
        meta_df = meta_df.append(this_meta, ignore_index=True)

        d = (pd.DataFrame({'Year C.E.': age_yrs_ce, proxyid: v[:].filled(np.nan)})
               .set_index('Year C.E.')
               .dropna(how='any'))
        data_df = data_df.join(d, how='outer')

    data_df = data_df.sort_index(ascending=False)

    return data_df, meta_df


def _lmr_df_from_nc_sites(fl, agemodel_iter=None, find_modern_seasonality=True):
    """Create LMR data and metadata dataframes from opened netCDF file group"""
    all_data_df, all_meta_df = lmr_da_dfs()

    for site_grp in fl.groups.values():
        try:
            site_data_df, site_meta_df = lmr_da_dfs(site_grp,
                                                    agemodel_iter=agemodel_iter,
                                                    find_modern_seasonality=find_modern_seasonality)
        except TypeError as e:
            errormsg = '{} raised - skipping file - {} - Mg/Ca variable likely missing cleaning protocol info'
            log.error(errormsg.format(e, site_grp.site_name))
            continue

        all_meta_df = all_meta_df.append(site_meta_df, ignore_index=True)
        all_data_df = all_data_df.join(site_data_df, how='outer')

    return all_meta_df, all_data_df


def nc2lmrdf(path_or_buffer, agemodel_iter=None,
             icevol_cor=True, find_modern_seasonality=True):
    """Read proxy netCDF and output to LMR DA-format dataframes.

    Parameters
    ----------
    path_or_buffer : str or netCDF4.Dataset
        Input proxy netCDF file path or opened buffer.
    agemodel_iter : int or None, optional
        Optional index of age-model iteration to use in output -- if multiple
        age-model iterations are available.
    icevol_cor : bool, optional
        Do you want to apply ice-volume correction to the d18O foraminiferal
        records? This is done with the `erebusfalls` package. `True` by default.
    find_modern_seasonality : bool, optional
        Do you want to estimate sample seasonality using the modern record?
        Sets annual seasonality if False. Default is True.

    Returns
    -------
    all_data_df : pd.Dataframe
        LMR proxy data.
    all_meta_df : pd.Dataframe
        LMR proxy metadata.

    Also see
    --------
    `proxysiphon.nc2lmrh5`

    """
    # Open proxy site netCDF file and collect needed data or read directly from
    # obj if it's already open.
    if isinstance(path_or_buffer, str):
        with netCDF4.Dataset(filename=path_or_buffer, mode='r') as fl:
            all_meta_df, all_data_df = _lmr_df_from_nc_sites(fl,
                                                             agemodel_iter=agemodel_iter,
                                                             find_modern_seasonality=find_modern_seasonality)
    else:
        all_meta_df, all_data_df = _lmr_df_from_nc_sites(path_or_buffer,
                                                         agemodel_iter=agemodel_iter,
                                                         find_modern_seasonality=find_modern_seasonality)

    all_data_df = all_data_df.sort_index(ascending=False)

    if icevol_cor:
        # Ice-volume correction to d18O foram proxies:
        all_data_df, all_meta_df = icevol_correction(all_data_df, all_meta_df)

    return all_data_df, all_meta_df


def nc2lmrh5(path_or_buffer, h5file, agemodel_iter=None, icevol_cor=True,
             find_modern_seasonality=True):
    """Read proxy netCDF and output to LMR DA-format HDF5 file.

    Parameters
    ----------
    path_or_buffer : str or netCDF4.Dataset
        Input proxy netCDF file path or opened buffer.
    h5file : str
        Path to output HDF5 file. Not written if None.
    agemodel_iter : int or None, optional
        Optional index of age-model iteration to use in output -- if multiple
        age-model iterations are available.
    icevol_cor : bool, optional
        Do you want to apply ice-volume correction to the d18O foraminiferal
        records? This is done with the `erebusfalls` package. `True` by default.
    find_modern_seasonality : bool, optional
        Do you want to estimate sample seasonality using the modern record?
        Sets annual seasonality if False. Default is True.

    Returns
    -------
    Nothing is returned, output is written to `h5file`.

    Also see
    --------
    `proxysiphon.nc2lmrdf`
    """
    all_data_df, all_meta_df = nc2lmrdf(path_or_buffer, agemodel_iter=agemodel_iter,
                                        icevol_cor=icevol_cor,
                                        find_modern_seasonality=find_modern_seasonality)

    # Write to H5 file.
    log.debug('writing to HDF5 file: {}'.format(h5file))
    all_meta_df.to_hdf(h5file, key='meta', mode='w', format='table',
                       complevel=9, complib='blosc')
    all_data_df.to_hdf(h5file, key='proxy', mode='r+', format='table',
                       complevel=9, complib='blosc')
