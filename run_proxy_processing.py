#! /usr/bin/env python

# Read in raw proxy files and process - cleaning and fitting age models.
#
# Can run from Bash shell with:
#
# python ./run_proxyprocessing.py \
#     --proxydir data/ncdc_proxies_parsed \
#     --daproxyfile data/data_assimilation/proxies/DADT_v0.0.0_Proxies.h5 \
#     --qcplotdir data/proxyagemodel/plots_proxyqc \
#     --agemodeldir data/proxyagemodel/pickledagemodels \
#     --simdir data/proxyagemodel/mcmcagemodel_proxies \
#     --agemodelconfigdir data/proxyagemodel/config \
#     --nsims 1000
#
# Otherwise you're going to want to import and run through `main()`.

import logging
import os
import pickle
import argparse
from glob import glob

import matplotlib
matplotlib.use('Agg')  # Needed to make QC plots on server.

import yaml
import numpy as np
import pandas as pd

from proxysiphon import read_ncdc, get_deltar_online, fit_agedepthmodel, date_proxy, qcreport


log = logging.getLogger(__name__)


def proxymedian_dfs(obj=None, proxy_medians=None):
    """Return proxy median data and metadata pandas df needed for DA input

    If no args are passed, return empty (template) pandas dataframes.
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
    all_meta = meta_template.copy()
    if obj is None and proxy_medians is None:
        return data_template.copy(), meta_template.copy()

    vars = list_proxyvars(obj)
    latlon = [obj.site_information.northernmost_latitude,
              obj.site_information.westernmost_longitude]
    elevation = obj.site_information.elevation
    data_df = data_template.copy()
    meta_df = meta_template.copy()

    for v in vars:
        d_subset = proxy_medians[['age', v]].dropna()
        youngest_ce = np.max(1950 - d_subset['age'])
        oldest_ce = np.min(1950 - d_subset['age'])
        proxyid = obj.site_name + ':' + str(v)
        this_meta = meta_template.copy()
        this_meta['Proxy ID'] = [proxyid]
        this_meta['Site'] = [obj.site_name]
        this_meta['Lat (N)'] = latlon[0]
        this_meta['Lon (E)'] = latlon[1]
        this_meta['Archive type'] = ['Marine sediments']
        this_meta['Proxy measurement'] = [v]
        this_meta['Resolution (yr)'] = [(youngest_ce - oldest_ce) / len(d_subset)]
        this_meta['Reference'] = [str(None)]
        this_meta['Databases'] = ['[DTDA]']
        this_meta['Seasonality'] = ['[1,2,3,4,5,6,7,8,9,10,11,12]']
        this_meta['Elev'] = elevation
        this_meta['Oldest (C.E.)'] = oldest_ce
        this_meta['Youngest (C.E.)'] = youngest_ce
        meta_df = meta_df.append(this_meta, ignore_index=True)

        d_subset = d_subset.rename(columns={v: proxyid, 'age': 'Year C.E.'})
        d_subset['Year C.E.'] = 1950 - d_subset['Year C.E.']
        d_subset = d_subset.set_index('Year C.E.')
        data_df = data_df.join(d_subset, how='outer')

    data_df = data_df.sort_index(ascending=False)
    return data_df, meta_df


def list_proxyvars(obj):
    """Return list of proxy variable names for dataframe or ncdcrecord"""
    try:
        proxy_vars = obj.columns
    except AttributeError:
        proxy_vars = obj.data.df.columns
    proxy_vars = list(proxy_vars)

    try:
        proxy_vars.remove('depth')
    except ValueError:
        pass

    try:
        proxy_vars.remove('age')
    except ValueError:
        pass

    try:
        proxy_vars.remove('original_age')
    except ValueError:
        pass

    return proxy_vars


def has_data(obj):
    """Check if has populated data information"""
    try:
        d = obj.data.df.copy()
    except KeyError:
        return False
    if len(d) > 0:
        result = True
    else:
        result = False
    return result


def has_chron(obj):
    """Check if has populated chronology information"""
    try:
        chron = obj.chronology_information.df.copy()
    except KeyError:
        return False
    if len(chron) > 1: # Needs to have more than one date.
        result = True
    else:
        result = False
    return result


def get_recent_date(obj):
    """Get the most recent date from record object.

    First try "Collection_Year" in "Data_Collection" section. If can't find, then
    "Date" in "Contribution_Date" section. Otherwise returns None.
    """
    col_year = obj.data_collection.collection_year
    if col_year is not None:
        return col_year
    # pub_year = obj.publication
    # if pub_year is not None:
    #     return pub_year
    con_date = obj.contribution_date.date
    try:
        return con_date.year
    except AttributeError:
        return None


def main(proxydir, daproxyfile, qcplotdir, agemodeldir,
         simdir, nsims, agemodelconfigdir=None):
    os.makedirs(agemodeldir, exist_ok=True)
    os.makedirs(qcplotdir, exist_ok=True)
    os.makedirs(simdir, exist_ok=True)

    glob_str = '*.txt'
    log.info('New run')
    targets = glob(os.path.join(proxydir, glob_str))
    targets.sort()

    proxysimmetadata_all = pd.DataFrame(columns=['Site', 'Lat (N)', 'Lon (E)'])
    alldata_df, allmeta_df = proxymedian_dfs()

    for fl in targets:
        log.info('Starting proxy file: {}'.format(fl))

        try:
            fldata = read_ncdc(fl)
            fldata.site_name = ''.join(os.path.basename(fl).split('.')[:-1])  # Filename without extension.
        except Exception as e:
            log.exception('Exception while reading proxy file:\n{}'.format(e))
            continue

        site_config = None
        try:
            if agemodelconfigdir is not None:
                config_path = os.path.join(agemodelconfigdir, fldata.site_name + '.yaml')
                try:
                    log.debug('Trying YAML config in {}'.format(config_path))
                    with open(config_path, 'r') as f:
                        site_config = yaml.safe_load(f)
                except FileNotFoundError:
                    log.debug('Found no config in {}'.format(config_path))
        except Exception as e:
            log.exception('Exception while reading configuration file:\n{}'.format(e))
            continue

        if not has_data(fldata):
            log.warning('Ditching file {} - no data information'.format(fl))
            continue

        latlon = [fldata.site_information.northernmost_latitude,
                  fldata.site_information.westernmost_longitude]
        if None in latlon:
            log.warning('No latlon data in {} - skipping file'.format(fl))
            continue
        log.debug('Found latlon: ' + str(latlon))

        data = fldata.data.df.copy()
        chron = fldata.chronology_information.df.copy()

        agemodel = None
        deltar = None
        deltar_error = None
        proxys_median = data.copy()

        try:
            if has_chron(fldata) and len(data['depth'].dropna().unique()) > 1:
                try:
                    if data.depth.isnull().all():
                        log.warning('Ditching file {} - data has no values in datas `depth` field'.format(fl))
                        continue
                except AttributeError:
                    log.warning('Ditching file {} - data has no `depth` field in data'.format(fl))
                    continue

                if 'depth' not in data.columns:
                    log.warning('Ditching file {} - data has no `depth` field'.format(fl))
                    continue

                n_unique_deltarerror = len(chron['delta_R_1s_err'].dropna().unique())
                n_unique_deltar = len(chron['delta_R'].dropna().unique())
                if (n_unique_deltarerror > 1) and (n_unique_deltar > 1):
                    log.info('Found multi-depth deltar and deltar_error. Using file values')
                elif n_unique_deltar > 1:
                    log.info('Found multi-depth deltar, without deltar_error. Using carbonferret deltar_error.')
                    deltar_error = get_deltar_online(latlon)[1]
                else:
                    deltar, deltar_error = get_deltar_online(latlon)
                    log.debug('deltar(deltar_error): {}({})'.format(deltar, deltar_error))

                pickle_path = os.path.join(agemodeldir, fldata.site_name + '.pickle')
                try:
                    log.debug('Trying pickled AgeDepthModel in {}'.format(pickle_path))
                    with open(pickle_path, 'rb') as f:
                        agemodel = pickle.load(f)
                    log.info('Pickled AgeDepthModel loaded')
                except FileNotFoundError:
                    log.debug('Did not find pickled AgeDepthModel')
                    try:
                        mcmc_kws = site_config['mcmc_kws']
                    except TypeError:
                        mcmc_kws = None
                    try:
                        agemodel, chrondata = fit_agedepthmodel(chron, data,
                                                                deltar, deltar_error,
                                                                minyr=1950 - get_recent_date(fldata),
                                                                mcmc_kws=mcmc_kws)
                        with open(pickle_path, 'wb') as f:
                            pickle.dump(agemodel, f, pickle.HIGHEST_PROTOCOL)
                        log.debug('Pickled age model dumped to {}'.format(pickle_path))
                    except Exception as e:
                        agemodel = None
                        log.exception('Error in agemodel fitting. Continuing without agemodel.\n{}'.format(e))
                if agemodel is not None:
                    log.debug('Dating proxies with age model')
                    proxys_median, proxys_ensemble = date_proxy(agemodel, data, nsims)

                    # Output metadata and agemodel simulations to hard file.
                    try:
                        log.debug('Compiling proxy simulations')
                        proxys_ensemble['Site'] = fldata.site_name
                        if 'original_age' in proxys_ensemble.columns:
                            proxys_ensemble.drop('original_age', axis=1, inplace=True)
                        proxys_ensemble = pd.melt(proxys_ensemble, id_vars=['Site', 'depth', 'mciter', 'age'])
                        proxys_ensemble = proxys_ensemble[proxys_ensemble.notnull()]
                        out_path = os.path.join(simdir, 'mcmc-' + fldata.site_name + '.h5')
                        proxys_ensemble.to_hdf(out_path, key='sim', mode='w', complevel=9, complib='blosc')

                        log.debug('Proxy simulations written to {}'.format(out_path))

                        log.debug('Appending metadata')
                        proxymetadata_this = pd.DataFrame({'Site': [fldata.site_name],
                                                           'Lat (N)': np.array([latlon[0]], dtype='float64'),
                                                           'Lon (E)': np.array([latlon[1]], dtype='float64')})
                        proxysimmetadata_all = proxysimmetadata_all.append(proxymetadata_this, ignore_index=True)
                        log.debug('Metadata appended')

                    except Exception as e:
                        log.exception('Exception while compiling proxysims and metadata:\n{}'.format(e))

            else:
                log.info('File {} - no chronology information or no depth in data'.format(fl))

        except Exception as e:
            log.exception('Exception in agemodel fitting:\n{}'.format(e))

        try:
            file_datadf, file_metadf = proxymedian_dfs(fldata, proxys_median)
            allmeta_df = allmeta_df.append(file_metadf, ignore_index=True)
            alldata_df = alldata_df.join(file_datadf, how='outer')

        except Exception as e:
            log.exception('Exception in proxy median metadata:\n{}'.format(e))

        try:
            qcplot_path = os.path.join(qcplotdir, fldata.site_name + '.pdf')
            proxy_vars = list_proxyvars(data)
            qcreport.plot.write_reportpdf(qcplot_path, fldata, proxy_vars, latlon,
                                          proxys_median=proxys_median,
                                          agemodel=agemodel)
            log.debug('Done with QC report plot')
        except Exception as e:
            log.exception('Exception in QC report:\n{}'.format(e))
            continue

        log.debug('Done with file')

    metadata_path = os.path.join(simdir, 'metadata.h5')
    proxysimmetadata_all.to_hdf(metadata_path, key='meta', mode='w', complevel=9, complib='blosc')
    log.debug('Proxy sim metadata written to {}'.format(metadata_path))

    allmeta_df.to_hdf(daproxyfile, key='meta', mode='w', format='table',
                      complevel=9, complib='blosc')
    log.debug('Proxy median metadata written to {}'.format(daproxyfile))
    alldata_df = alldata_df.sort_index(ascending=False)
    alldata_df.to_hdf(daproxyfile, key='proxy', mode='r+', format='table',
                      complevel=9, complib='blosc')
    log.debug('Proxy median data written to {}'.format(daproxyfile))

    log.info('Run done')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='NCDC proxy file processing')
    parser.add_argument('--proxydir', metavar='PROXYPATH', nargs=1, 
                        help='directory to search for NOAA NCDC proxy .txt files to be parsed')
    parser.add_argument('--daproxyfile', metavar='DAPROXYPATH', nargs=1,
                        help='Path that DA proxy file dataframe HDF5 file will be written to')
    parser.add_argument('--qcplotdir', metavar='QCPLOTPATH', nargs=1,
                        help='directory QC plots will be written to')
    parser.add_argument('--agemodeldir', metavar='AGEMODELPATH', nargs=1, 
                        help='directory that pickled AgeDepthModels will be written to')
    parser.add_argument('--agemodelconfigdir', metavar='AGEMODELCONFIGPATH',
                        nargs=1, default=None,
                        help='directory to check for YAML config files')
    parser.add_argument('--simdir', metavar='SIMPATH', nargs=1, 
                        help='directory that mcmc agemodel proxy simulations will be written to')
    parser.add_argument('--nsims', type=int, metavar='N', nargs=1,
                        default=1000,
                        help='number of random mcmc agemodel proxy simulations to save')
    parser.add_argument('--log', metavar='LOGPATH', nargs=1,
                        default=['proxyprocessing.log'],
                        help='file path that log will be written to')
    args = parser.parse_args()

    logging.basicConfig(level=logging.DEBUG,
                        filename=args.log[0],
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    logging.getLogger('urllib3.connectionpool').setLevel(logging.INFO)
    logging.getLogger('chardet.charsetprober').setLevel(logging.INFO)

    main(proxydir=args.proxydir[0], daproxyfile=args.daproxyfile[0],
         qcplotdir=args.qcplotdir[0], agemodeldir=args.agemodeldir[0],
         agemodelconfigdir=args.agemodelconfigdir[0],
         simdir=args.simdir[0], nsims=args.nsims[0])
