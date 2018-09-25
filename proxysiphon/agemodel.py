import logging

import numpy as np
import pandas as pd
import scipy.stats as stats
import snakebacon as sb
import carbonferret as cf


log = logging.getLogger(__name__)


def remove_outliers(df, col_name, iqr_min=10):
    """Remove outliers from dataframe copy based on col_name"""
    out = df.copy()
    iqr = stats.iqr(df.loc[:, col_name])
    if iqr < iqr_min:
        return out
    upper = np.percentile(df.loc[:, col_name], 75)
    lower = np.percentile(df.loc[:, col_name], 25)
    out = out[out.loc[:, col_name] < (upper + 1.5 * iqr)]
    out = out[out.loc[:, col_name] > (lower - 1.5 * iqr)]
    return out


def get_deltar_online(latlon, max_distance=3000):
    """Use carbonferret to grab an estimate ΔR from internet"""
    max_distance = int(max_distance)
    nearby = cf.find_near(lat=latlon[0], lon=latlon[1], n=10)

    nearby = nearby[nearby['distance (km)'] <= max_distance]
    # nearby = remove_outliers(nearby, 'DeltaR')
    # nearby = remove_outliers(nearby, 'DeltaRErr')
    log.debug('ΔR and ΔRσ from {} samples'.format(len(nearby)))

    w = 1/len(nearby)
    deltar_mean = nearby['DeltaR'].mean()
    # Use pooled or combined variance of carbonferret deltaR distributions.
    var = (w * ((nearby['DeltaR'] - deltar_mean)**2 + nearby['DeltaRErr']**2)).sum()
    sigma = np.sqrt(var)
    return tuple([deltar_mean, sigma])


def fit_agedepthmodel(chron, pdata, deltar=None, deltar_error=None, minyr=None, mcmc_kws=None):
    log.debug('Fitting new age model.')
    if minyr is None:
        minyr = -1000

    chron = pd.DataFrame(chron).copy()
    pdata = pd.DataFrame(pdata).copy()

    chron['depth'] = chron.loc[:, ('depth_top', 'depth_bottom')].mean(axis=1)
    chron['labid'] = chron.Labcode
    chron['age'] = chron.loc[:, '14C_date']
    chron['error'] = chron.loc[:, '14C_1s_err']

    chron['cc'] = 2  # This is marine calibration curve, use for all depths.

    chron.sort_values('depth', inplace=True)

    # Cleanup dataframe so works with snakebacon.
    if deltar is not None:
        chron['delta_R'] = deltar
    if deltar_error is not None:
        chron['delta_R_1s_err'] = deltar_error

    if chron.other_date.notnull().any():
        # Check that can't have 14C_date and other_date at same depth. Same for 1s_error.
        assert ~chron.loc[:, ['14C_date', 'other_date']].notnull().all(axis=1).any()
        assert ~chron.loc[:, ['14C_1s_err', 'other_1s_err']].notnull().all(axis=1).any()
        other_msk = ~chron['14C_date'].notnull()
        # `other_dates` should not have delta_R values
        chron.loc[other_msk, 'delta_R'] = 0
        chron.loc[other_msk, 'delta_R_1s_err'] = 0
        # Move `other_dates` and `errors` to chron.age and chron.error.
        chron.loc[other_msk, 'age'] = chron.loc[other_msk, 'other_date']
        chron.loc[other_msk, 'error'] = chron.loc[other_msk, 'other_1s_err']
        chron.loc[other_msk, 'cc'] = 0  # Using constant calibration curve for non-14C dates
        # Drop rows with `other_date` but no `other_1s_error`.
        chron.drop(chron[other_msk & chron['other_date'].notnull() & ~chron['other_1s_err'].notnull()].index, inplace=True)

    # Have any NaNs in age, depth or error?
    assert chron.loc[:, ['age', 'depth', 'error']].notnull().all().all()

    coredates = sb.ChronRecord(chron)

    d_min = np.min([x.min() for x in [chron.depth, pdata.depth]])
    d_max = np.max([x.max() for x in [chron.depth, pdata.depth]])
    sug_acc_mean = sb.suggest_accumulation_rate(coredates)
    # n_segs = np.ceil((d_max - d_min) / 5)  # Num segments in mcmc, ~ 5cm, rounded up.

    guesses = np.random.randn(2) * coredates.error[:2] + coredates.age[:2]  # TODO(brews): Check whether need sqrt error here.
    guesses[guesses < minyr] = minyr  # Line #70 of Bacon.R warns that otherwise twalk MCMC will not run.

    # if n_segs > 200 or n_segs < 5:
    # n_segs = np.ceil((d_max - d_min) / 10)
    # assert (n_segs < 500) and (n_segs > 5)  # Sanity check for extremely long or short MCMC runs.

    mcmc_params = dict(depth_min=d_min, depth_max=d_max,
                       cc=chron.cc.values,
                       d_r=chron.delta_R.values,
                       d_std=chron.delta_R_1s_err.values,
                       t_a=[3], t_b=[4], k=50,#n_segs,
                       minyr=minyr, maxyr=50000,
                       th01=guesses[0], th02=guesses[1],
                       acc_mean=sug_acc_mean, acc_shape=1.5,
                       mem_strength=4, mem_mean=0.7)

    if mcmc_kws is not None:
        mcmc_params.update(mcmc_kws)

    log.debug('MCMC parameters: {}'.format(mcmc_params))

    agemodel = sb.AgeDepthModel(coredates, mcmc_kws=mcmc_params)
    log.debug('Age model fit done')
    return agemodel, coredates, mcmc_params



def date_proxy(admodel, pdata, nsims):
    pdata = pd.DataFrame(pdata).copy()
    if 'age' in pdata.columns:
        # Rename to avoid error with sb.ProxyRecord
        log.debug('Renaming age column in proxy data')
        cols = list(pdata.columns.values)
        cols[cols.index('age')] = 'original_age'
        pdata.columns = cols

    orig_pdata = sb.ProxyRecord(pdata)
    pdata_median = admodel.date(orig_pdata, 'median').to_pandas()
    pdata_ensemble = admodel.date(orig_pdata, 'ensemble', nsims).to_pandas()
    pdata_median['age'] = pdata_median['age'].round()
    pdata_ensemble['age'] = pdata_ensemble['age'].round()
    return pdata_median, pdata_ensemble
