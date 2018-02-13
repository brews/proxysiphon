import datetime
import copy
import logging

import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pylab as plt
import cartopy.feature as cfeature
import cartopy.crs as ccrs


log = logging.getLogger(__name__)


def write_reportpdf(qcplot_path, fldata, proxy_vars, latlon, deltar=None, deltar_error=None, proxys_median=None, agemodel=None):
    log.debug('Writing QC report plots')

    n_vars = len(proxy_vars)
    n_cols = n_vars + 1

    if deltar is None and agemodel is not None:
        deltar = agemodel.mcmcsetup.mcmc_kws['d_r'][0]
        deltar_error = agemodel.mcmcsetup.mcmc_kws['d_std'][0]

    with PdfPages(qcplot_path) as pdf:

        fig = plt.figure(figsize=(6.5, 9))

        ax2 = plt.subplot2grid((n_cols, 2), (0, 1),
                               projection=ccrs.Robinson(central_longitude=latlon[1]))
        ax3 = plt.subplot2grid((n_cols, 2), (0, 0))

        ax2 = site_map(latlon, ax=ax2)

        ax3 = summary_text(fldata, deltar_used=deltar, deltar_error_used=deltar_error, ax=ax3)

        for i, var in enumerate(proxy_vars):
            this_ax = plt.subplot2grid((n_cols, 1), (1 + i, 0), colspan=2)
            this_ax = proxyvar_timeseries(proxys_median, var, ax=this_ax)
            this_ax.xaxis.label.set_visible(False)
            this_ax.title.set_visible(False)
            if i == 0:
                this_ax.title.set_visible(True)
                this_ax.set_title('Proxy variables')
        this_ax.xaxis.label.set_visible(True)
        this_ax.set_xlabel('Age (cal yr BP)')
        fig.tight_layout()
        pdf.savefig(bbox_inches='tight')
        plt.close()

        fig = plt.figure(figsize=(6.5, 9))
        ax1 = plt.subplot2grid((3, 2), (0, 0))
        ax2 = plt.subplot2grid((3, 2), (0, 1))
        ax4 = plt.subplot2grid((3, 2), (1, 0), rowspan=2, colspan=2)

        if agemodel is not None:
            ax1 = sedrate(agemodel, ax=ax1)
            ax2 = sedmemory(agemodel, ax=ax2)
            ax4 = agedepth(agemodel, proxys_median, ax=ax4)

        fig.tight_layout()
        pdf.savefig(bbox_inches='tight')
        plt.close()

        fig = plt.figure(figsize=(6.5, 9))
        ax1 = plt.subplot(1, 1, 1)
        if agemodel is not None:
            ax1 = agedepth(agemodel, proxys_median, maxage=25000, ax=ax1)
        fig.tight_layout()
        pdf.savefig(bbox_inches='tight')
        plt.close()

        log.debug('QC report plot saved to {}'.format(qcplot_path))


def proxyvar_timeseries(proxys_median, var, ax=None):
    if ax is None:
        ax = plt.gca()
    # This is handled poorly on my part.
    if 'original_age' in proxys_median.columns:
        # ie...if has new agemodel data
        old_age = proxys_median.loc[:, ('original_age', var)].dropna()
        new_age = proxys_median.loc[:, ('age', var)].dropna()
        ax.plot(new_age.age, new_age[var], '.', color='C3', label='MCMC median age')
        ax.plot(new_age.age, new_age[var], color='C3', linewidth=0.5, label='_nolegend_')
        ax.plot(old_age.original_age, old_age[var], 'x', color='C0', label='File age')
        ax.plot(old_age.original_age, old_age[var], color='C0', linewidth=0.5, label='_nolegend_')
        log.debug('Found new agemodel proxy timeseries')
    else:
        old_age = proxys_median.loc[:, ('age', var)].dropna()
        ax.plot(old_age.age, old_age[var], 'x', color='C0', label='File age')
        ax.plot(old_age.age, old_age[var], color='C0', linewidth=0.5, label='_nolegend_')
        log.debug('Assuming no new agemodel proxy timeseries')

    if 'd18o' in var or 'd18O' in var:
        ax.invert_yaxis()
    ax.set_ylabel(var)
    ax.grid()
    return ax


def agedepth(agemodel, proxy_median, prior_dwidth=30, maxage=None, ax=None):
    """Plot comparing fit agemodel with age in the original proxy data"""
    if ax is None:
        ax = plt.gca()

    # This hack is going to break at some point.
    if maxage is None:
        maxage = agemodel.mcmcsetup.mcmc_kws['maxyr']
    if np.max(agemodel.age_median() > maxage):
        agemodel = copy.deepcopy(agemodel)
        too_old = ~(agemodel.age_median() > maxage)
        agemodel._depth = agemodel.depth[too_old]
        agemodel._age_ensemble = agemodel.age_ensemble[too_old]

    ax = agemodel.plot(ax=ax)
    ax.images[0].set_cmap(plt.cm.Greys)
    for l in ax.lines:
        l.set_color('C3')
    ax.plot(proxy_median.depth, proxy_median.original_age, 'C0', label='File age model')
    ax = agemodel.plot_prior_dates(dwidth=prior_dwidth, ax=ax)
    ax.collections[0].set_color('k')
    ax.collections[0].set_zorder(10)
    ax.autoscale_view()
    ax.set_title('Age model')
    return ax


def site_map(latlon, ax=None):
    """Simple global site map, given latlon"""
    if ax is None:
        ax = plt.gca(projection=ccrs.Robinson(central_longitude=latlon[1]))
    ax.set_global()
    ax.add_feature(cfeature.LAND, facecolor='#B0B0B0')
    ax.outline_patch.set_linewidth(0.5)
    ax.plot(latlon[1], latlon[0], 'o', color='C0', transform=ccrs.Geodetic())
    return ax


def summary_text(fldata, deltar_used=None, deltar_error_used=None, ax=None):
    if ax is None:
        ax = plt.gca()

    latlon = [fldata.site_information.northernmost_latitude,
              fldata.site_information.westernmost_longitude]

    elevation = fldata.site_information.elevation

    # These checks are not well written.
    try:
        deltar_original = np.int(np.round(fldata.chronology_information.df['delta_R'].values[0]))
        if np.isnan(deltar_original):
            deltar_original = None
    except (AttributeError, ValueError, IndexError) as e:
        deltar_original = None

    try:
        deltar_std_original = np.int(np.round(fldata.chronology_information.df['delta_R_1s_err'].values[0]))
        if np.isnan(deltar_std_original):
            deltar_std_original = None
    except (AttributeError, ValueError, IndexError) as e:
        deltar_std_original = None

    try:
        deltar_used = np.int(np.round(deltar_used))
    except (AttributeError, TypeError):
        deltar_used = None

    try:
        deltar_error_used = np.int(np.round(deltar_error_used, 2))
    except (AttributeError, TypeError):
        deltar_error_used = None


    text_template = 'Latitude: {}°\nLongitude: {}°\nElevation: {} m ' \
                    '\n\nΔR: {}\nΔRσ: {}\nFile ΔR: {}\nFile ΔRσ: {}'
    text_str = text_template.format(latlon[0], latlon[1], elevation,
                                    deltar_used, deltar_error_used,
                                    deltar_original, deltar_std_original)
    ax.text(0.05, 0.9, text_str, verticalalignment='top',
            horizontalalignment='left', transform=ax.transAxes)
    ax.set_title(fldata.site_name + '\n' + str(datetime.datetime.now().isoformat()))
    ax.axis('off')
    return ax


def sedrate(agemodel, ax=None):
    """Plot prior and posterior sediment rates for agemodel"""
    if ax is None:
        ax = plt.gca()
    ax = agemodel.plot_sediment_rate(ax)
    ax.lines[-1].set_color('C3')
    return ax


def sedmemory(agemodel, ax=None):
    """Plot prior and posterior sediment memory for agemodel"""
    if ax is None:
        ax = plt.gca()
    ax = agemodel.plot_sediment_memory(ax)
    ax.lines[-1].set_color('C3')
    return ax
