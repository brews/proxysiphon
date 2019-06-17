import logging
import datetime
from copy import deepcopy
import unidecode
import numpy as np
import netCDF4
from proxysiphon.agemodel import get_deltar_online, fit_agedepthmodel, date_proxy


log = logging.getLogger(__name__)


def _normalize_to_ascii_array(a, dtype='S50'):
    """Normalize sequence of UTF-8 string to np.Array of ASCII"""
    normed = np.array([unidecode.unidecode(str(x)) for x in a], dtype=dtype)
    return normed


class RedateMixin:
    """Mixins to redate LGM proxy records"""

    def swapin_custom_deltar(self, d_r=None, d_std=None):
        """Swap-in custom proxy ΔR mean and standard-deviation values.

        ``d_r`` is put into self.chronology_information.df column "delta_R" and
        ``d_std`` is put into column "delta_R_1s_err". Any existing values in
        these columns are first moved to a new column "*_original", but only if
        the "*_original" columns don't already exist. If "*_original" columns
        already exist, then current values in "delta_R" or "delta_R_1s_err" are
        discarded in the returned proxy record.

        Parameters
        ----------
        d_r : ndarray, scalar or None
            Carbon reservoir (ΔR) mean value to swap-in. Ignored if None.
        d_std : ndarray, scalar or None
            ΔR standard-deviation value to swap-in. Ignored
            if None.

        Returns
        -------
        out : Modified copy of proxy record.
        """
        out = self.copy()

        if d_r is None and d_std is None:
            return out

        # If we're plugging in a new value, we preserve the original but moving it to a new *_original
        # columns.
        if d_r is not None:
            if 'delta_R_original' not in out.chronology_information.df.columns:
                out.chronology_information.df['delta_R_original'] = out.chronology_information.df['delta_R']
            out.chronology_information.df['delta_R'] = d_r
            out.chronology_information.df['delta_R'].replace(to_replace='None', value=np.nan, inplace=True)

        if d_std is not None:
            if 'delta_R_1s_err_original' not in out.chronology_information.df.columns:
                out.chronology_information.df['delta_R_1s_err_original'] = out.chronology_information.df[
                    'delta_R_1s_err']
            out.chronology_information.df['delta_R_1s_err'] = d_std
            out.chronology_information.df['delta_R_1s_err'].replace(to_replace='None', value=np.nan, inplace=True)

        return out

    def copy(self):
        """Return deep copy of self"""
        return deepcopy(self)

    def redate(self, **kwargs):
        """Redate proxy data with (snake)bacon

        Parameters
        ----------
        nsims : int
            Number of age-model simulations to retain.
        **kwargs :
            Passed on to ``self._fit_agemodel``.

        Returns
        -------
        Updated copy of self. The ``snakebacon.AgeDepthModel`` used to redate
        the record is monkeypatched into
        ``self.chronology_information.bacon_agemodel`` in this returned copy.
        """
        rec = self.copy()

        agemodel, mcmc_kwargs = rec._fit_agemodel(**kwargs)
        a_median, a_ensemble = rec._date_sampledepths(agemodel, kwargs.pop('nsims', None))
        rec.chronology_information.bacon_agemodel = agemodel
        rec.data.age_median = a_median
        rec.data.age_ensemble = a_ensemble
        return rec

    def _fit_agemodel(self, **kwargs):
        """Fit snakebacon model to NcdcRecord
        """
        chron_df = self.chronology_information.df.copy()
        data_df = self.data.df.copy()
        myr = 1950 - self.recent_date()
        deltar = self.chronology_information.df['delta_R'].values
        deltar_error = self.chronology_information.df['delta_R_1s_err'].values

        agemodel, _, mcmc_kwargs = fit_agedepthmodel(chron_df, data_df, deltar, deltar_error,
                                                     minyr=myr, mcmc_kws=kwargs)
        return agemodel, mcmc_kwargs

    def _date_sampledepths(self, agemodel, nsims=1000):
        """Date self proxy sample depths given a snakebacon.AgeModel

        Parameters
        ----------
        agemodel : snakebacon.AgeModel
            Age model to date from.
        nsims : scalar or None, optional
            Number of draws to include in returned ensemble. ``None`` returns
            1000 draws.

        Returns
        -------
        p_median : pandas.DataFrame
            Median age per depth.
        p_ensemble :
            Ensemble of ages per depth.
        """
        if nsims is None:
            nsims = 1000

        data_df = self.data.df.copy()
        p_median, p_ensemble = date_proxy(agemodel, data_df, nsims)

        p_median = (p_median[['depth', 'age']].rename(columns={'age': 'age_median'})
                    .set_index('depth'))
        p_ensemble = (p_ensemble[['depth', 'age', 'mciter']].rename(columns={'mciter': 'draw'})
                      .pivot(index='depth', columns='draw', values='age'))
        return p_median, p_ensemble

    def recent_date(self):
        """Get the most recent date from self metadata.

        First try "Collection_Year" in "Data_Collection" section. If can't
        find, then earliest year in publications list.
        """
        out = None
        col_year = self.data_collection.collection_year
        if col_year is not None:
            out = col_year
        else:
            # Return earliest year in publications.
            out = min([int(p.published_date_or_year) for p in self.publication if p.published_date_or_year is not None])
        return out

    def update_deltar(self):
        """Set self.ChronologyInformation.df delta_R and delta_r_1s_error, returning an updated copy

        If these variables are changed/updated, the originals are moved to column names *_original.
        """
        x = self.copy()
        latlon = (float(x.site_information.northernmost_latitude),
                  float(x.site_information.easternmost_longitude))

        chron_df = x.chronology_information.df.copy()

        delta_r_used = None
        delta_r_1s_err_used = None

        n_unique_deltarerror = len(chron_df['delta_R'].dropna().unique())
        n_unique_deltar = len(chron_df['delta_R_1s_err'].dropna().unique())

        if (n_unique_deltarerror > 1) and (n_unique_deltar > 1):
            log.info('Found multi-depth deltar and deltar_error. Using NcdcRecords original values')
        elif n_unique_deltar > 1:
            log.info('Found multi-depth deltar, without deltar_error. Using carbonferret "delta_R_1s_err".')
            delta_r_1s_err_used = get_deltar_online(latlon)[1]
        else:
            delta_r_used, delta_r_1s_err_used = get_deltar_online(latlon)
            log.debug('deltar(deltar_error): {}({})'.format(delta_r_used, delta_r_1s_err_used))

        if delta_r_used is not None:
            chron_df['delta_R_original'] = chron_df['delta_R'].copy()
            chron_df['delta_R'] = delta_r_used

        if delta_r_1s_err_used is not None:
            chron_df['delta_R_1s_err_original'] = chron_df['delta_R_1s_err'].copy()
            chron_df['delta_R_1s_err'] = delta_r_1s_err_used

        x.chronology_information.df = chron_df.copy()
        return x

    def average_duplicate_datadepths(self):
        """Average obs at the same depth in self.data.df, return copy of self
        """
        if not self.data.df['depth'].duplicated().any():
            log.debug('Found no duplicate depths in Data')
            # Return if no duplicates
            return self

        x = self.copy()
        log.debug('Found duplicate depths in Data')
        x.data.df = x.data.df.groupby('depth', as_index=False).mean()
        return x

    def has_chron(self):
        """True if self has populated chronology information"""
        result = False
        try:
            chron = self.chronology_information.df
        except KeyError:
            return result
        if len(chron) > 1:  # Needs to have more than one date.
            result = True
        return result

    def has_data(self):
        """True if self has populated data information"""
        result = False
        try:
            d = self.data.df
        except KeyError:
            return result
        if len(d) > 0:
            result = True
        return result

    def slice_datadepths(self, shallow=None, deep=None):
        """Cut self.data.df to given depths, return updated copy of self

        The cut is inclusive for `shallow` and `deep` limits. Both `shallow` and
        `deep` must be positive. If no values are given for `shallow` and `deep`,
        the cuts are at the deepest and shallowest samples in
        `self.chronology_information`.
        """
        if len(self.data.df) < 1 or 'depth' not in self.data.df.columns:
            return self

        out = self.copy()

        if shallow is None:
            shallow = out.data.df['depth'].min()
        if deep is None:
            deep = out.data.df['depth'].max()

        assert shallow >= 0 and deep >= 0, 'cut depths must be positive'

        out.data.df = out.data.df[(out.data.df['depth'] >= shallow) & (out.data.df['depth'] <= deep)]

        return out


class NetcdfMixin:
    """Mixins to export LGM proxy records to a NetCDF file"""

    @staticmethod
    def _variable_attributes(varname):
        """Get dict of netCDF4 variable attributes for a given NcdcRecord Data column name"""
        proxy_map = {'d13c': ('d13C', 'per mil'),
                     'd18o': ('d18O', 'per mil'),
                     'mgca': ('Mg/Ca', 'per mil'),
                     'percent': ('Percent foraminifera', '%'),
                     'tex86': ('TEX86', 'index'),
                     'uk37': ("UK'37", 'index')}
        # Second part of the column name...
        foraminifera_map = {'bulloides': 'Globigerina bulloides',
                            'crassaformis': 'Globorotalia crassaformis',
                            'dutertrei': 'Neogloboquadrina dutertrei',
                            'inflata': 'Globoconella inflata',
                            'mabahethi': 'Cibicides mabahethi',
                            'marginata': 'Bulimina marginata',
                            'menardii': 'Globorotalia menardii',
                            'obliquiloculata': 'Pulleniatina obliquiloculata',
                            'pachyderma': 'Neogloboquadrina pachyderma sinistral',
                            'pachysin': 'Neogloboquadrina pachyderma sinistral',
                            'pachyderma_d': 'Neogloboquadrina incompta',
                            'peregrina': 'Uvigerina peregrina',
                            'quinqueloba': 'Turborotalita quinqueloba',
                            'ruber': 'Globigerinoides ruber white',
                            'ruber_lato': 'Globigerinoides ruber white',
                            'ruber_pink': 'Globigerinoides ruber pink',
                            'ruber_stricto': 'Globigerinoides ruber white',
                            'sacculifer': 'Trilobatus sacculifer',
                            'truncatulinoides': 'Globorotalia pachytheca',
                            'tumida': 'Globorotalia tumida',
                            'acicula': 'Creseis acicula',
                            }
        varname = varname.lower()
        try:
            if '_' not in varname:
                out = {'long_name': proxy_map[varname][0], 'units': proxy_map[varname][1]}
            else:
                proxy, foram = varname.split('_', 1)
                out = {'long_name': proxy_map[proxy][0], 'units': proxy_map[proxy][1],
                       'foraminifera_type': foraminifera_map[foram]}
        except KeyError:  # Variable name not found.
            out = {}
        return out

    def _attach_site_ncgroup(self, parent):
        """Create the NetCDF4 site group and populate it

        This is run whenever self.to_netcdf() is called.

        Parameters
        ----------
        parent : netcdf4.Group or netcdf4.Dataset
            Object the created site group will use as a parent.

        Returns
        -------
        this_site : netcdf4.Group
            Reference to the created chronology group.
        """
        site_name = self.site_information.site_name.strip()
        # Normalize site group name to ASCII and replace whitespace characters.
        grp_name = unidecode.unidecode(site_name.lower()).replace(' ', '_')

        # Create and populate site group
        this_site = parent.createGroup(grp_name)
        this_site.site_name = str(site_name)
        this_site.comment = str(self.description)
        this_site.latitude = float(self.site_information.northernmost_latitude)
        this_site.longitude = float(self.site_information.easternmost_longitude)
        this_site.elevation = int(self.site_information.elevation)
        try:
            this_site.collection_year = int(self.data_collection.collection_year)
        except TypeError:  # If collection year doesn't exist
            pass
        publication_str = '\n\n'.join([x.to_citationstr() for x in self.publication])
        this_site.references = publication_str
        return this_site

    def _attach_chronology_ncgroup(self, parent):
        chron = parent.createGroup('chronology')
        chron.createDimension('depth_top', None)
        chron.createDimension('str_dim', 50)

        labcode = chron.createVariable('labcode', 'S1', ('depth_top', 'str_dim'))
        labcode.long_name = 'Lab sample code'
        labcode.missing_value = 'NA'
        labcode._Encoding = 'ascii'
        labcode[:] = _normalize_to_ascii_array(self.chronology_information.df['Labcode'].fillna('NA'))

        depth_top = chron.createVariable('depth_top', 'f4', ('depth_top',))
        depth_top.long_name = 'Sample top depth'
        depth_top.axis = 'Z'
        depth_top.positive = 'down'
        depth_top.units = 'cm'
        depth_top.missing_value = np.nan
        depth_top[:] = self.chronology_information.df['depth_top'].values

        depth_bottom = chron.createVariable('depth_bottom', 'f4', ('depth_top',))
        depth_bottom.long_name = 'Sample bottom depth'
        depth_bottom.units = 'cm'
        depth_bottom.positive = 'down'
        depth_bottom.missing_value = np.nan
        depth_bottom[:] = self.chronology_information.df['depth_bottom'].values

        mat_dated = chron.createVariable('mat_dated', 'S1', ('depth_top', 'str_dim'))
        mat_dated.long_name = 'Material dated'
        mat_dated.missing_value = 'NA'
        mat_dated._Encoding = 'ascii'
        mat_dated[:] = _normalize_to_ascii_array(self.chronology_information.df['mat_dated'].fillna('NA'))

        c14_date = chron.createVariable('c14_date', 'f4', ('depth_top',))
        c14_date.long_name = '14C date'
        c14_date.units = 'RC yr BP'
        c14_date.missing_value = np.nan
        c14_date[:] = self.chronology_information.df['14C_date'].values

        c14_1s_err = chron.createVariable('c14_1s_err', 'f4', ('depth_top',))
        c14_1s_err.long_name = '14C 1-sigma error'
        c14_date.units = 'RC yr BP'
        c14_1s_err.missing_value = np.nan
        c14_date[:] = self.chronology_information.df['14C_1s_err'].values

        delta_r = chron.createVariable('delta_r', 'f4', ('depth_top',))
        delta_r.long_name = 'delta R'
        delta_r.missing_value = np.nan
        delta_r[:] = self.chronology_information.df['delta_R'].values

        if 'delta_R_original' in self.chronology_information.df.columns:
            delta_r_orig = chron.createVariable('delta_r_original', 'f4', ('depth_top',))
            delta_r_orig.long_name = 'Original delta R'
            delta_r_orig.description = 'Carbon reservoir correction (delta R) value(s) given in the orignal proxy site data set.'
            delta_r_orig.missing_value = np.nan
            delta_r_orig[:] = self.chronology_information.df['delta_R_original'].values

        delta_r_1s_error = chron.createVariable('delta_r_1s_error', 'f4', ('depth_top',))
        delta_r_1s_error.missing_value = np.nan
        delta_r_1s_error.long_name = 'delta R 1-sigma error'
        delta_r_1s_error[:] = self.chronology_information.df['delta_R_1s_err'].values

        if 'delta_R_1s_err_original' in self.chronology_information.df.columns:
            delta_r_1s_error_orig = chron.createVariable('delta_r_1s_error_original', 'f4', ('depth_top',))
            delta_r_1s_error_orig.long_name = 'Original delta R 1-sigma error'
            delta_r_1s_error_orig.description = 'Carbon reservoir correction 1-sigma error value(s) given in the orignal proxy site data set.'
            delta_r_1s_error_orig.missing_value = np.nan
            delta_r_1s_error_orig[:] = self.chronology_information.df['delta_R_1s_err_original'].values

        other_date = chron.createVariable('other_date', 'f4', ('depth_top',))
        other_date.missing_value = np.nan
        other_date.long_name = 'Other date'
        other_date[:] = self.chronology_information.df['other_date'].values

        other_1s_err = chron.createVariable('other_1s_err', 'f4', ('depth_top',))
        other_1s_err.long_name = 'Other date 1-sigma error'
        other_1s_err.missing_value = np.nan
        other_1s_err[:] = self.chronology_information.df['other_1s_err'].values

        other_type = chron.createVariable('other_type', 'S1', ('depth_top', 'str_dim'))
        other_type.long_name = 'Other date type'
        other_type.missing_value = 'NA'
        other_type._Encoding = 'ascii'
        other_type[:] = _normalize_to_ascii_array(self.chronology_information.df['other_type'].fillna('NA'))
        return chron

    def _attach_data_ncgroup(self, parent):
        """Create and populate data group"""
        data = parent.createGroup('data')
        data.createDimension('depth', None)
        depth = data.createVariable('depth', 'f4', ('depth',), zlib=True)
        depth.long_name = 'Sample depth'
        depth.positive = 'down'
        depth.axis = 'Z'
        depth[:] = self.data.df['depth'].values

        file_depth_unit = str(self.variables['depth'].units)
        if file_depth_unit == '' or file_depth_unit is None:
            depth.units = 'cm'
        else:
            depth.units = file_depth_unit


        age_original = data.createVariable('age_original', 'f4', ('depth',),
                                           zlib=True)
        age_original.missing_value = np.nan
        age_original.long_name = 'Original age'
        age_original[:] = self.data.df['age'].values

        file_age_unit = str(self.variables['age'].units)
        if file_age_unit == '' or file_age_unit is None:
            age_original.units = 'cal years BP'
        else:
            age_original.units = file_age_unit

        if hasattr(self.data, 'age_ensemble') and hasattr(self.data, 'age_median'):
            data.createDimension('draw', self.data.age_ensemble.shape[1])

            age_median = data.createVariable('age_median', 'f4', ('depth',),
                                             zlib=True)
            age_median.units = 'cal years BP'
            age_median.long_name = 'Median age'
            age_median.missing_value = np.nan
            age_median[:] = self.data.age_median['age_median'].values

            agedraw = data.createVariable('age_ensemble', 'f4', ('depth', 'draw'),
                                          zlib=True)
            agedraw.units = 'cal years BP'
            agedraw.long_name = 'Age ensemble'
            agedraw.missing_value = np.nan
            agedraw[:] = self.data.age_ensemble.values

        for col in list(self.data.df.columns):
            col_name = col.lower()

            if col_name in ['depth', 'age']:
                continue

            var = data.createVariable(col_name, 'f4', ('depth',), zlib=True)
            var.missing_value = np.nan

            # Add more attributes to variable.
            attrib_dict = self._variable_attributes(col_name)
            for k, v in attrib_dict.items():
                setattr(var, k, v)

            # Overwrite units attributes with whatever units are given in the NcdcRecord.
            var.units = str(self.variables[col].units)
            var.comments = str(self.variables[col].detail)

            # Grab Mg/Ca cleaning information from "Data Collection Information - Notes"
            if 'mgca' in col_name:
                cleaning_note = str(self.data_collection.notes)
                if 'mg_bcp' in cleaning_note:
                    var.mgca_cleaning_protocol = 'Barker cleaning with hydrogen peroxide'
                elif 'mg_red' in cleaning_note:
                    var.mgca_cleaning_protocol = 'Fully reductive cleaning'
                else:
                    var.mgca_cleaning_protocol = 'NA'
            var[:] = self.data.df[col].values
        return data

    def _attach_ncgroups(self, fl):
        """Dump contents into netCDF4.Dataset

        This is run whenever self.to_netcdf() is called.

        This runs all of the self._attach_*ncgroup() methods and attaches them
        to a netcdf4.Dataset object.
        """
        site_group = self._attach_site_ncgroup(fl)
        # Create and populate chronology group, if chronology_information exists
        if not self.chronology_information.df.empty:
            self._attach_chronology_ncgroup(site_group)
        self._attach_data_ncgroup(site_group)

    def to_netcdf(self, path_or_buffer):
        """Write NcdcRecord contents to a netCDF file
        """
        if isinstance(path_or_buffer, str):
            # Append to file, if it exists, if doesn't exist, create file.
            try:
                with netCDF4.Dataset(filename=path_or_buffer, mode='a', format='NETCDF4') as fl:
                    self._attach_ncgroups(fl)
            except FileNotFoundError:
                with netCDF4.Dataset(filename=path_or_buffer, mode='w', format='NETCDF4') as fl:
                    self._attach_ncgroups(fl)

        else:
            self._attach_ncgroups(path_or_buffer)


class QcPlotMixin:
    """Mixins to add QC plot methods LGM proxy records"""

    @staticmethod
    def _ax_setup(*args, **kwargs):
        try:
            import matplotlib.pylab as plt
        except ModuleNotFoundError:
            raise ModuleNotFoundError('matplotlib needs to be installed for plots')

        return plt.gca(*args, **kwargs)

    def plot_datavariable(self, variable, ax=None):
        """
        Plot a variable in the record data as timeseries.

        The plot also compares the variables from a redated age model with
        the original age model, if the record has been redated.

        Parameters
        ----------
        variable : str
            Name of variable to plot. Must be in ``self.data``.
        ax : :class:`mpl.axes.Axes` or None, optional
            Existing axes to plot onto.

        Returns
        -------
        ax : :class:`mpl.axes.Axes`

        """
        if ax is None:
            ax = self._ax_setup()

        if variable not in self.data.df.columns:
            raise KeyError('{} not found'.format(variable))

        if hasattr(self.data, 'age_median'):
            proxy_df = (self.data.df.set_index('depth')
                        .join(self.data.age_median, lsuffix='__', sort=True))

            new_age = proxy_df.loc[:, ('age_median', variable)].dropna()
            ax.plot(new_age.loc[:, 'age_median'], new_age.loc[:, variable], '.',
                    color='C3', label='MCMC median age')
            ax.plot(new_age.loc[:, 'age_median'], new_age.loc[:, variable],
                    color='C3', linewidth=0.5, label='_nolegend_')
            log.debug('Found new agemodel proxy timeseries')
        else:
            proxy_df = self.data.df.set_index('depth')
            log.debug('Assuming no new agemodel proxy timeseries')

        old_age = proxy_df.loc[:, ('age', variable)].dropna()
        ax.plot(old_age.loc[:, 'age'], old_age.loc[:, variable], 'x',
                color='C0', label='File age')
        ax.plot(old_age.loc[:, 'age'], old_age.loc[:, variable],
                color='C0', linewidth=0.5, label='_nolegend_')

        if 'd18o' in variable.lower():
            ax.invert_yaxis()
        ax.set_ylabel(variable)
        ax.grid()

        return ax

    def plot_sitemap(self, ax=None):
        """
        Plot sample site map using cartopy.

        Requires ``cartopy`` to be installed. Uses the site's northmost latitude
        and easternmost longitude to plot. So, we assume
        ``self.site_information.northernmost_latitude`` and
        ``self.site_information.easternmost_longitude`` are populated.

        Parameters
        ----------
        ax : :class:`mpl.axes.Axes` or None, optional
            Existing axes to plot onto.

        Returns
        -------
        ax : :class:`mpl.axes.Axes`

        """
        try:
            import cartopy.crs as ccrs
            import cartopy.feature as cfeature
        except ModuleNotFoundError:
            raise ModuleNotFoundError('cartopy needs to be installed for mapping')

        latlon = (self.site_information.northernmost_latitude,
                  self.site_information.easternmost_longitude)

        if ax is None:
            ax = self._ax_setup(projection=ccrs.Robinson(central_longitude=latlon[1]))

        ax.set_global()
        ax.add_feature(cfeature.LAND, facecolor='#B0B0B0')
        ax.outline_patch.set_linewidth(0.5)
        ax.plot(latlon[1], latlon[0], 'o', color='C0', transform=ccrs.Geodetic())

        return ax

    def plot_sedrate(self, ax=None):
        """
        Plot prior and posterior sediment rates for record agemodel.

        Requires ``self.chronology_information.bacon_agemodel` to be populated
        with a ``snakebacon.AgeDepthModel``-like instance. You can do this with
        :method:`self.redateredated()`, for example.

        Parameters
        ----------
        ax : :class:`mpl.axes.Axes` or None, optional
            Existing axes to plot onto.

        Returns
        -------
        ax : :class:`mpl.axes.Axes`

        """
        if ax is None:
            ax = self._ax_setup()

        ax = self.chronology_information.bacon_agemodel.plot_sediment_rate(ax)
        ax.lines[-1].set_color('C3')

        return ax

    def plot_sedmemory(self, ax=None):
        """
        Plot prior and posterior sediment memory for record agemodel.

        Requires ``self.chronology_information.bacon_agemodel` to be populated
        with a ``snakebacon.AgeDepthModel``-like instance. You can do this with
        :method:`self.redateredated()`, for example.

        Parameters
        ----------
        ax : :class:`mpl.axes.Axes` or None, optional
            Existing axes to plot onto.

        Returns
        -------
        ax : :class:`mpl.axes.Axes`

        """
        if ax is None:
            ax = self._ax_setup()

        ax = self.chronology_information.bacon_agemodel.plot_sediment_memory(ax)
        ax.lines[-1].set_color('C3')

        return ax

    def plot_agedepth(self, maxage=None, depthlines=None, prior_dwidth=30, ax=None):
        """
        Plot age models in relation to core depth.

        Parameters
        ----------
        maxage : float, int, or None, optional
            Cutoff age for the plot age model.
        depthlines : iterable or None, optional
            Depths to plot a vertical indicator line.
        prior_dwidth : int, optional
            Passed to :method:`snakebacon.AgeDepthModel.plot_prior_dates`.
        ax : :class:`mpl.axes.Axes` or None, optional
            Existing axes to plot onto.

        Returns
        -------
        ax : :class:`mpl.axes.Axes`

        """
        if ax is None:
            ax = self._ax_setup()

        agemodel = self.chronology_information.bacon_agemodel
        data_df = self.data.df
        # Hack to copy and crop the age model to a certain age.
        if maxage is not None:
            agemodel = deepcopy(agemodel)
            too_old = ~(agemodel.age_median() > maxage)
            agemodel._depth = agemodel.depth[too_old]
            agemodel._age_ensemble = agemodel.age_ensemble[too_old]

            # data_df = self.data.df.copy()
            # data_df = data_df.loc[data_df['age'] <= maxage, ('depth', 'age')]

        if len(agemodel._depth) > 0:
            # Skip age model plotting if maxage cut-out all samples.
            ax = agemodel.plot(ax=ax)
            ax.collections[-1].set_cmap('Greys')

            for l in ax.lines:
                l.set_color('C3')

            ax.plot(data_df.loc[:, 'depth'], data_df.loc[:, 'age'], 'C0',
                    label='File age model')
            ax = agemodel.plot_prior_dates(dwidth=prior_dwidth, ax=ax)
            ax.collections[-1].set_color('k')
            ax.collections[-1].set_zorder(10)

            ax.autoscale_view()
        else:
            log.warning('No data age-depth data to plot')

        try:
            for d in depthlines:
                ax.axvline(x=d, color='C4', linestyle='-.', zorder=1.5)
        except TypeError: # if depthlines is None
            pass

        ax.set_title('Age model')

        return ax

    def plot_deltar(self, ax=None):
        """
        Plot a description of the site carbon reservoir information.

        Parameters
        ----------
        ax : :class:`mpl.axes.Axes` or None, optional
            Existing axes to plot onto.

        Returns
        -------
        ax : :class:`mpl.axes.Axes`

        """
        if ax is None:
            ax = self._ax_setup()

        latlon = (self.site_information.northernmost_latitude,
                  self.site_information.easternmost_longitude)

        elevation = self.site_information.elevation

        # These checks are not well written.
        try:
            deltar_original = np.int(np.round(self.chronology_information.df['delta_R_original'].values[0]))
            if np.isnan(deltar_original):
                deltar_original = None
        except (KeyError, ValueError, IndexError) as e:
            deltar_original = None

        try:
            deltar_std_original = np.int(np.round(self.chronology_information.df['delta_R_1s_err_original'].values[0]))
            if np.isnan(deltar_std_original):
                deltar_std_original = None
        except (KeyError, ValueError, IndexError) as e:
            deltar_std_original = None

        try:
            deltar_used = np.int(np.round(self.chronology_information.df['delta_R'].values[0]))
        except (KeyError, TypeError):
            deltar_used = None

        try:
            deltar_error_used = np.int(np.round(self.chronology_information.df['delta_R_1s_err'].values[0]))
        except (KeyError, TypeError):
            deltar_error_used = None

        text_template = 'Latitude: {}°\nLongitude: {}°\nElevation: {} m ' \
                        '\n\nΔR: {}\nΔRσ: {}\nFile ΔR: {}\nFile ΔRσ: {}'
        text_str = text_template.format(latlon[0], latlon[1], elevation,
                                        deltar_used, deltar_error_used,
                                        deltar_original, deltar_std_original)

        ax.text(0.05, 0.9, text_str, verticalalignment='top',
                horizontalalignment='left', transform=ax.transAxes)
        ax.set_title('{}\n{}'.format(self.site_information.site_name, datetime.date.today().isoformat()))
        ax.axis('off')

        return ax

    def to_qcpdf(self, pdfpath, proxy_vars=None, plot_agedepth_kws=None):
        """
        Write quality-control report PDF.

        Requires :module:`matplotlib` and :module:`cartopy` to be installed.

        Parameters
        ----------
        pdfpath : str
            Path to write PDF to.
        proxy_vars : iterable or None, optional
            Name of series to include in time series plot. Attempts to use all
            available proxies if ``None``.
        plot_agedepth_kws : dict or None, optional
            Key-word arguments to pass to both call to ``self.plot_agedepth``
            when plotting.
        """
        log.debug('Writing QC report plots')

        try:
            from matplotlib.backends.backend_pdf import PdfPages
            import matplotlib.pylab as plt
        except ModuleNotFoundError:
            raise ModuleNotFoundError('matplotlib needs to be installed for plots')

        try:
            import cartopy.crs as ccrs
        except ModuleNotFoundError:
            raise ModuleNotFoundError('cartopy needs to be installed for mapping')

        if plot_agedepth_kws is None:
            plot_agedepth_kws = {}

        # Find "non-dimension" data variables to plot, if none were passed.
        not_proxy = ['age', 'age_median', 'depth', 'age_ensemble']
        if proxy_vars is None:
            # TODO(brews): This logic might be good candidate for a more general method.
            proxy_vars = [str(x) for x in self.data.df.columns if str(x).lower() not in not_proxy]

        latlon = (self.site_information.northernmost_latitude,
                  self.site_information.easternmost_longitude)

        n_vars = len(proxy_vars)
        n_cols = n_vars + 1

        has_baconagemodel = hasattr(self.chronology_information, 'bacon_agemodel')

        with PdfPages(pdfpath) as pdf:
            fig = plt.figure(figsize=(6.5, 9))

            ax2 = plt.subplot2grid((n_cols, 2), (0, 1),
                                   projection=ccrs.Robinson(central_longitude=latlon[1]))
            ax3 = plt.subplot2grid((n_cols, 2), (0, 0))

            self.plot_sitemap(ax=ax2)

            self.plot_deltar(ax=ax3)

            for i, varstr in enumerate(proxy_vars):
                this_ax = plt.subplot2grid((n_cols, 1), (1 + i, 0), colspan=2)
                this_ax = self.plot_datavariable(varstr, ax=this_ax)
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

            if has_baconagemodel:
                self.plot_sedrate(ax=ax1)
                self.plot_sedmemory(ax=ax2)
                self.plot_agedepth(maxage=50000, ax=ax4, **plot_agedepth_kws)

            fig.tight_layout()
            pdf.savefig(bbox_inches='tight')
            plt.close()

            fig = plt.figure(figsize=(6.5, 9))
            ax1 = plt.subplot(1, 1, 1)
            if has_baconagemodel:
                self.plot_agedepth(maxage=25000, ax=ax1, **plot_agedepth_kws)
            fig.tight_layout()
            pdf.savefig(bbox_inches='tight')
            plt.close()

            log.debug('QC report plot saved to {}'.format(pdfpath))
