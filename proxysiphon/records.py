from dataclasses import dataclass, field
import unidecode
from chardet import detect as chdetect
import numpy as np
from pandas import DataFrame
import netCDF4
from proxysiphon.proxychimp import Guts


def read_ncdc(filepath_or_buffer, encoding=None):
    """Read NOAA NCDC txt file

    Parameters
    ----------
    filepath_or_buffer
    encoding : str or None, optional
        File encoding. Default is None which attempts to guess the encoding with
        `chardet.detect`.

    Returns
    -------
    out : NcdcRecord
    """
    with open(filepath_or_buffer, 'rb') as fl:
        flbytes = fl.read()
    if encoding is None:
        encoding = chdetect(flbytes)['encoding']
    g = Guts(flbytes.decode(encoding))
    return g.to_ncdcrecord()


def read_lgm(filepath_or_buffer, encoding=None):
    """Read NOAA NCDC txt file for LGM proxies

    Parameters
    ----------
    filepath_or_buffer
    encoding : str or None, optional
        File encoding. Default is None which attempts to guess the encoding with
        `chardet.detect`.

    Returns
    -------
    out : NcdcRecord
    """
    out = read_ncdc(filepath_or_buffer, encoding=None)
    return LgmRecord(**out.__dict__)


def read_petm(filepath_or_buffer, encoding=None):
    """Read NOAA NCDC txt file for PETM proxies

    Parameters
    ----------
    filepath_or_buffer
    encoding : str or None, optional
        File encoding. Default is None which attempts to guess the encoding with
        `chardet.detect`.

    Returns
    -------
    out : NcdcRecord
    """
    out = read_ncdc(filepath_or_buffer, encoding=None)
    return PetmRecord(**out.__dict__)


def _normalize_to_ascii_array(a, dtype='S50'):
    """Normalize sequence of UTF-8 string to np.Array of ASCII"""
    normed = np.array([unidecode.unidecode(str(x)) for x in a], dtype=dtype)
    return normed


@dataclass
class SiteInformation:
    """Proxy site information"""
    site_name: str = None
    location: str = None
    country: str = None
    northernmost_latitude: float = None
    southernmost_latitude: float = None
    easternmost_longitude: float = None
    westernmost_longitude: float = None
    elevation: float = None


@dataclass
class DataCollection:
    """Proxy site data collection information"""
    collection_name: str = None
    first_year: float = None
    last_year: float = None
    time_unit: str = None
    core_length: str = None
    notes: str = None
    collection_year: int = None


@dataclass
class VariableInfo:
    """Proxy site Data Variable information"""
    what: str
    material: str
    error: str
    units: str
    seasonality: str
    archive: str
    detail: str
    method: str
    datatype: str


@dataclass
class ChronologyInformation:
    """Proxy site chronology information"""
    df: DataFrame = field(default_factory=DataFrame)


@dataclass
class Data:
    """Proxy site data variables"""
    df: DataFrame = field(default_factory=DataFrame)


@dataclass
class Publication:
    """Proxy site publication"""
    authors: str = None
    published_date_or_year: int = None
    published_title: str = None
    journal_name: str = None
    volume: str = None
    edition: str = None
    issue: str = None
    pages: str = None
    report_number: str = None
    doi: str = None
    online_resource: str = None
    full_citation: str = None
    abstract: str = None

    def to_citationstr(self):
        """Citation str of publication"""
        if self.full_citation is not None:
            return str(self.full_citation)

        out = '{authors} ({year}): {title}.'.format(authors=self.authors,
                                                    year=self.published_date_or_year,
                                                    title=self.published_title)
        # This is a bit lazy...
        if self.journal_name is not None:
            out += ' {},'.format(self.journal_name)
        if self.edition is not None:
            out += ' {},'.format(self.edition)
        if self.volume is not None:
            out += ' {},'.format(self.volume)
        if self.issue is not None:
            out += ' {},'.format(self.issue)
        if self.pages is not None:
            out += ' {},'.format(self.pages)
        if self.report_number is not None:
            out += ' {},'.format(self.report_number)
        if self.doi is not None:
            out += ' doi:{}'.format(self.doi)
        if self.online_resource is not None:
            out += ' {}'.format(self.online_resource)

        return out


@dataclass
class NcdcRecord:
    """Proxy site NCDC record"""
    chronology_information: ChronologyInformation = None
    data: Data = None
    data_collection: DataCollection = None
    description: str = None
    original_source_url: str = None
    publication: list = field(default_factory=list)
    site_information: SiteInformation = None
    variables: dict = field(default_factory=dict)

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
        """Create the NetCDF4 site chronology group and populate it

        This is run whenever self.to_netcdf() is called.

        Parameters
        ----------
        parent : netcdf4.Group or netcdf4.Dataset
            Object the created chronology group will use as a parent.

        Returns
        -------
        chron : netcdf4.Group
            Reference to the created chronology group.
        """
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
        """Create the NetCDF4 site data group and populate it

        This is run whenever self.to_netcdf() is called.

        Parameters
        ----------
        parent : netcdf4.Group or netcdf4.Dataset
            Object the created data group will use as a parent.

        Returns
        -------
        data : netcdf4.Group
            Reference to the created chronology group.
        """
        data = parent.createGroup('data')
        data.createDimension('depth', None)
        depth = data.createVariable('depth', 'f4', ('depth',), zlib=True)
        depth.long_name = 'Sample depth'
        depth.positive = 'down'
        depth.axis = 'Z'
        depth[:] = self.data.df['depth'].values
        for k, v in self.variables['depth'].__dict__.items():
            setattr(depth, k, v)

        for col in list(self.data.df.columns):
            col_name = col.lower()

            if col_name == 'depth':
                continue

            var = data.createVariable(col_name, 'f4', ('depth',), zlib=True)
            var.missing_value = np.nan
            var[:] = self.data.df[col].values
            # Add more attributes to variable.
            for k, v in self.variables[col].__dict__.items():
                setattr(var, k, v)

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
            with netCDF4.Dataset(filename=path_or_buffer, mode='w', format='NETCDF4') as fl:
                self._attach_ncgroups(fl)
        else:
            self._attach_ncgroups(path_or_buffer)


class LgmRecord(NcdcRecord):
    """Proxy site LGM NCDC record"""

    @staticmethod
    def _variable_attributes(varname):
        """Variable attributes for proxy data

        This is an alternative to assigning metadata to self.data.df columns
        with self.variables. This is used to write data attributes when dumping
        to a NetCDF file.

        Parameters
        ----------
        varname : str
            String of data variables, as in self.data.df.columns.

        Returns
        -------
        out : dict
            A dictionary of metadata {'key': 'value'} for the input variable. An
            empty dictionary is returned if no information is available.
        """
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
        # Create and populate data group
        data = parent.createGroup('data')
        data.createDimension('depth', None)
        depth = data.createVariable('depth', 'f4', ('depth',), zlib=True)
        depth.long_name = 'Sample depth'
        depth.positive = 'down'
        depth.unints = 'cm'
        depth.axis = 'Z'
        depth[:] = self.data.df['depth'].values

        age_original = data.createVariable('age_original', 'f4', ('depth',),
                                           zlib=True)
        age_original.units = 'cal years BP'
        age_original.missing_value = np.nan
        age_original.long_name = 'Original age'
        age_original[:] = self.data.df['age'].values

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


class PetmRecord(LgmRecord):
    """Proxy site PETM NCDC record"""
    pass
