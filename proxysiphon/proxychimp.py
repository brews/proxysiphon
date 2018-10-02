import logging
import collections
import datetime
from io import BytesIO
import pandas as pd
from chardet import detect as chdetect
import proxysiphon.records as records


DIVIDER = '#------------------------'
DATALINE_TRIGGER = '# Data lines follow (have no #)'
HEADINGS = ['Contribution_Date', 'Title', 'Investigators',
            'NOTE: Please cite original publication, online resource and date accessed when using this data.',
            'Description and Notes', 'Publication', 'Funding_Agency', 
            'Site Information', 'Data_Collection', 'Species', 
            'Chronology_Information', 'Variables', 'Data']
CHRON_HEADER = '# Labcode\tdepth_top\tdepth_bottom\tmat_dated\t14C_date\t14C_1s_err\tdelta_R\tdelta_R_1s_err\tother_date\tother_1s_err\tother_type\t'
MISSINGVALUE_LABEL = '# Missing Value: '


log = logging.getLogger(__name__)


def grab_latlon(g):
    """Grab latitude from proxychimp.Guts"""
    victim = g.pull_section('Site Information')[0]
    fields = dict(Northernmost_Latitude=None,
                  Southernmost_Latitude=None,
                  Easternmost_Longitude=None,
                  Westernmost_Longitude=None)
    for ln in victim:
        if any([x in ln for x in fields.keys()]):
            ln_split = ln.split(':')
            v = float(ln_split[-1])
            k = ln_split[0][1:].lstrip().rstrip()
            assert k in fields.keys()
            fields[k] = v
    return fields['Northernmost_Latitude'], fields['Westernmost_Longitude']

def grab_elevation(g):
    """Grab site elevation from proxychimp.Guts"""
    victim = g.pull_section('Site Information')[0]
    elevation = None
    for ln in victim:
        if 'Elevation:' in ln:
            ln_split = ln.split(':')
            elevation = float(ln_split[-1])
    return elevation

def grab_collection_year(g):
    """Get collection year from Data_Collection section of proxychimp.Guts"""
    victim = g.pull_section('Data_Collection')[0]
    yr = None
    for ln in victim:
        if 'Collection_Year' in ln:
            ln_split = ln.split(':')
            yr = int(ln_split[-1])
    return yr

def grab_publication_year(g):
    """Get publication year from proxychimp.Guts"""
    victim = g.pull_section('Publication')
    yr = []
    for p in victim:
        for ln in p:
            if 'Published_Date_or_Year' in ln:
                ln_split = ln.split(':')
                yr.append(int(ln_split[-1]))
    return yr

def grab_contribution_date(g):
    """Get contribution date from proxychimp.Guts"""
    victim = g.pull_section('Contribution_Date')[0]
    d = None
    for ln in victim:
        if 'Date' in ln:
            ln_split = ln.split(':')
            d_str = ln_split[-1].split('-')
            assert len(d_str) == 3
            d = datetime.date(int(d_str[0]), int(d_str[1]), int(d_str[2]))
    return d


def find_values(x, k, sep=':', fun=None):
    """Find values in a list of strings containing [k][sep][values]

    Parameters
    ----------
    x : iterable
        Iterable of strings representing lines in a file.
    k : str
        A key or pattern that to search x for.
    sep : str
        A separator that divides the `k` pattern from the target value.
    fun : function-like, optional
        Function used to format target value before returning.

    Returns
    -------
    The target value is returned after removing excess whitespace from left
    and right of the value. If the key or target pattern is not found,
    then None is returned. If an empty string, or whitespace is the found
    value, then None is returned.
    """
    out = None
    for l in x:
        if k in l and sep in l:
            val = l.split(sep, maxsplit=1)[1:][0].lstrip().rstrip()
            if val != '':
                out = val
    if fun is not None and out is not None:
        out = fun(out)
    return out


class Guts:
    """Ghetto and error-tolerant way to parse sections of NCDC proxy text files
    """

    def __init__(self, str_or_path):
        self._section_headings = ['# ' + s for s in HEADINGS]

        self.path = None
        self.encodingguess = None
        self.filestr = None
        try:
            with open(str_or_path, 'rb') as fl:
                flbytes = fl.read()
            self.path = str(str_or_path)
            self.encodingguess = chdetect(flbytes)
            # TODO(brews): Not sure we need ~both~ filestr, and lines + self.data + self.description.
            self.filestr = flbytes.decode(self.encodingguess['encoding'])
        except (FileNotFoundError, OSError):
            self.filestr = str(str_or_path)
        lines = self.filestr.splitlines()

        self.data = []
        self.description = []
        self.data_beginline = None
        self._divide_portions(lines)

        self.sectionindex = None
        self._write_sectionindex()

    def _divide_portions(self, lines):
        """Divide guts into 'description' and 'data' sections
        """
        prev_description = True
        dataline_flag = False
        for n, ln in enumerate(lines):
            if ln[0] == '#' or dataline_flag is False:
                self.description.append(ln)
                if DATALINE_TRIGGER in ln:
                    dataline_flag = True
                    log.debug('Found dataline flag on line {0}'.format(n))
            elif dataline_flag:
                self.data.append(ln)
                if prev_description is True:
                    self.data_beginline = n
                    log.debug('Data portion begins on line {0}'.format(n))
                prev_description = False

    def _write_sectionindex(self):
        """Populate the section index"""
        all_keys = []
        all_start = []
        all_stop = []
        prev_divider = False

        for idx, ln in enumerate(self.description):  # Lines after 6-line title.
            if prev_divider is True and any([h == ln.rstrip() for h in self._section_headings]):
                key = ln[1:].rstrip().lstrip()

                # If already have start idx for other section, append end idx
                # for that section
                if len(all_start) > 0:
                    all_stop.append(idx - 1)

                all_keys.append(key)
                all_start.append(idx)
                prev_divider = False
            if DIVIDER in ln:
                prev_divider = True
        all_stop.append(self.data_beginline)

        section_map = collections.defaultdict(list)
        for (k, start, stop) in zip(all_keys, all_start, all_stop):
            section_map[k].append((start, stop))
        section_map.default_factory = None
        self.sectionindex = section_map

    def pull_section(self, section):
        """Grab a list of list of line strings from the file description for each 'section'"""
        try:
            out = [self.description[slice(*idx)] for idx in self.sectionindex[section]]
        except KeyError:
            raise KeyError('section key "{}" not found'.format(section))
        return out

    def available_sections(self):
        """Get a list of available sections in the file"""
        return list(self.sectionindex.keys())

    def yank_data_df(self):
        """Get 'data' information as dataframe"""
        lines = [x.rstrip() for x in self.data]
        data_bytes = '\n'.join(lines).encode('utf-8')
        missingvalues = self.guess_missingvalues()
        df = pd.read_table(BytesIO(data_bytes), na_values=missingvalues)
        # df = pd.read_table(BytesIO(data_bytes), na_values=[-999, 'NaN'])
        return df

    def yank_chron_df(self, section_name='Chronology_Information', missingvalues=None):
        """Get chronology information as pandas.DataFrame"""
        if missingvalues is None:
            missingvalues = [-999, 'NaN']
        section = self.pull_section(section_name)[0]
        start_idx = section.index(CHRON_HEADER)
        g_chrond = section[start_idx:]
        g_chrond_cleaned = [x[2:].rstrip() for x in g_chrond]  # Removes the '# ' and ending white space.
        data_bytes = '\n'.join(g_chrond_cleaned).encode('utf-8')
        df = pd.read_table(BytesIO(data_bytes), na_values=missingvalues)
        return df

    def guess_missingvalues(self):
        """Guess data section missing values"""
        section = self.pull_section('Data')[0]
        out = None
        for ln in section:
            if MISSINGVALUE_LABEL in ln:
                out = ln.split(MISSINGVALUE_LABEL)[1]
        log.debug('Guessed missing value(s): {}'.format(out))
        return out

    def has_data(self):
        """Check if has populated data information"""
        try:
            d = self.yank_data_df()
        except KeyError:
            return False
        if len(d) > 0:
            result = True
        else:
            result = False
        return result

    def has_chron(self):
        """Check if has populated chronology information"""
        try:
            chron = self.yank_chron_df()
        except KeyError:
            return False
        if len(chron) > 0:
            result = True
        else:
            result = False
        return result

    def has_deltar(self):
        """Check if has populated delta R chronology information"""
        try:
            chron = self.yank_chron_df()
        except KeyError:
            return False
        if any(chron.delta_R.notnull()):
            result = True
        else:
            result = False
        return result

    def has_deltar_error(self):
        """Check if has populated delta R error chronology information"""
        try:
            chron = self.yank_chron_df()
        except KeyError:
            return False
        if any(chron.delta_R_1s_err.notnull()):
            result = True
        else:
            result = False
        return result

    def has_datacolumn(self, name):
        """Check if name is in data section columns"""
        try:
            data = self.yank_data_df()
        except KeyError:
            return False
        if name in data.columns:
            result = True
        else:
            result = False
        return result

    def yank_original_source_url(self):
        """Get string of original source URL

        Returns
        -------
        out : str or None

        Raises
        ------
        AssertionError
            If more than one section is found in the file data.
        """
        target_section = 'NOTE: Please cite original publication, online ' \
                         'resource and date accessed when using this data.'
        target_key = 'Original_Source_URL:'

        sections = self.pull_section(target_section)
        assert len(sections) < 2, 'More than one section found'
        section = sections[0]

        try:
            out = find_values(section, target_key, fun=str)
        except AttributeError:
            # target key not found in string
            out = None

        return out

    def yank_data_collection(self):
        """Get data collection information

        Returns
        -------
        out : dict or None

        Raises
        ------
        AssertionError
            If more than one section is found in the file data.
        """
        target_section = 'Data_Collection'
        # List of tuples, tuples give (dict_key, source_key, type_fun)
        target_keys = [('collection_name', 'Collection_Name:', str),
                       ('first_year', 'First_Year:', float),
                       ('last_year', 'Last_Year:', float),
                       ('time_unit', 'Time_Unit:', str),
                       ('core_length', 'Core_Length:', str),
                       ('notes', 'Notes:', str),
                       ('collection_year', 'Collection_Year:', int)]
        out = {k[0]:None for k in target_keys}

        sections = self.pull_section(target_section)
        assert len(sections) < 2, 'More than one section found'
        section = sections[0]

        for dict_key, source_key, type_fun in target_keys:
            out[dict_key] = find_values(section, source_key, fun=type_fun)

        return out

    def yank_description_and_notes(self):
        """Get string of description and notes

        Returns
        -------
        out : str or None

        Raises
        ------
        AssertionError
            If more than one section is found in the file data.
        """
        target_section = 'Description and Notes'
        target_key = 'Description:'

        sections = self.pull_section(target_section)
        assert len(sections) < 2, 'More than one section found'
        section = sections[0]

        try:
            out = find_values(section, target_key, fun=str)
        except AttributeError:
            # target key not found in string
            out = None

        return out

    def yank_publication(self):
        """Get list of publication information

        Returns
        -------
        out : list[dict]
        """
        target_section = 'Publication'
        # List of tuples, tuples give (dict_key, source_key, type_fun)
        target_keys = [('authors', '# Authors:', str),
                       ('published_date_or_year', '# Published_Date_or_Year:', int),
                       ('published_title', '# Published_Title:', str),
                       ('journal_name', '# Journal_Name:', str),
                       ('volume', '# Volume:', str),
                       ('edition', '# Edition:', str),
                       ('issue', '# Issue:', str),
                       ('pages', '# Pages:', str),
                       ('report_number', '# Report Number:', str),
                       ('doi', '# DOI:', str),
                       ('online_resource', '# Online_Resource:', str),
                       ('full_citation', '# Full_Citation:', str),
                       ('abstract', '# Abstract:', str)]
        dict_template = {k[0]:None for k in target_keys}

        out = []
        sections = self.pull_section(target_section)

        for section in sections:

            this_pub = dict_template.copy()
            for dict_key, source_key, type_fun in target_keys:
                this_pub[dict_key] = find_values(section, source_key, fun=type_fun)
            out.append(this_pub)

        return out

    def yank_variables(self):
        """Get variable section information

        Returns
        -------
        out : dict

        Raises
        ------
        AssertionError
            If more than one section is found in the file data.
        """
        target_section = 'Variables'
        sections = self.pull_section(target_section)
        assert len(sections) < 2, 'More than one section found'
        section = sections[0]

        out = {}
        for ln in section:

            # Skip line if not data variables line.
            if ln[:3] != '## ':
                continue

            var_name, components_group = ln[3:].split('\t')
            out[var_name] = tuple(components_group.split(','))

        return out

    def yank_site_information(self):
        """Get site information

        Returns
        -------
        out : dict

        Raises
        ------
        AssertionError
            If more than one section is found in the file data.
        """
        target_section = 'Site Information'
        # List of tuples, tuples give (dict_key, source_key, type_fun)
        target_keys = [('site_name', '# Site_Name:', str),
                       ('location', '# Location:', str),
                       ('country', '# Country:', str),
                       ('northernmost_latitude', '# Northernmost_Latitude:', float),
                       ('southernmost_latitude', '# Southernmost_Latitude:', float),
                       ('easternmost_longitude', '# Easternmost_Longitude:', float),
                       ('westernmost_longitude', '# Westernmost_Longitude:', float),
                       ('elevation', '# Elevation:', float)]
        out = {k[0]:None for k in target_keys}

        sections = self.pull_section(target_section)
        assert len(sections) < 2, 'More than one section found'
        section = sections[0]

        for dict_key, source_key, type_fun in target_keys:
            out[dict_key] = find_values(section, source_key, fun=type_fun)

        return out

    def to_ncdcrecord(self):
        """to NcdcRecord instance"""
        chron = records.ChronologyInformation(df=self.yank_chron_df())
        d = records.Data(df=self.yank_data_df())
        d_collection = records.DataCollection(**self.yank_data_collection())
        description = self.yank_description_and_notes()
        orig_url = self.yank_original_source_url()
        pubs = [records.Publication(**p) for p in self.yank_publication()]
        site_info = records.SiteInformation(**self.yank_site_information())

        vars = {}
        file_vars = self.yank_variables()
        for k,v in file_vars.items():
            vars[k] = records.VariableInfo(*v)

        out = records.NcdcRecord(chronology_information=chron,
                                 data=d,
                                 data_collection=d_collection,
                                 description=description,
                                 original_source_url=orig_url,
                                 publication=pubs,
                                 site_information=site_info,
                                 variables=vars)
        return out
