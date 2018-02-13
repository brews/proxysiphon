import logging
import datetime
from io import BytesIO

import pandas as pd
from chardet import detect as chdetect


DIVIDER = '#------------------------'
DATALINE_TRIGGER = '# Data lines follow (have no #)'
HEADINGS = ['Contribution_Date', 'Title', 'Investigators', 
            'Description and Notes', 'Publication', 'Funding_Agency', 
            'Site Information', 'Data_Collection', 'Species', 
            'Chronology_Information', 'Variables', 'Data']
CHRON_HEADER = '# Labcode\tdepth_top\tdepth_bottom\tmat_dated\t14C_date\t14C_1s_err\tdelta_R\tdelta_R_1s_err\tother_date\tother_1s_err\tother_type\t'
MISSINGVALUE_LABEL = '# Missing Value: '


log = logging.getLogger(__name__)


def grab_latlon(g):
    """Grab latitude from proxychimp.Guts"""
    victim = g.pull_section('Site Information')
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
    victim = g.pull_section('Site Information')
    elevation = None
    for ln in victim:
        if 'Elevation:' in ln:
            ln_split = ln.split(':')
            elevation = float(ln_split[-1])
    return elevation

def grab_collection_year(g):
    """Get collection year from Data_Collection section of proxychimp.Guts"""
    victim = g.pull_section('Data_Collection')
    yr = None
    for ln in victim:
        if 'Collection_Year' in ln:
            ln_split = ln.split(':')
            yr = int(ln_split[-1])
    return yr

def grab_publication_year(g):
    """Get publication year from proxychimp.Guts"""
    victim = g.pull_section('Publication')
    yr = None
    for ln in victim:
        if 'Published_Date_or_Year' in ln:
            ln_split = ln.split(':')
            yr = int(ln_split[-1])
    return yr

def grab_contribution_date(g):
    """Get contribution date from proxychimp.Guts"""
    victim = g.pull_section('Contribution_Date')
    d = None
    for ln in victim:
        if 'Date' in ln:
            ln_split = ln.split(':')
            d_str = ln_split[-1].split('-')
            assert len(d_str) == 3
            d = datetime.date(int(d_str[0]), int(d_str[1]), int(d_str[2]))
    return d

def plot_sitemap(g):
    """Plot site location on global map from proxychimp.Guts"""
    pass


def plot_agedepth(g):
    """Plot given age-depth relationship from proxychimp.Guts data"""
    pass



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
        prev_divider = False
        for idx, ln in enumerate(self.description):  # Lines after 6-line title.
            if prev_divider is True and any([h == ln.rstrip() for h in self._section_headings]):
                key = ln[1:].rstrip().lstrip()
                all_keys.append(key)
                all_start.append(idx)
                prev_divider = False
            if DIVIDER in ln:
                prev_divider = True
        all_stop = [x - 1 for x in all_start[1:]]
        all_stop.append(self.data_beginline)
        self.sectionindex = {k: (v0, v1) for (k, v0, v1) in zip(all_keys, all_start, all_stop)}

    def pull_section(self, section):
        """Grab a list of line strings from the file description"""
        return self.description[slice(*self.sectionindex[section])]

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
        """Get chronology information as dataframe"""
        if missingvalues is None:
            missingvalues = [-999, 'NaN']
        section = self.pull_section(section_name)
        start_idx = section.index(CHRON_HEADER)
        g_chrond = section[start_idx:]
        g_chrond_cleaned = [x[2:].rstrip() for x in g_chrond]  # Removes the '# ' and ending white space.
        data_bytes = '\n'.join(g_chrond_cleaned).encode('utf-8')
        df = pd.read_table(BytesIO(data_bytes), na_values=missingvalues)
        return df

    def guess_missingvalues(self):
        """Guess data section missing values"""
        section = self.pull_section('Data')
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
