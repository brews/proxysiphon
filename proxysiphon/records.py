import datetime
from chardet import detect as chdetect

import proxysiphon.proxychimp as proxychimp


def read_ncdc(filepath_or_buffer, encoding=None):
    """Read NOAA NCDC txt file"""
    with open(filepath_or_buffer, 'rb') as fl:
        flbytes = fl.read()
    if encoding is None:
        encoding = chdetect(flbytes)['encoding']
    filestr = flbytes.decode(encoding)
    return NcdcRecord(filestr)


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


class SiteInformation:
    def __init__(self, s):
        self.site_name = None
        self.location = None
        self.country = None
        self.northernmost_latitude = None
        self.southernmost_latitude = None
        self.easternmost_longitude = None
        self.westernmost_longitude = None
        self.elevation = None

        s = str(s)
        s = s.splitlines()
        self.site_name = find_values(s, '# Site_Name:', fun=str)
        self.location = find_values(s, '# Location:', fun=str)
        self.country = find_values(s, '# Country:', fun=str)
        self.northernmost_latitude = find_values(s, '# Northernmost_Latitude:', fun=float)
        self.southernmost_latitude = find_values(s, '# Southernmost_Latitude:', fun=float)
        self.easternmost_longitude = find_values(s, '# Easternmost_Longitude:', fun=float)
        self.westernmost_longitude = find_values(s, '# Westernmost_Longitude:', fun=float)
        self.elevation = find_values(s, '# Elevation:', fun=float)


class DataCollection:
    def __init__(self, s):
        self.collection_name = None
        self.first_year = None
        self.last_year = None
        self.time_unit = None
        self.core_length = None
        self.notes = None
        self.collection_year = None

        s = str(s)
        s = s.splitlines()
        self.collection_name = find_values(s, 'Collection_Name:', fun=str)
        self.first_year = find_values(s, 'First_Year:', fun=float)
        self.last_year = find_values(s, 'Last_Year:', fun=float)
        self.time_unit = find_values(s, 'Time_Unit:', fun=str)
        self.core_length = find_values(s, 'Core_Length:', fun=str)
        self.notes = find_values(s, 'Notes:', fun=str)
        self.collection_year = find_values(s, 'Collection_Year:', fun=int)


class ChronologyInformation:
    def __init__(self, s):
        self.df = None
        if len(s) == 2:
            s_describe = s[0]
            self.df = s[1].copy()


class Data:
    def __init__(self, s):
        self.df = None
        if len(s) == 2:
            s_describe = s[0]
            self.df = s[1].copy()


class Publication:
    def __init__(self, s):
        self.doi = None

        s = str(s)
        s = s.splitlines()
        self.authors = find_values(s, '# Authors:', fun=str)
        self.published_date_or_year = find_values(s, '# Published_Date_or_Year:',
                                                  fun=int)
        self.published_title = find_values(s, '# Published_Title:', fun=str)
        self.journal_name = find_values(s, '# Journal_Name:', fun=str)
        self.volume = find_values(s, '# Volume:', fun=str)
        self.edition = find_values(s, '# Edition:', fun=str)
        self.issue = find_values(s, '# Issue:', fun=str)
        self.pages = find_values(s, '# Pages:', fun=str)
        self.report_number = find_values(s, '# Report Number:', fun=str)
        self.doi = find_values(s, '# DOI:', fun=str)
        self.online_resource = find_values(s, '# Online_Resource:', fun=str)
        self.full_citation = find_values(s, '# Full_Citation:', fun=str)
        self.abstract = find_values(s, '# Abstract:', fun=str)


class ContributionDate:
    def __init__(self, s):
        self.date = None

        s = str(s)
        s = s.splitlines()
        try:
            v = find_values(s, 'Date:', fun=str).split('-')
            assert len(v) == 3
            self.date = datetime.date(year=int(v[0]), month=int(v[1]), day=int(v[2]))
        except AttributeError:
            # Date information not found, as v is None, as is self.date
            pass


class Description:
    def __init__(self, s):
        self.description = None

        s = str(s)
        s = s.splitlines()
        try:
            self.description = find_values(s, 'Description:', fun=str)
        except AttributeError:
            # Description information  not found, should be None.
            pass


class SourceUrl:
    def __init__(self, s):
        self.original_source_url = None

        s = str(s)
        s = s.splitlines()
        try:
            self.original_source_url = find_values(s, 'Original_Source_URL:', fun=str)
        except AttributeError:
            # Source URL information  not found, should be None.
            pass


class NcdcRecord:
    def __init__(self, s):
        g = proxychimp.Guts(s)

        def stringer(guts, sec, cls=None):
            out = []
            try:
                secs = guts.pull_section(sec)
                for sec in secs:
                    if cls is not None:
                        this_data = cls('\n'.join(sec))
                    else:
                        this_data = '\n'.join(sec)
                    out.append(this_data)
            except KeyError:
                out.append(None)
            return out

        self.original_source_url = stringer(g, 'NOTE: Please cite original publication, online resource and date accessed when using this data.', SourceUrl)
        self.contribution_date = ContributionDate(stringer(g, 'Contribution_Date')[0])
        self.description = Description(stringer(g, 'Description_and_Notes')[0])
        self.publication = stringer(g, 'Publication', Publication)
        self.site_information = SiteInformation(stringer(g, 'Site Information')[0])
        self.chronology_information = ChronologyInformation([stringer(g, 'Chronology_Information')[0], g.yank_chron_df()])
        self.data = Data([stringer(g, 'Data')[0], g.yank_data_df()])
        self.data_collection = DataCollection(stringer(g, 'Data_Collection')[0])
