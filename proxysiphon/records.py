from chardet import detect as chdetect

import proxysiphon.proxychimp as proxychimp


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
    """
    Attributes
    ----------
    site_name : str or None
    location : str or None
    country : str or None
    northernmost_latitude : float or None
    southernmost_latitude : float or None
    easternmost_longitude : float or None
    westernmost_longitude : float or None
    elevation : float or None
    """
    def __init__(self, s=None, **kwargs):
        """
        Parameters
        ----------
        s : str or None, optional
            A string of 'Site Information' section from NOAA
            NCDC file. If `None`, then information is taken from `**kwargs`.
        **kwargs : optional
            Keywords value pairs are assigned to object attributes. Overwrites
            values in `s`.
        """
        self.site_name = None
        self.location = None
        self.country = None
        self.northernmost_latitude = None
        self.southernmost_latitude = None
        self.easternmost_longitude = None
        self.westernmost_longitude = None
        self.elevation = None

        if s is not None:
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

        allowed_keys = ['site_name', 'location', 'country', 'northernmost_latitude',
                        'southernmost_latitude', 'easternmost_longitude',
                        'westernmost_longitude', 'elevation']
        self.__dict__.update((k, v) for k, v in kwargs.items() if k in allowed_keys)


class DataCollection:
    """
    Attributes
    ----------
    collection_name : str or None
    first_year : float or None
    last_year : float or None
    time_unit : str or None
    core_length : str or None
    notes : str or None
    collection_year int or None
    """
    def __init__(self, s=None, **kwargs):
        """
        Parameters
        ----------
        s : str or None, optional
            'Data Collection' section string from NOAA NCDC file.  If `None`,
            then information is taken from `**kwargs`.
        **kwargs : optional
            Keywords value pairs are assigned to object attributes. Overwrites
            values in `s`.
        """
        self.collection_name = None
        self.first_year = None
        self.last_year = None
        self.time_unit = None
        self.core_length = None
        self.notes = None
        self.collection_year = None

        if s is not None:
            s = str(s)
            s = s.splitlines()
            self.collection_name = find_values(s, 'Collection_Name:', fun=str)
            self.first_year = find_values(s, 'First_Year:', fun=float)
            self.last_year = find_values(s, 'Last_Year:', fun=float)
            self.time_unit = find_values(s, 'Time_Unit:', fun=str)
            self.core_length = find_values(s, 'Core_Length:', fun=str)
            self.notes = find_values(s, 'Notes:', fun=str)
            self.collection_year = find_values(s, 'Collection_Year:', fun=int)

        allowed_keys = ['collection_name', 'first_year', 'last_year',
                        'time_unit', 'core_length', 'notes', 'collection_year']
        self.__dict__.update((k, v) for k, v in kwargs.items() if k in allowed_keys)


class ChronologyInformation:
    """
    Attributes
    ----------
    df : pandas.DataFrame or None
        DataFrame of chronology information.
    """
    def __init__(self, s=None, **kwargs):
        """
        Parameters
        ----------
        s : sequence or None, optional
            First member is 'Chronology Information' section string from NOAA
            NCDC file. Second is a pandas.DataFrame of chronology information
            table.
        **kwargs : optional
            Keywords value pairs are assigned to object attributes. Overwrites
            values in `s`.
        """
        self.df = None

        if s is not None:
            assert len(s) == 2, 'len(s) must equal 2'
            # s_describe = s[0]
            self.df = s[1].copy()

        allowed_keys = ['df']
        self.__dict__.update((k, v) for k, v in kwargs.items() if k in allowed_keys)


class Data:
    """
    Attributes
    ----------
    df : pandas.DataFrame or None
        DataFrame of core data.
    """
    def __init__(self, s=None, **kwargs):
        """
        Parameters
        ----------
        s : sequence or None, optional
            First member is 'Data' section string from NOAA NCDC file. Second
            is a pandas.DataFrame of proxy data table.
        **kwargs : optional
            Keywords value pairs are assigned to object attributes. Overwrites
            values in `s`.
        """
        self.df = None

        if s is not None:
            assert len(s) == 2, 'len(s) must equal 2'
            # s_describe = s[0]
            self.df = s[1].copy()

        allowed_keys = ['df']
        self.__dict__.update((k, v) for k, v in kwargs.items() if k in allowed_keys)


class Publication:
    """
    Attributes
    ----------
    authors : str or None
    published_date_or_year : int or None
    published_title: str or None
    journal_name : str or None
    volume : str or None
    edition : str or None
    issue : str or None
    pages : str or None
    report_number : str or None
    doi : str or None
    online_resource : str or None
    full_citation : str or None
    abstract : str or None
    """
    def __init__(self, s=None, **kwargs):
        """
        Parameters
        ----------
        s : str
            Single 'Publication' section string from NOAA NCDC file.
        **kwargs : optional
            Keywords value pairs are assigned to object attributes. Overwrites
            values in `s`.
        """
        self.authors = None
        self.published_date_or_year = None
        self.published_title = None
        self.journal_name = None
        self.volume = None
        self.edition = None
        self.issue = None
        self.pages = None
        self.report_number = None
        self.doi = None
        self.online_resource = None
        self.full_citation = None
        self.abstract = None

        if s is not None:
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

        allowed_keys = ['authors', 'published_date_or_year', 'published_title',
                        'journal_name', 'volume', 'edition', 'issue', 'pages',
                        'report_number', 'doi', 'online_resource',
                        'full_citation', 'abstract']
        self.__dict__.update((k, v) for k, v in kwargs.items() if k in allowed_keys)


def fetch_description(s):
    """Grab 'Description' string from section string `s`
    """
    s = str(s)
    s = s.splitlines()

    out = None

    try:
        out = find_values(s, 'Description:', fun=str)
    except AttributeError:
        # Description information  not found, should be None.
        pass
    return out


def fetch_sourceurl(s):
    """Grab 'Original_Source_URL' string from section string `s`
    """
    s = str(s)
    s = s.splitlines()

    out = None

    try:
        out = find_values(s, 'Original_Source_URL:', fun=str)
    except AttributeError:
        # Source URL information not found, should be None.
        pass
    return out


class NcdcRecord:
    """

    Attributes
    ----------
    original_source_url : sequence
        Sequence of strings giving URLs for data source.
    description : str or None
        String with additional description and notes.
    publication : sequence of Publication
        Publications citing proxy material. Empty sequence if no publications
        given.
    site_information : SiteInformation
        Information on proxy site.
    chronology_information : ChronologyInformation
        Time-constraints on proxy.
    data : Data
        Actual proxy data for site.
    data_collection : DataCollection
        Notes and additional information on collection name and when the data
        was collected.
    """
    def __init__(self, s):
        """
        Parameters
        ----------
        s : str
            String from NOAA NCDC file.
        """
        g = proxychimp.Guts(s)

        def stringer(guts, sec, cls=None):
            # Should prob move this func outside of the class constructor.
            out = []
            secs = guts.pull_section(sec)
            for sec in secs:
                if cls is not None:
                    this_data = cls('\n'.join(sec))
                else:
                    this_data = '\n'.join(sec)
                out.append(this_data)
            return out

        self.original_source_url = stringer(g, 'NOTE: Please cite original publication, online resource and date accessed when using this data.', fetch_sourceurl)
        self.description = stringer(g, 'Description and Notes', fetch_description)[0]
        self.publication = stringer(g, 'Publication', Publication)
        self.site_information = SiteInformation(stringer(g, 'Site Information')[0])
        self.chronology_information = ChronologyInformation([stringer(g, 'Chronology_Information')[0], g.yank_chron_df()])
        self.data = Data([stringer(g, 'Data')[0], g.yank_data_df()])
        self.data_collection = DataCollection(stringer(g, 'Data_Collection')[0])
