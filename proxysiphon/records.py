from dataclasses import dataclass, field
from chardet import detect as chdetect
from pandas import DataFrame
from proxysiphon.proxychimp import Guts
import proxysiphon.lgm as lgm


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


class LgmRecord(lgm.QcPlotMixin, lgm.RedateMixin, lgm.NetcdfMixin, NcdcRecord):
    """Proxy site LGM NCDC record"""
    pass


class PetmRecord(LgmRecord):
    """PETM proxy site NCDC record"""
    pass
