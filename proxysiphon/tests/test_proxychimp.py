from tempfile import NamedTemporaryFile
import pytest
import pandas as pd

from proxysiphon import proxychimp


datapayload = ['#------------------------',
           '# Contribution_Date',
           '#',
           '#------------------------',
           '# Title',
           '#',
           '#------------------------',
           '# Chronology_Information',
           '#',
           '# Labcode\tdepth_top\tdepth_bottom\tmat_dated\t14C_date\t14C_1s_err\tdelta_R\tdelta_R_1s_err\tother_date\tother_1s_err\tother_type\t',
           '# 152757	5	6	G. ruber or mixed planktonic	475	30	188	73	-999	-999	-999	',
           '#',
           '#------------------------',
           '# Data',
           '# Data lines follow (have no #)',
           '# Missing Value: -999',
           'depth\tage\tbacon',
           '1\t2\t3',
           '4\t5\t6']


@pytest.fixture(scope='module')
def chron_nodeltaR_nodata_guts():
    payload = ['#------------------------',
               '# Contribution_Date',
               '#',
               '#------------------------',
               '# Title',
               '#',
               '#------------------------',
               '# Chronology_Information',
               '#',
               '# Labcode\tdepth_top\tdepth_bottom\tmat_dated\t14C_date\t14C_1s_err\tdelta_R\tdelta_R_1s_err\tother_date\tother_1s_err\tother_type\t',
               '# 152757	5	6	G. ruber or mixed planktonic	-999	-999	-999	-999	-999	-999	-999	',
               '#',
               '#------------------------',
               '# Data',
               '# Data lines follow (have no #)',
               '# Missing Value: -999',
               'depth\tage\tbacon']
    with NamedTemporaryFile('wb') as tf:
        tf.write('\n'.join(payload).encode('utf-8'))
        tf.flush()
        g = proxychimp.Guts(tf.name)
    return g


@pytest.fixture(scope='module')
def chron_guts():
    payload = datapayload
    with NamedTemporaryFile('wb') as tf:
        tf.write('\n'.join(payload).encode('utf-8'))
        tf.flush()
        g = proxychimp.Guts(tf.name)
    return g


@pytest.fixture(scope='module')
def dumb_guts():
    payload = ['#------------------------',
               '# NOTE: Please cite original publication, online resource and date accessed when using this data.',
               '# If there is no publication information, please cite Investigator,', '#',
               '# Description/Documentation lines begin with #', '# Data lines have no #', '#',
               '# Online_Resource: http://www.ncdc.noaa.gov/paleo/study/',
               '# Online_Resource: http://www1.ncdc.noaa.gov/pub/data/paleo/', '#',
               '# Original_Source_URL: https://www.ncdc.noaa.gov/paleo-search/study/2622',
               '#------------------------',
               '# Contribution_Date',
               '#',
               '#------------------------',
               '# Title',
               '#',
               '#------------------------',
               '# Data_Collection', '#   Collection_Name: P178-15P',
               '#   First_Year: 39485', '#   Last_Year: -18',
               '#   Time_Unit: cal yr BP', '#   Core_Length: ',
               '#   Notes: mg_red', '#   Collection_Year: 1923',
               '#------------------------',
               '# Site Information', '# Site_Name: P178-15P',
               '# Location: Arabian Sea', '# Country: ',
               '# Northernmost_Latitude: 11.955', '# Southernmost_Latitude: 11.955',
               '# Easternmost_Longitude: 44.3', '# Westernmost_Longitude: 44.3',
               '# Elevation: -869',
               '#------------------------',
               '# Description and Notes',
               '#        Description: d18O sacc from 2003 paper, mg/ca sacc from 2002 paper, alkeno',
               '#------------------------',
               '# Publication',
               '# Authors: Tierney, Jessica E.; Pausata, Francesco S. R.; deMenocal, Peter B.',
               '# Published_Date_or_Year: 2016',
               '# Published_Title: Deglacial Indian monsoon failure and North Atlantic stadials linked by Indian Ocean surface cooling',
               '# Journal_Name: Nature Geoscience', '# Volume: 9',
               '# Edition: ', '# Issue: ', '# Pages: 46-50',
               '# Report Number: ', '# DOI: 10.1038/ngeo2603',
               '# Online_Resource: ', '# Full_Citation: ', '# Abstract:',
               '#------------------------',
               '# Chronology_Information',
               '#',
               '# Labcode\tdepth_top\tdepth_bottom\tmat_dated\t14C_date\t14C_1s_err\tdelta_R\tdelta_R_1s_err\tother_date\tother_1s_err\tother_type\t',
               '#',
               '#---------------------------------------',
               '# Variables',
               '# Data variables follow that are preceded by "##" in columns one and two.',
               '# Variables list, one per line, shortname-tab-longname components (9 components: what, material, error, units, seasonality, archive, detail, method, C or N for Character or Numeric data)',
               '## depth	,,,cm,,,,,',
               '## age	,,,cal yr BP,,,,,',
               '## bacon	,,,index,,,,,',
               '#------------------------',
               '# Data',
               '# Data lines follow (have no #)',
               '# Missing Value: -999',
               'depth\tage\tbacon',
               '1\t2\t3',
               '4\t5\t6']
    with NamedTemporaryFile('wb') as tf:
        tf.write('\n'.join(payload).encode('utf-8'))
        tf.flush()
        g = proxychimp.Guts(tf.name)
    return g


def test_find_values():
    lines = ['apple: 1\n', 'bee: 1:2\n', 'charlie: 1.5\n', 'dingo-2\t\n', 'echo: ']
    assert proxychimp.find_values(lines, 'NotAKey') is None
    assert proxychimp.find_values(lines, 'apple') == '1'
    assert proxychimp.find_values(lines, 'bee') == '1:2'
    assert proxychimp.find_values(lines, 'charlie', fun=float) == 1.5
    assert proxychimp.find_values(lines, 'bee', sep='-', fun=int) is None
    assert proxychimp.find_values(lines, 'dingo', sep='-', fun=int) == 2
    assert proxychimp.find_values(lines, 'echo') is None


def test_str_guts__init_():
    filestr = '\n'.join(datapayload)
    g = proxychimp.Guts(filestr)
    goal_dict = {'Labcode': [152757], 'depth_top': [5], 'depth_bottom': [6],
                 'mat_dated': ['G. ruber or mixed planktonic'],
                 '14C_date': [475],'14C_1s_err': [30],
                 'delta_R': [188], 'delta_R_1s_err': [73],
                 'other_date': [pd.np.nan], 'other_1s_err': [pd.np.nan],
                 'other_type': [pd.np.nan]}
    goal = pd.DataFrame(goal_dict)
    df = g.yank_chron_df()
    for k in goal_dict.keys():
        pd.testing.assert_series_equal(goal[k], df[k])


def test__divide_portions(dumb_guts):
    goal_data = ['depth\tage\tbacon', '1\t2\t3', '4\t5\t6']
    goal_description_ends = ('#------------------------', '# Missing Value: -999')
    goal_beginline = 70
    assert goal_data == dumb_guts.data
    assert goal_description_ends[0] == dumb_guts.description[0]
    assert goal_description_ends[-1] == dumb_guts.description[-1]
    assert goal_beginline == dumb_guts.data_beginline


def test__write_sectionindex(dumb_guts):
    goal_index = {'Chronology_Information': [(55, 59)],
                  'Contribution_Date': [(12, 14)],
                  'Data': [(67, 70)],
                  'Title': [(15, 17)]}
    assert 10 == len(dumb_guts.sectionindex.keys())
    for k, goal in goal_index.items():
        assert goal == dumb_guts.sectionindex[k]


def test_pull_section(dumb_guts):
    goal = ['# Data', '# Data lines follow (have no #)',
            '# Missing Value: -999']
    assert dumb_guts.pull_section('Data')[0] == goal


def test_pull_section_badname(dumb_guts):
    with pytest.raises(KeyError):
        dumb_guts.pull_section('data')

def test_available_sections(dumb_guts):
    goal = ['Contribution_Date', 'Title', 'Data_Collection',
            'Site Information', 'Description and Notes', 'Publication',
            'Chronology_Information', 'Variables', 'Data']
    # Skip first section as it is very long:
    assert goal == dumb_guts.available_sections()[1:]


def test_guess_missingvalues(dumb_guts):
    goal = '-999'
    assert goal == dumb_guts.guess_missingvalues()


def test_yank_data_df(dumb_guts):
    goal_dict = {'depth': [1, 4], 'age': [2, 5], 'bacon': [3, 6]}
    goal = pd.DataFrame(goal_dict)
    df = dumb_guts.yank_data_df()
    for k in goal_dict.keys():
        pd.testing.assert_series_equal(goal[k], df[k])


def test_yank_chron_df(chron_guts):
    goal_dict = {'Labcode': [152757], 'depth_top': [5], 'depth_bottom': [6],
                 'mat_dated': ['G. ruber or mixed planktonic'],
                 '14C_date': [475],'14C_1s_err': [30],
                 'delta_R': [188], 'delta_R_1s_err': [73],
                 'other_date': [pd.np.nan], 'other_1s_err': [pd.np.nan],
                 'other_type': [pd.np.nan]}
    goal = pd.DataFrame(goal_dict)
    df = chron_guts.yank_chron_df()
    for k in goal_dict.keys():
        pd.testing.assert_series_equal(goal[k], df[k])


def test_yank_publication(dumb_guts):
    v1 = dumb_guts.yank_publication()[0]
    assert 'Tierney, Jessica E.; Pausata, Francesco S. R.; deMenocal, Peter B.' == v1['authors']
    assert '10.1038/ngeo2603' == v1['doi']


def test_yank_site_infromation(dumb_guts):
    v = dumb_guts.yank_site_information()
    assert 'P178-15P' == v['site_name']
    assert 'Arabian Sea' == v['location']
    assert v['country'] is None
    assert 11.955 == v['northernmost_latitude']
    assert 11.955 == v['southernmost_latitude']
    assert 44.3 == v['easternmost_longitude']
    assert 44.3 == v['westernmost_longitude']
    assert -869 == v['elevation']


def test_yank_data_collection(dumb_guts):
    v = dumb_guts.yank_data_collection()
    assert 'P178-15P' == v['collection_name']
    assert 39485 == v['first_year']
    assert -18 == v['last_year']
    assert 'cal yr BP' == v['time_unit']
    assert v['core_length'] is None
    assert 'mg_red' == v['notes']
    assert 1923 == v['collection_year']


def test_pull_variables(dumb_guts):
    v = dumb_guts.yank_variables()
    assert 'cm' == v['depth'][3]


def test_yank_description_and_notes(dumb_guts):
    v = dumb_guts.yank_description_and_notes()
    assert 'd18O sacc from 2003 paper, mg/ca sacc from 2002 paper, alkeno' == v


def test_yank_original_source_url(dumb_guts):
    v = dumb_guts.yank_original_source_url()
    assert 'https://www.ncdc.noaa.gov/paleo-search/study/2622' == v


def test_has_chron(dumb_guts, chron_guts):
    assert dumb_guts.has_chron() is False
    assert chron_guts.has_chron() is True


def test_has_deltar(dumb_guts, chron_guts, chron_nodeltaR_nodata_guts):
    assert chron_nodeltaR_nodata_guts.has_deltar() is False
    assert dumb_guts.has_deltar() is False
    assert chron_guts.has_deltar() is True


def test_has_datacolumn(dumb_guts, chron_nodeltaR_nodata_guts):
    assert chron_nodeltaR_nodata_guts.has_datacolumn('boogers') is False
    assert dumb_guts.has_datacolumn('depth') is True
