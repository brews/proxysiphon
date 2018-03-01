import pytest
import pandas as pd

from proxysiphon import records, proxychimp


@pytest.fixture(scope='module')
def testchrondatalines():
    lines = ['#------------------------',
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
    return '\n'.join(lines)


def test_find_values():
    lines = ['apple: 1\n', 'bee: 1:2\n', 'charlie: 1.5\n', 'dingo-2\t\n', 'echo: ']
    assert records.find_values(lines, 'NotAKey') is None
    assert records.find_values(lines, 'apple') == '1'
    assert records.find_values(lines, 'bee') == '1:2'
    assert records.find_values(lines, 'charlie', fun=float) == 1.5
    assert records.find_values(lines, 'bee', sep='-', fun=int) is None
    assert records.find_values(lines, 'dingo', sep='-', fun=int) == 2
    assert records.find_values(lines, 'echo') is None


def test_SiteInformation():
    testlines = """# Site Information
                # Site_Name: P178-15P
                # Location: Arabian Sea
                # Country: 
                # Northernmost_Latitude: 11.955
                # Southernmost_Latitude: 11.955
                # Easternmost_Longitude: 44.3
                # Westernmost_Longitude: 44.3
                # Elevation: -869"""
    v = records.SiteInformation(testlines)
    assert v.site_name == 'P178-15P'
    assert v.location == 'Arabian Sea'
    assert v.country is None
    assert v.northernmost_latitude == 11.955
    assert v.southernmost_latitude == 11.955
    assert v.easternmost_longitude == 44.3
    assert v.westernmost_longitude == 44.3
    assert v.elevation == -869


def test_DataCollection():
    testlines = """# Data_Collection
                #   Collection_Name: P178-15P
                #   First_Year: 39485
                #   Last_Year: -18
                #   Time_Unit: cal yr BP
                #   Core_Length: 
                #   Notes: mg_red
                #   Collection_Year: 1923"""
    v = records.DataCollection(testlines)
    assert v.collection_name == 'P178-15P'
    assert v.first_year == 39485
    assert v.last_year == -18
    assert v.time_unit == 'cal yr BP'
    assert v.core_length is None
    assert v.notes == 'mg_red'
    assert v.collection_year == 1923


def test_ChronologyInformation(testchrondatalines):
    goal_dict = {'Labcode': [152757], 'depth_top': [5], 'depth_bottom': [6],
                 'mat_dated': ['G. ruber or mixed planktonic'],
                 '14C_date': [475],'14C_1s_err': [30],
                 'delta_R': [188], 'delta_R_1s_err': [73],
                 'other_date': [pd.np.nan], 'other_1s_err': [pd.np.nan],
                 'other_type': [pd.np.nan]}
    goal = pd.DataFrame(goal_dict)
    g = proxychimp.Guts(testchrondatalines)
    df = records.ChronologyInformation(['\n'.join(g.pull_section('Data')), g.yank_chron_df()]).df
    for k in goal_dict.keys():
        pd.testing.assert_series_equal(goal[k], df[k])


def test_Data(testchrondatalines):
    goal_dict = {'depth': [1, 4], 'age': [2, 5], 'bacon': [3, 6]}
    goal = pd.DataFrame(goal_dict)
    g = proxychimp.Guts(testchrondatalines)
    df = records.Data(['\n'.join(g.pull_section('Data')), g.yank_data_df()]).df
    for k in goal_dict.keys():
        pd.testing.assert_series_equal(goal[k], df[k])


def test_Publication():
    testlines = """# Publication
                # Authors: Tierney, Jessica E.; Pausata, Francesco S. R.; deMenocal, Peter B.
                # Published_Date_or_Year: 2016
                # Published_Title: Deglacial Indian monsoon failure and North Atlantic stadials linked by Indian Ocean surface cooling
                # Journal_Name: Nature Geoscience
                # Volume: 9
                # Edition: 
                # Issue: 
                # Pages: 46-50
                # Report Number: 
                # DOI: 10.1038/ngeo2603
                # Online_Resource: 
                # Full_Citation: 
                # Abstract: """
    v = records.Publication(testlines)
    assert v.doi == '10.1038/ngeo2603'


def test_ContributionDate():
    testlines = """# Contribution_Date
                #    Date: 2017-03-31"""
    v = records.ContributionDate(testlines)
    assert v.date.year == 2017
    assert v.date.month == 3
    assert v.date.day == 31


def test_NcdcRecord_chron_data(testchrondatalines):
    v = records.NcdcRecord(testchrondatalines)
    goal_datadict = {'depth': [1, 4], 'age': [2, 5], 'bacon': [3, 6]}
    goal_data = pd.DataFrame(goal_datadict)
    df_data = v.data.df.copy()
    for k in goal_datadict.keys():
        pd.testing.assert_series_equal(goal_data[k], df_data[k])

    goal_chrondict = {'Labcode': [152757], 'depth_top': [5], 'depth_bottom': [6],
                 'mat_dated': ['G. ruber or mixed planktonic'],
                 '14C_date': [475],'14C_1s_err': [30],
                 'delta_R': [188], 'delta_R_1s_err': [73],
                 'other_date': [pd.np.nan], 'other_1s_err': [pd.np.nan],
                 'other_type': [pd.np.nan]}
    goal_chron = pd.DataFrame(goal_chrondict)
    df_chron = v.chronology_information.df.copy()
    for k in goal_chrondict.keys():
        pd.testing.assert_series_equal(goal_chron[k], df_chron[k])
