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
               '# Contribution_Date',
               '#',
               '#------------------------',
               '# Title',
               '#',
               '#------------------------',
               '# Chronology_Information',
               '#',
               '# Labcode\tdepth_top\tdepth_bottom\tmat_dated\t14C_date\t14C_1s_err\tdelta_R\tdelta_R_1s_err\tother_date\tother_1s_err\tother_type\t',
               '#',
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
    goal_beginline = 15
    assert dumb_guts.data == goal_data
    assert dumb_guts.description[0] == goal_description_ends[0]
    assert dumb_guts.description[-1] == goal_description_ends[-1]
    assert dumb_guts.data_beginline == goal_beginline


def test__write_sectionindex(dumb_guts):
    goal_index = {'Chronology_Information': (7, 11),
                  'Contribution_Date': (1, 3),
                  'Data': (12, 15),
                  'Title': (4, 6)}
    assert len(dumb_guts.sectionindex.keys()) == 4
    for k, goal in goal_index.items():
        assert dumb_guts.sectionindex[k] == goal


def test_pull_section(dumb_guts):
    goal = ['# Data', '# Data lines follow (have no #)', '# Missing Value: -999']
    assert dumb_guts.pull_section('Data') == goal


def test_pull_section_badname(dumb_guts):
    with pytest.raises(KeyError):
        dumb_guts.pull_section('data')


def test_available_sections(dumb_guts):
    goal = ['Contribution_Date', 'Title', 'Chronology_Information', 'Data']
    assert dumb_guts.available_sections() == goal


def test_guess_missingvalues(dumb_guts):
    goal = '-999'
    assert dumb_guts.guess_missingvalues() == goal


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
