import glob
import os

import pandas as pd
import pytest
import six


@pytest.fixture
def data_folder():
    """Absolute path to where test data is stored"""
    return os.path.join(os.path.abspath(os.path.dirname(__file__)), 'data')


@pytest.fixture
def sj_out_tab(tmpdir):
    s = """chr1    76   299   1       2       1       0       1       39
chr1    201   299   1       1       1       0       1       10
chr1    201  249  1       1       0       0       1       22
chr1    201  799  1       1       1       19      20      43
chr1    201  799  1       1       0       8       15      41
chr1    155832  164262  1       1       1       61      3       46
chr1    156087  156200  1       1       0       1       14      44
chr1    329977  334128  1       1       1       0       2       14
chr1    569184  569583  1       1       0       0       1       17
chr1    655581  659737  1       1       1       0       2       14
chr1    661725  662046  1       1       0       0       1       22
chr1    668587  671992  1       1       0       0       4       28
"""
    df = pd.read_table(six.StringIO(s), header=None, sep='\s+')
    filename = '{0}/SJ.out.tab'.format(tmpdir)
    df.to_csv(filename, index=False, header=False, sep='\t')
    return filename


@pytest.fixture
def sj_filenames(data_folder):
    return glob.glob('{}/*SJ.out.tab'.format(data_folder))


@pytest.fixture
def splice_junctions(sj_filenames):
    from outrigger.star import read_multiple_sj_out_tab

    return read_multiple_sj_out_tab(sj_filenames)


def test_read_sj_out_tab(sj_out_tab):
    pass


def test_int_to_intron_motif():
    pass


def test_read_multiple_sj_out_tab():
    pass


def test_sj_count_to_metadata():
    pass
