import os

import pandas as pd
import pandas.util.testing as pdt
import pytest
import six


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


def test_read_sj_out_tab(sj_out_tab, simulated_unprocessed):
    from outrigger.io.star import read_sj_out_tab

    test = read_sj_out_tab(sj_out_tab)
    csv = os.path.join(simulated_unprocessed, 'true_splice_junctions.csv')
    true = pd.read_csv(csv)
    assert (test.junction_start < test.junction_stop).all()
    pdt.assert_frame_equal(test, true)


def test_int_to_intron_motif():
    from outrigger.io.star import int_to_junction_motif

    ints = [0, 1, 2, 3, 4, 5, 6]
    test = [int_to_junction_motif(i) for i in ints]
    true = ['non-canonical', 'GT/AG', 'GT/AG', 'GC/AG', 'GC/AG', 'AT/AC',
            'AT/AC']
    assert test == true


@pytest.fixture(params=[True, False])
def ignore_multimapping(request):
    return request.param


@pytest.fixture
def splice_junction_csv(ignore_multimapping, tasic2016_intermediate):
    """Different file depending on whether multimapping is True"""
    template = os.path.join(tasic2016_intermediate,
                            'index', 'star',
                            'splice_junctions_ignore_multimapping{}.csv')
    return template.format(str(ignore_multimapping))


def test_read_multiple_sj_out_tab(sj_filenames, ignore_multimapping,
                                  splice_junction_csv):
    from outrigger.io.star import read_multiple_sj_out_tab
    from outrigger.common import READS

    # Read csv file and convert to numeric
    true = pd.read_csv(splice_junction_csv)
    true = true.convert_objects()

    test = read_multiple_sj_out_tab(
        sj_filenames, ignore_multimapping=ignore_multimapping)
    assert READS in test
    pdt.assert_frame_equal(test, true)


def test_make_metadata(junction_metadata, junction_reads):
    from outrigger.io.star import make_metadata

    true = junction_metadata
    test = make_metadata(junction_reads)
    pdt.assert_frame_equal(test, true)
