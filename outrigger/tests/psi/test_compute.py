import os

import numpy as np
import pytest
import pandas as pd
import pandas.util.testing as pdt
import six

idx = pd.IndexSlice


@pytest.fixture
def junction12():
    """Junction between exon_cols 1 and 2"""
    return "junction:chr1:176-224:+"


@pytest.fixture
def junction23():
    """Junction between exon_cols 2 and 3"""
    return 'junction:chr1:251-299:+'


@pytest.fixture
def junction13():
    """Junction between exon_cols 1 and 3"""
    return 'junction:chr1:176-299:+'


@pytest.fixture
def isoform1_junctions():
    return ['junction13']


@pytest.fixture
def isoform2_junctions():
    return ['junction12', 'junction23']


@pytest.fixture(params=['isoform1', 'isoform2',
                        pytest.mark.xfail('keyerror_isoform')])
def isoform_junctions(request, isoform1_junctions, isoform2_junctions):
    if request.param == 'isoform1':
        return isoform1_junctions
    elif request.param == 'isoform2':
        return isoform2_junctions
    else:
        return 'keyerror_isoform'


@pytest.fixture
def event_name():
    return 'exon:chr1:150-175:+@exon:chr1:225-250:+@exon:chr1:300-350:+'


@pytest.fixture(params=[100, 2, np.nan],
                ids=['enough reads', 'not enough reads', 'not there'])
def junction12_reads(request):
    return request.param


@pytest.fixture(params=[100, 2, np.nan],
                ids=['enough reads', 'not enough reads', 'not there'])
def junction23_reads(request):
    return request.param


@pytest.fixture(params=[100, 2, np.nan],
                ids=['enough reads', 'not enough reads', 'not there'])
def junction13_reads(request):
    return request.param


@pytest.fixture
def reads_col():
    return 'reads'


@pytest.fixture
def splice_junction_reads(junction12, junction12_reads,
                          junction23, junction23_reads,
                          junction13, junction13_reads, reads_col):
    s = """sample_id,junction,{6}
sample1,{0},{1}
sample1,{2},{3}
sample1,{4},{5}""".format(junction12, junction12_reads,
                          junction23, junction23_reads,
                          junction13, junction13_reads, reads_col)
    data = pd.read_csv(six.StringIO(s), comment='#')
    data = data.dropna()
    data = data.set_index(
        ['junction', 'sample_id'])
    data = data.sort_index()
    return data


@pytest.fixture
def junction_locations(junction12, junction23, junction13):
    return pd.Series({'junction12': junction12, 'junction23': junction23,
                      'junction13': junction13, 'illegal_junctions': np.nan})


@pytest.fixture
def junction_to_reads(junction12, junction12_reads,
                      junction23, junction23_reads,
                      junction13, junction13_reads):
    """Helper function for testing"""
    return pd.Series({junction12: junction12_reads,
                      junction23: junction23_reads,
                      junction13: junction13_reads})


def test_maybe_get_isoform_reads(splice_junction_reads, junction_locations,
                                 isoform_junctions, junction_to_reads,
                                 reads_col):
    from outrigger.psi.compute import maybe_get_isoform_reads
    test = maybe_get_isoform_reads(splice_junction_reads, junction_locations,
                                   isoform_junctions, reads_col=reads_col)
    junctions = junction_locations[isoform_junctions]
    reads = junction_to_reads[junctions.values]
    reads = reads.dropna()
    if reads.empty:
        true = pd.Series()
    else:
        true = reads.copy()
        true.name = reads_col
        true.index = pd.MultiIndex.from_arrays(
            [reads.index, ['sample1']*len(reads.index)],
            names=['junction', 'sample_id'])
    pdt.assert_series_equal(test, true)


@pytest.fixture
def exons_to_junctions(splice_type, simulated_outrigger_index):
    # strand_str = 'positive' if strand == "+" else 'negative'

    folder = os.path.join(simulated_outrigger_index, splice_type)
    basename = 'events_positive_strand.csv'
    filename = os.path.join(folder, basename)
    return pd.read_csv(filename, index_col=0)


@pytest.fixture
def events(splice_type):
    if splice_type == 'se':
        # if strand == '+':
        return ['isoform1=junction:chr1:176-299:+|isoform2=junction:chr1:176-199:+@exon:chr1:200-250:+@junction:chr1:251-299:+',  # noqa
                'isoform1=junction:chr1:176-299:+|isoform2=junction:chr1:176-224:+@exon:chr1:225-250:+@junction:chr1:251-299:+',  # noqa
                'isoform1=junction:chr1:176-299:+|isoform2=junction:chr1:176-224:+@exon:chr1:225-275:+@junction:chr1:276-299:+']  # noqa
    if splice_type == 'mxe':
        # if strand == '+':
        return ['isoform1=junction:chr1:176-299:+@exon:chr1:300-350:+@junction:chr1:351-399:+|isoform2=junction:chr1:176-224:+@exon:chr1:225-250:+@junction:chr1:251-399:+']  # noqa


@pytest.fixture
def illegal_junctions(splice_type):
    if splice_type == 'se':
        return None
    if splice_type == 'mxe':
        return 'junction23'


def test__single_event_psi(splice_junction_reads, isoform1_junctions,
                           isoform2_junctions):
    from outrigger.psi.compute import _single_event_psi

    test = _single_event_psi(event_id, event_df, splice_junction_reads,
                             isoform1_junctions, isoform2_junctions)



def test_calculate_psi():
    from outrigger.psi.compute import calculate_psi


def test_psi(splice_junction_reads, junction12_reads, junction23_reads,
             junction13_reads, exons_to_junctions, capsys, events):
    from outrigger.psi.compute import calculate_psi, MIN_READS

    reads12 = junction12_reads if junction12_reads >= MIN_READS else 0
    reads23 = junction23_reads if junction23_reads >= MIN_READS else 0
    reads13 = junction13_reads if junction13_reads >= MIN_READS else 0

    if reads12 == 0 or reads23 == 0:
        isoform2_reads = 0
    else:
        isoform2_reads = reads12 + reads23
    isoform1_reads = reads13

    # This tests whether both are greater than zero
    if isoform1_reads or isoform2_reads:
        true_psi = isoform2_reads/(isoform2_reads + 2.*isoform1_reads)
    else:
        true_psi = np.nan

    other_isoform1_psi = 0. if isoform1_reads > 0 else np.nan

    test = calculate_psi(exons_to_junctions, splice_junction_reads,
                         isoform1_junctions=['junction13'],
                         isoform2_junctions=['junction12', 'junction23'])
    out, err = capsys.readouterr()
    assert 'Iterating over' in out

    s = """sample_id,{0}# noqa
sample1,{2},{1},{2}""".format(','.join(events), true_psi, other_isoform1_psi)

    true = pd.read_csv(six.StringIO(s), index_col=0, comment='#')
    true = true.dropna(axis=1)

    if true.empty:
        true = pd.DataFrame(index=splice_junction_reads.index.levels[1])

    pdt.assert_frame_equal(test, true)
