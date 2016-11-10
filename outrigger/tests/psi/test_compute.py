import os

import numpy as np
import pytest
import pandas as pd
import pandas.util.testing as pdt
import six

idx = pd.IndexSlice


@pytest.fixture
def dummy_junction12():
    """Junction between exon_cols 1 and 2"""
    return "junction:chr1:176-224:+"


@pytest.fixture
def dummy_junction23():
    """Junction between exon_cols 2 and 3"""
    return 'junction:chr1:251-299:+'


@pytest.fixture
def dummy_junction13():
    """Junction between exon_cols 1 and 3"""
    return 'junction:chr1:176-299:+'


@pytest.fixture
def dummy_isoform1_junctions(splice_type):
    return ['junction13']


@pytest.fixture
def dummy_isoform2_junctions(splice_type):
    return ['junction12', 'junction23']


@pytest.fixture(params=['isoform1', 'isoform2',
                        pytest.mark.xfail('keyerror_isoform')])
def dummy_isoform_junctions(request, dummy_isoform1_junctions,
                            dummy_isoform2_junctions):
    if request.param == 'isoform1':
        return dummy_isoform1_junctions
    elif request.param == 'isoform2':
        return dummy_isoform2_junctions
    else:
        return 'keyerror_isoform'


@pytest.fixture
def dummy_event_name():
    return 'exon:chr1:150-175:+@exon:chr1:225-250:+@exon:chr1:300-350:+'


@pytest.fixture(params=[100, 2, np.nan],
                ids=['enough reads', 'not enough reads', 'not there'])
def dummy_junction12_reads(request):
    return request.param


@pytest.fixture(params=[100, 2, np.nan],
                ids=['enough reads', 'not enough reads', 'not there'])
def dummy_junction23_reads(request):
    return request.param


@pytest.fixture(params=[100, 2, np.nan],
                ids=['enough reads', 'not enough reads', 'not there'])
def dummy_junction13_reads(request):
    return request.param


@pytest.fixture
def reads_col():
    return 'reads'


@pytest.fixture
def dummy_splice_junction_reads(dummy_junction12, dummy_junction12_reads,
                          dummy_junction23, dummy_junction23_reads,
                          dummy_junction13, dummy_junction13_reads, reads_col):
    """Completely fake dataset for sanity checking"""
    s = """sample_id,junction,{6}
sample1,{0},{1}
sample1,{2},{3}
sample1,{4},{5}""".format(dummy_junction12, dummy_junction12_reads,
                          dummy_junction23, dummy_junction23_reads,
                          dummy_junction13, dummy_junction13_reads, reads_col)
    data = pd.read_csv(six.StringIO(s), comment='#')
    data = data.dropna()
    data = data.set_index(
        ['junction', 'sample_id'])
    data = data.sort_index()
    return data


@pytest.fixture
def dummy_junction_locations(dummy_junction12, dummy_junction23,
                              dummy_junction13):
    return pd.Series({'junction12': dummy_junction12,
                      'junction23': dummy_junction23,
                      'junction13': dummy_junction13,
                      'illegal_junctions': np.nan})


@pytest.fixture
def dummy_junction_to_reads(dummy_junction12, dummy_junction12_reads,
                      dummy_junction23, dummy_junction23_reads,
                      dummy_junction13, dummy_junction13_reads):
    """Helper function for testing"""
    return pd.Series({dummy_junction12: dummy_junction12_reads,
                      dummy_junction23: dummy_junction23_reads,
                      dummy_junction13: dummy_junction13_reads})


def test_maybe_get_isoform_reads(dummy_splice_junction_reads,
                                 dummy_junction_locations,
                                 dummy_isoform_junctions,
                                 dummy_junction_to_reads,
                                 reads_col):
    from outrigger.psi.compute import maybe_get_isoform_reads
    test = maybe_get_isoform_reads(dummy_splice_junction_reads,
                                   dummy_junction_locations,
                                   dummy_isoform_junctions,
                                   reads_col=reads_col)
    junctions = dummy_junction_locations[dummy_isoform_junctions]
    reads = dummy_junction_to_reads[junctions.values]
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
def dummy_events(splice_type):
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


def test_dummy_calculate_psi(dummy_splice_junction_reads,
                             dummy_junction12_reads,
                             dummy_junction23_reads,
                             dummy_junction13_reads,
                             dummy_isoform1_junctions,
                             dummy_isoform2_junctions,
                             dummy_exons_to_junctions, capsys,
                             dummy_events):
    from outrigger.psi.compute import calculate_psi, MIN_READS

    reads12 = dummy_junction12_reads if \
        dummy_junction12_reads >= MIN_READS else 0
    reads23 = dummy_junction23_reads if \
        dummy_junction23_reads >= MIN_READS else 0
    reads13 = dummy_junction13_reads if \
        dummy_junction13_reads >= MIN_READS else 0

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

    test = calculate_psi(dummy_exons_to_junctions, dummy_splice_junction_reads,
                         isoform1_junctions=dummy_isoform1_junctions,
                         isoform2_junctions=dummy_isoform2_junctions,
                         n_jobs=1)
    out, err = capsys.readouterr()
    assert 'Iterating over' in out

    s = """sample_id,{0}# noqa
sample1,{2},{1},{2}""".format(','.join(dummy_events), true_psi,
                              other_isoform1_psi)

    true = pd.read_csv(six.StringIO(s), index_col=0, comment='#')
    true = true.dropna(axis=1)

    if true.empty:
        true = pd.DataFrame(index=dummy_splice_junction_reads.index.levels[1])

    pdt.assert_frame_equal(test, true)


# --- Test with real data --- #
@pytest.fixture
def event_id(splice_type):
    if splice_type == 'se':
        return 'isoform1=junction:chr10:128491034-128492058:-|isoform2=junction:chr10:128491765-128492058:-@novel_exon:chr10:128491720-128491764:-@junction:chr10:128491034-128491719:-'
    elif splice_type == 'mxe':
        return 'isoform1=junction:chr2:136763622-136770056:+@exon:chr2:136770057-136770174:+@junction:chr2:136770175-136773894:+|isoform2=junction:chr2:136763622-136769742:+@exon:chr2:136769743-136769860:+@junction:chr2:136769861-136773894:+' # noqa


@pytest.fixture
def event_df_csv(splice_type, tasic2016_intermediate_psi):
    return os.path.join(tasic2016_intermediate_psi,
                        '{splice_type}_event_df.csv'.format(
                            splice_type=splice_type))


@pytest.fixture
def splice_junction_reads_csv(tasic2016_intermediate_psi):
    return os.path.join(tasic2016_intermediate_psi,
                        'splice_junction_reads.csv')

@pytest.fixture
def splice_junction_reads(splice_junction_reads_csv):
    return pd.read_csv(splice_junction_reads_csv, index_col=[0, 1])


@pytest.fixture
def single_event_psi_csv(tasic2016_intermediate_psi, splice_type):
    return os.path.join(tasic2016_intermediate_psi,
                        '{splice_type}_event_psi.csv'.format(
                            splice_type=splice_type))

@pytest.fixture
def event_annotation_csv(splice_type, tasic2016_outrigger_output_index):
    return os.path.join(tasic2016_outrigger_output_index, splice_type,
                        'events.csv')


@pytest.fixture
def event_annotation(event_annotation_csv):
    return pd.read_csv(event_annotation_csv, index_col=0)


@pytest.fixture
def isoform1_junctions(splice_type):
    if splice_type == 'se':
        return ['junction13']
    if splice_type == 'mxe':
        return ['junction13', 'junction34']


@pytest.fixture
def isoform2_junctions(splice_type):
    if splice_type == 'se':
        return ['junction12', 'junction23']
    if splice_type == 'mxe':
        return ['junction12', 'junction24']


@pytest.fixture
def psi_csv(splice_type, tasic2016_outrigger_output_psi):
    return os.path.join(tasic2016_outrigger_output_psi, splice_type, 'psi.csv')


@pytest.fixture
def psi_df(psi_csv):
    return pd.read_csv(psi_csv, index_col=0)

def test__single_event_psi(event_id, event_df_csv, splice_junction_reads,
                           isoform1_junctions, isoform2_junctions,
                           single_event_psi_csv):
    from outrigger.psi.compute import _single_event_psi
    event_df = pd.read_csv(event_df_csv, index_col=0)

    test = _single_event_psi(event_id, event_df, splice_junction_reads,
                             isoform1_junctions, isoform2_junctions)
    true = pd.read_csv(single_event_psi_csv, index_col=0, squeeze=True)
    pdt.assert_series_equal(test, true)


def test__maybe_parallelize_psi(event_annotation, splice_junction_reads,
                                isoform1_junctions, isoform2_junctions, psi_df,
                                capsys, n_jobs):
    from outrigger.psi.compute import _maybe_parallelize_psi

    test = _maybe_parallelize_psi(event_annotation, splice_junction_reads,
                                  isoform1_junctions, isoform2_junctions,
                                  n_jobs=n_jobs)
    true = [column for name, column in psi_df.iteritems()]

    out, err = capsys.readouterr()

    if n_jobs == 1:
        assert 'Iterating' in out
    else:
        assert 'Parallelizing' in out
    pdt.assert_equal(test, true)


def test_calculate_psi(event_annotation, splice_junction_reads,
                       isoform1_junctions, isoform2_junctions, psi_df):
    from outrigger.psi.compute import calculate_psi

    test = calculate_psi(event_annotation, splice_junction_reads,
                         isoform1_junctions, isoform2_junctions)
    true = psi_df
    pdt.assert_frame_equal(test, true)
