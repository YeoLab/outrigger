import os

import numpy as np
import pytest
import pandas as pd
import pandas.util.testing as pdt
import six

idx = pd.IndexSlice



@pytest.fixture(params=({'isoform1': 99, 'isoform2': 99, 'result': (99, 99)},
                       {'isoform1': 2, 'isoform2': 99, 'result': (2, 99)},
                       {'isoform1': 99, 'isoform2': 2, 'result': (99, 2)},
                       {'isoform1': 99, 'isoform2': 0, 'result': (99, 0)},
                       {'isoform1': 2, 'isoform2': 26, 'result': (np.nan,
                                                                  np.nan)},
                       {'isoform1': 0, 'isoform2': 15, 'result': (0, 15)},
                       {'isoform1': 0, 'isoform2': 99, 'result': (0, 99)},
                       {'isoform1': 2, 'isoform2': 2, 'result':
                           (np.nan, np.nan)}))
def maybe_sufficient_isoforms(request):
    index = ['sample1']

    isoform1 = pd.Series(request.param['isoform1'], index=index)
    isoform2 = pd.Series(request.param['isoform2'], index=index)

    expected1 = pd.Series(request.param['result'][0], index=index).dropna()
    expected2 = pd.Series(request.param['result'][1], index=index).dropna()

    expected1 = expected1.astype(int)
    expected2 = expected2.astype(int)
    return isoform1, isoform2, expected1, expected2


def test__remove_insufficient_reads(maybe_sufficient_isoforms):
    from outrigger.psi.compute import _remove_insufficient_reads

    isoform1, isoform2, expected1, expected2 = maybe_sufficient_isoforms

    test1, test2 = _remove_insufficient_reads(isoform1, isoform2)
    pdt.assert_series_equal(test1, expected1)
    pdt.assert_series_equal(test2, expected2)


@pytest.fixture
def skipped_exon_junction_reads_for_rejecting_csv(simulated):
    return os.path.join(simulated, 'psi', 'skipped_exon_junctions_psi.csv')


@pytest.fixture
def skipped_exon_junction_reads_for_rejecting(
        skipped_exon_junction_reads_for_rejecting_csv):
    df = pd.read_csv(skipped_exon_junction_reads_for_rejecting_csv,
                       dtype={'junction12': int, 'junction13': int,
                              'junction23': int, 'psi': float},
                       comment='#', index_col=0)
    return df


def test__maybe_reject_events(skipped_exon_junction_reads_for_rejecting):
    from outrigger.psi.compute import _maybe_reject_events

    for i, row in skipped_exon_junction_reads_for_rejecting.iterrows():
        isoform1, isoform2 = _maybe_reject_events(
            row[['junction13']],
            row[['junction12',
                                                       'junction23']],
                                            n_junctions=3)
        if np.isnan(row['psi']):
            assert isoform1 is None
            assert isoform2 is None
        else:
            pdt.assert_series_equal(
                isoform1,
                row[['junction13']])
            pdt.assert_series_equal(
                isoform2,
                row[
                    ['junction12', 'junction23']])


@pytest.fixture
def dummy_junction12():
    """Junction between exons 1 and 2"""
    return "junction:chr1:176-224:+"


@pytest.fixture
def dummy_junction23():
    """Junction between exons 2 and 3"""
    return 'junction:chr1:251-299:+'


@pytest.fixture
def dummy_junction13():
    """Junction between exons 1 and 3"""
    return 'junction:chr1:176-299:+'


@pytest.fixture
def dummy_junction14():
    """Junction between exons 1 and 4"""
    return "junction:chr1:176-324:+"


@pytest.fixture
def dummy_junction24():
    """Junction between exons 2 and 4"""
    return 'junction:chr1:251-399:+'


@pytest.fixture
def dummy_junction34():
    """Junction between exons 3 and 4"""
    return 'junction:chr1:351-399:+'


@pytest.fixture
def dummy_junction_number_to_id(dummy_junction12, dummy_junction13,
                                dummy_junction23, dummy_junction14,
                                dummy_junction24, dummy_junction34):
    d = {'junction12': dummy_junction12, 'junction13': dummy_junction13,
         'junction23': dummy_junction23, 'junction14': dummy_junction14,
         'junction24': dummy_junction24, 'junction34': dummy_junction34}
    return d


@pytest.fixture
def dummy_isoform1_junction_numbers(splice_type):
    if splice_type == 'se':
        return ['junction13']
    if splice_type == 'mxe':
        return ['junction13', 'junction34']


@pytest.fixture
def dummy_isoform2_junction_numbers(splice_type):
    if splice_type == 'se':
        return ['junction12', 'junction23']
    if splice_type == 'mxe':
        return ['junction12', 'junction24']


@pytest.fixture
def dummy_isoform1_junction_ids(dummy_isoform1_junction_numbers,
                                dummy_junction_number_to_id):
    ids = [dummy_junction_number_to_id[j]
           for j in dummy_isoform1_junction_numbers]
    return ids


@pytest.fixture
def dummy_isoform2_junction_ids(dummy_isoform2_junction_numbers,
                                dummy_junction_number_to_id):
    ids = [dummy_junction_number_to_id[j]
           for j in dummy_isoform2_junction_numbers]
    return ids


@pytest.fixture
def dummy_legal_junction_numbers(dummy_isoform1_junction_numbers,
                                 dummy_isoform2_junction_numbers):
    return dummy_isoform1_junction_numbers + dummy_isoform2_junction_numbers


@pytest.fixture(params=[(100, 100), (2, 2),
                        (np.nan, np.nan), (100, 2),
                        (2, np.nan), (100, np.nan)],
                ids=['enough reads', 'not enough reads', 'not there',
                     'one not enough', 'not enough, one not there',
                     'one not there'])
def dummy_isoform1_reads(request, dummy_isoform1_junction_ids):
    reads = dict(zip(dummy_isoform1_junction_ids, request.param))
    return reads


@pytest.fixture(params=[(100, 100), (2, 2),
                        (np.nan, np.nan), (100, 2),
                        (2, np.nan), (100, np.nan)],
                ids=['enough reads', 'not enough reads', 'not there',
                     'one not enough', 'not enough, one not there',
                     'one not there'])
def dummy_isoform2_reads(request, dummy_isoform2_junction_ids,):
    reads = dict(zip(dummy_isoform2_junction_ids, request.param))
    return reads


@pytest.fixture
def dummy_isoform_reads(dummy_isoform1_reads, dummy_isoform2_reads):
    """Combine isoform1 and isoform2 reads into one dict"""
    reads = dummy_isoform1_reads.copy()
    reads.update(dummy_isoform2_reads)
    reads = pd.Series(reads)
    return reads


@pytest.fixture(params=['isoform1', 'isoform2',
                        pytest.mark.xfail('keyerror_isoform')])
def dummy_isoform_junctions(request, dummy_isoform1_junction_numbers,
                            dummy_isoform2_junction_numbers):
    if request.param == 'isoform1':
        return dummy_isoform1_junction_numbers
    elif request.param == 'isoform2':
        return dummy_isoform2_junction_numbers
    else:
        return 'keyerror_isoform'


@pytest.fixture
def dummy_splice_junction_reads(dummy_isoform_reads):
    """Completely fake dataset for sanity checking"""
    from outrigger.common import READS

    s = 'sample_id,junction,{reads}\n'.format(reads=READS)
    for junction_id, reads in dummy_isoform_reads.iteritems():
        s += 'sample1,{junction_id},{reads}\n'.format(junction_id=junction_id,
                                                      reads=reads)
    data = pd.read_csv(six.StringIO(s), comment='#')
    data = data.dropna()
    data = data.set_index(
        ['junction', 'sample_id'])
    data = data.sort_index()
    return data


@pytest.fixture
def illegal_junctions(splice_type, dummy_junction14, dummy_junction23):
    if splice_type == 'se':
        return np.nan
    if splice_type == 'mxe':
        return dummy_junction14 + '|' + dummy_junction23


@pytest.fixture
def dummy_junction_locations(dummy_legal_junction_numbers,
                             dummy_junction_number_to_id,
                             illegal_junctions):
    from outrigger.common import ILLEGAL_JUNCTIONS

    d = {junction_xy: dummy_junction_number_to_id[junction_xy]
         for junction_xy in dummy_legal_junction_numbers}
    d[ILLEGAL_JUNCTIONS] = illegal_junctions
    return pd.Series(d)

# @pytest.fixture
# def dummy_junction_to_reads(dummy_junction12, dummy_junction12_reads,
#                       dummy_junction23, dummy_junction23_reads,
#                       dummy_junction13, dummy_junction13_reads):
#     """Helper function for testing"""
#     return pd.Series({dummy_junction12: dummy_junction12_reads,
#                       dummy_junction23: dummy_junction23_reads,
#                       dummy_junction13: dummy_junction13_reads})


def test_maybe_get_isoform_reads(dummy_splice_junction_reads,
                                 dummy_junction_locations,
                                 dummy_isoform_junctions,
                                 dummy_isoform_reads,):
    from outrigger.psi.compute import _maybe_get_isoform_reads, READS

    test = _maybe_get_isoform_reads(dummy_splice_junction_reads,
                                    dummy_junction_locations,
                                    dummy_isoform_junctions, READS)
    junctions = dummy_junction_locations[dummy_isoform_junctions]
    reads = dummy_isoform_reads[junctions.values]
    reads = reads.dropna()
    if reads.empty:
        true = pd.Series()
    else:
        true = reads.copy()
        true.name = READS
        true.index = pd.MultiIndex.from_arrays(
            [reads.index, ['sample1']*len(reads.index)],
            names=['junction', 'sample_id'])
    pdt.assert_series_equal(test, true)


@pytest.fixture
def dummy_exons_to_junctions(splice_type, simulated_outrigger_index):
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


def test_dummy_calculate_psi(dummy_splice_junction_reads,
                             dummy_isoform_reads,
                             dummy_isoform1_junction_ids,
                             dummy_isoform2_junction_ids,
                             dummy_isoform1_junction_numbers,
                             dummy_isoform2_junction_numbers,
                             dummy_exons_to_junctions,
                             dummy_events, splice_type):
    from outrigger.psi.compute import calculate_psi, MIN_READS

    isoform_reads = dummy_isoform_reads.copy()
    # isoform_reads[isoform_reads < MIN_READS] = np.nan


    isoform1_reads = isoform_reads[dummy_isoform1_junction_ids].sum()

    isoform2_reads = isoform_reads[dummy_isoform2_junction_ids].sum()

    # This tests whether both are greater than zero
    if isoform1_reads or isoform2_reads:
        multiplier = float(len(dummy_isoform2_junction_ids)) / \
                     len(dummy_isoform1_junction_ids)
        true_psi = isoform2_reads/(isoform2_reads +
                                   multiplier * isoform1_reads)
    else:
        true_psi = np.nan

    other_isoform1_psi = 0. if isoform1_reads > 0 else np.nan

    if splice_type == 'se':
        s = """sample_id,{0}# noqa
sample1,{2},{1},{2}""".format(','.join(dummy_events), true_psi,
                              other_isoform1_psi)
    if splice_type == 'mxe':
        s = """sample_id,{0}# noqa
sample1,{1}""".format(','.join(dummy_events), true_psi)

    true = pd.read_csv(six.StringIO(s), index_col=0, comment='#')
    true = true.dropna(axis=1)

    if true.empty:
        true = pd.DataFrame(index=dummy_splice_junction_reads.index.levels[1])

    test = calculate_psi(dummy_exons_to_junctions, dummy_splice_junction_reads,
                         isoform1_junctions=dummy_isoform1_junction_numbers,
                         isoform2_junctions=dummy_isoform2_junction_numbers,
                         n_jobs=1)

    pdt.assert_frame_equal(test, true)



# --- Test with real data --- #
@pytest.fixture
def event_id(splice_type):
    if splice_type == 'se':
        return 'isoform1=junction:chr10:128491034-128492058:-|isoform2=junction:chr10:128491765-128492058:-@novel_exon:chr10:128491720-128491764:-@junction:chr10:128491034-128491719:-'  # noqa
    elif splice_type == 'mxe':
        return 'isoform1=junction:chr2:136763622-136770056:+@exon:chr2:136770057-136770174:+@junction:chr2:136770175-136773894:+|isoform2=junction:chr2:136763622-136769742:+@exon:chr2:136769743-136769860:+@junction:chr2:136769861-136773894:+'  # noqa


@pytest.fixture
def event_df_csv(splice_type, tasic2016_intermediate_psi):
    return os.path.join(tasic2016_intermediate_psi,
                        '{splice_type}_event_df.csv'.format(
                            splice_type=splice_type))


@pytest.fixture
def splice_junction_reads_csv(tasic2016_outrigger_junctions):
    return os.path.join(tasic2016_outrigger_junctions, 'reads.csv')


@pytest.fixture
def splice_junction_reads(splice_junction_reads_csv):
    df = pd.read_csv(splice_junction_reads_csv,
                     index_col=['junction_id', 'sample_id'])
    df = df.sort_index()
    return df


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
                           isoform1_junctions, isoform2_junctions, psi_df):
    from outrigger.psi.compute import _single_event_psi
    event_df = pd.read_csv(event_df_csv, index_col=0)

    test = _single_event_psi(event_id, event_df, splice_junction_reads,
                             isoform1_junctions, isoform2_junctions)
    true = psi_df[test.name]
    pdt.assert_series_equal(test, true)


def test__maybe_parallelize_psi(event_annotation, splice_junction_reads,
                                isoform1_junctions, isoform2_junctions, psi_df,
                                capsys, n_jobs):
    from outrigger.psi.compute import _maybe_parallelize_psi

    tests = _maybe_parallelize_psi(event_annotation, splice_junction_reads,
                                   isoform1_junctions, isoform2_junctions,
                                   n_jobs=n_jobs)
    tests = [t for t in tests if t is not None]
    trues = [column.dropna() for name, column in psi_df.iteritems()]

    out, err = capsys.readouterr()

    if n_jobs == 1:
        assert 'Iterating' in out
    else:
        assert 'Parallelizing' in out

    for test, true in zip(tests, trues):
        pdt.assert_series_equal(test, true)


def test_calculate_psi(event_annotation, splice_junction_reads,
                       isoform1_junctions, isoform2_junctions, psi_df):
    from outrigger.psi.compute import calculate_psi

    test = calculate_psi(event_annotation, splice_junction_reads,
                         isoform1_junctions, isoform2_junctions)
    true = psi_df
    pdt.assert_frame_equal(test, true)
