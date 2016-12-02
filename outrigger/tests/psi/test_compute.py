import os

import numpy as np
import pytest
import pandas as pd
import pandas.util.testing as pdt
import six

idx = pd.IndexSlice



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


@pytest.fixture
def mutually_exclusive_exon_junction_reads_for_rejecting_csv(simulated):
    return os.path.join(simulated, 'psi', 'mutually_exclusive_exon_junctions_psi.csv')


@pytest.fixture
def mutually_exclusive_exon_junction_reads_for_rejecting(
        mutually_exclusive_exon_junction_reads_for_rejecting_csv):
    df = pd.read_csv(mutually_exclusive_exon_junction_reads_for_rejecting_csv,
                     dtype={'junction12': int, 'junction13': int,
                            'junction23': int, 'psi': float},
                     comment='#', index_col=0)
    return df


@pytest.fixture
def junction_reads_for_rejecting(
    splice_type, skipped_exon_junction_reads_for_rejecting,
    mutually_exclusive_exon_junction_reads_for_rejecting):
    if splice_type == 'se':
        return skipped_exon_junction_reads_for_rejecting
    elif splice_type == 'mxe':
        return mutually_exclusive_exon_junction_reads_for_rejecting


def test__single_isoform_maybe_reject(junction_reads_for_rejecting,
                                      dummy_isoform1_junction_numbers,
                                      dummy_isoform2_junction_numbers):
    from outrigger.psi.compute import _single_isoform_maybe_reject

    n_junctions = len(dummy_isoform1_junction_numbers) \
                  + len(dummy_isoform2_junction_numbers)

    for i, row in junction_reads_for_rejecting.iterrows():
        # debug = True if row.name == 'junction13 ≥ 10' else False
        debug = False
        isoform1, isoform2, case = _single_isoform_maybe_reject(
            row[dummy_isoform1_junction_numbers],
            row[dummy_isoform2_junction_numbers],
            n_junctions=n_junctions, debug=debug)
        try:
            if np.isnan(row['mean_psi']):
                assert isoform1 is None
                assert isoform2 is None
            else:
                pdt.assert_series_equal(isoform1,
                                        row[dummy_isoform1_junction_numbers])
                pdt.assert_series_equal(isoform2,
                                        row[dummy_isoform2_junction_numbers])
            print(row.name, 'passed')
        except AssertionError:
            raise AssertionError('The junction configuration [{title}] did '
                                 'not pass (Test: {test_case}, '
                                 'True: {true_case})'.format(
                title=row.name, test_case=case, true_case=row['case']))


@pytest.fixture(params=[({'junction12': 1000, 'junction23': 20}, 'unequal'),
                        ({'junction12': 100, 'junction23': 20}, 'similar'),
                        ({'junction12': 20, 'junction23': 1000}, 'unequal'),
                        ({'junction23': 1000}, 'one junction'),
                        ({'junction12': 20, 'junction23': 100}, 'similar')])
def isoform_read_coverage(request):
    return pd.Series(request.param[0]), request.param[1]


def test__check_unequal_read_coverage(isoform_read_coverage):
    from outrigger.psi.compute import _single_sample_check_unequal_read_coverage

    isoform, equality = isoform_read_coverage

    test = _single_sample_check_unequal_read_coverage(isoform)
    if equality == 'unequal':
        assert test is None
    else:
        pdt.assert_series_equal(isoform, test)


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
def incompatible_junctions(splice_type, dummy_junction14, dummy_junction23):
    if splice_type == 'se':
        return np.nan
    if splice_type == 'mxe':
        return dummy_junction14 + '|' + dummy_junction23


@pytest.fixture
def dummy_junction_locations(dummy_legal_junction_numbers,
                             dummy_junction_number_to_id,
                             incompatible_junctions):
    from outrigger.common import INCOMPATIBLE_JUNCTIONS

    d = {junction_xy: dummy_junction_number_to_id[junction_xy]
         for junction_xy in dummy_legal_junction_numbers}
    d[INCOMPATIBLE_JUNCTIONS] = incompatible_junctions
    return pd.Series(d)


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




# --- Test with real data --- #
@pytest.fixture
def event_id(splice_type):
    if splice_type == 'se':
        return 'isoform1=junction:chr10:128491034-128492058:-|isoform2=junction:chr10:128491765-128492058:-@novel_exon:chr10:128491720-128491764:-@junction:chr10:128491034-128491719:-'  # noqa
    elif splice_type == 'mxe':
        return 'isoform1=junction:chr2:136763622-136770056:+@exon:chr2:136770057-136770174:+@junction:chr2:136770175-136773894:+|isoform2=junction:chr2:136763622-136769742:+@exon:chr2:136769743-136769860:+@junction:chr2:136769861-136773894:+'  # noqa


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
def reads2d(splice_junction_reads_csv):
    from outrigger.common import SAMPLE_ID, JUNCTION_ID, READS
    df = pd.read_csv(splice_junction_reads_csv)
    df2d = df.pivot(index=SAMPLE_ID, columns=JUNCTION_ID, values=READS)
    df2d = df2d.fillna(0)
    return df2d


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
def event_df_csv(splice_type, tasic2016_intermediate_psi):
    return os.path.join(tasic2016_intermediate_psi,
                        '{splice_type}_event_df.csv'.format(
                            splice_type=splice_type))

@pytest.fixture
def single_event_summary_csv(splice_type, tasic2016_intermediate_psi):
    return os.path.join(tasic2016_intermediate_psi, splice_type, 'summary.csv')


def test__single_event_psi(event_id, event_df_csv, reads2d,
                           isoform1_junctions, isoform2_junctions,
                           single_event_summary_csv):
    from outrigger.psi.compute import _single_event_psi
    event_df = pd.read_csv(event_df_csv, index_col=0)

    true = pd.read_csv(single_event_summary_csv)

    test = _single_event_psi(event_id, event_df, reads2d,
                             isoform1_junctions, isoform2_junctions)
    pdt.assert_frame_equal(test, true)


@pytest.fixture
def psi_csv(splice_type, tasic2016_outrigger_output_psi):
    return os.path.join(tasic2016_outrigger_output_psi, splice_type, 'psi.csv')


@pytest.fixture
def psi_df(psi_csv):
    return pd.read_csv(psi_csv, index_col=0)


@pytest.fixture
def summary_csv(splice_type, tasic2016_outrigger_output_psi):
    return os.path.join(tasic2016_outrigger_output_psi, splice_type,
                        'summary.csv')


@pytest.fixture
def summary_df(summary_csv):
    return pd.read_csv(summary_csv)


def test__maybe_parallelize_psi(event_annotation, reads2d,
                                isoform1_junctions, isoform2_junctions,
                                capsys, n_jobs, summary_df):
    from outrigger.psi.compute import _maybe_parallelize_psi

    tests = _maybe_parallelize_psi(event_annotation, reads2d,
                                   isoform1_junctions, isoform2_junctions,
                                   n_jobs=n_jobs)
    tests = [t for t in tests if t is not None]
    trues = [df for name, df in summary_df.groupby('event_id')]

    out, err = capsys.readouterr()

    if n_jobs == 1:
        assert 'Iterating' in out
    else:
        assert 'Parallelizing' in out

    for test, true in zip(tests, trues):
        # Only need to make sure the values are the same, the index
        # doesn't matter
        true.index = test.index

        # Make sure columns are in same order
        true = true[test.columns]
        pdt.assert_frame_equal(test, true)


def test_calculate_psi(event_annotation, reads2d,
                       isoform1_junctions, isoform2_junctions,
                       psi_df, summary_df):
    from outrigger.psi.compute import calculate_psi

    true_psi = psi_df
    true_summary = summary_df

    test_psi, test_summary = calculate_psi(event_annotation, reads2d,
                                           isoform1_junctions,
                                           isoform2_junctions)
    # When psi is written to CSV, only the index name is preserved, not the
    # column names so need to get rid of it for these comparisons
    test_psi.columns.name = None
    pdt.assert_frame_equal(test_psi, true_psi)
    pdt.assert_frame_equal(test_summary, true_summary)
