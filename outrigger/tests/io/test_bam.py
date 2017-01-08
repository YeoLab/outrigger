import collections
import glob
import os

import pandas as pd
import pandas.util.testing as pdt
import pysam
import pytest


@pytest.fixture
def bamfile(tasic2016_bam):
    return os.path.join(
        tasic2016_bam,
        'CAV_LP_Ipsi_tdTpos_cell_1-SRR2140356-GSM1840944_R1.polyATrim.'
        'adapterTrim.rmRep.sorted.rg.subset.sorted.bam')


@pytest.fixture
def single_bam_final_junction_reads_table_csv(bamfile,
                                              tasic2016_intermediate_bam,
                                              ignore_multimapping, stranded):
    """csv for bamfile (singular)"""
    basename = os.path.basename(bamfile)
    strandedness = '' if stranded else '_unstranded'
    suffix = '.outrigger_junction_reads_ignore-multimapping{}{}.csv'.format(
        ignore_multimapping, strandedness)
    basename2 = basename + suffix
    csv = os.path.join(tasic2016_intermediate_bam, basename2)
    return csv


@pytest.fixture
def suffix_template(bamfile):
    return os.path.basename(bamfile).replace(
        '.bam', '.junction_reads_{}mapped{}.csv')


@pytest.fixture
def uniquely_csv(suffix_template, tasic2016_intermediate_bam, stranded):
    strandedness = '' if stranded else '_unstranded'
    return os.path.join(tasic2016_intermediate_bam,
                        suffix_template.format('uniquely', strandedness))


@pytest.fixture
def multi_csv(suffix_template, tasic2016_intermediate_bam, stranded):
    strandedness = '' if stranded else '_unstranded'
    return os.path.join(tasic2016_intermediate_bam,
                        suffix_template.format('multi', strandedness))


def read_intermediate_junctions(csv, stranded):
    index_col = [0, 1, 2, 3] if stranded else [0, 1, 2]
    return pd.read_csv(csv, index_col=index_col,
                       squeeze=True, header=None)


@pytest.fixture
def multi(multi_csv, stranded):
    return read_intermediate_junctions(multi_csv, stranded).to_dict()


@pytest.fixture
def uniquely(uniquely_csv, stranded):
    return read_intermediate_junctions(uniquely_csv, stranded).to_dict()


@pytest.fixture
def uniquely_summed_csv(bamfile, tasic2016_intermediate_bam):
    basename = os.path.basename(bamfile)
    basename = basename.replace('.bam', '.uniquely_mapped_summed.csv')

    csv = os.path.join(tasic2016_intermediate_bam, basename)
    return csv


@pytest.fixture(params=[None, 'multi', 'uniquely'])
def empty(request):
    return request.param


@pytest.fixture
def single_bam_combined_uniquely_multi_csv(tasic2016_intermediate_bam, bamfile,
                                           ignore_multimapping, empty,
                                           stranded):
    basename = os.path.basename(bamfile)
    empty_suffix = '' if empty is None else '_empty_' + empty
    stranded_suffix = '' if stranded else '_unstranded'
    basename = basename.replace(
        '.bam', '.junction_reads_ignore-multimapping{}{}{}.csv')
    basename = basename.format(ignore_multimapping, empty_suffix,
                               stranded_suffix)
    return os.path.join(tasic2016_intermediate_bam, basename)


@pytest.fixture
def multiple_bams_reads_table_csvs(tasic2016_intermediate_bam,
                                   ignore_multimapping, stranded):
    """csvs for bamfiles (multiple)"""
    strandedness = '' if stranded else '_unstranded'
    suffix = '*.outrigger_junction_reads_ignore-multimapping{}{}.csv'.format(
        ignore_multimapping, strandedness)
    csvs = glob.glob(os.path.join(tasic2016_intermediate_bam, suffix))
    return csvs


def test__report_read_positions(bamfile, stranded):
    from outrigger.io.bam import _report_read_positions

    bam = pysam.AlignmentFile(bamfile, 'rb')

    test = collections.Counter()

    for read in bam:
        _report_read_positions(read, test, stranded)
        break
    bam.close()

    if stranded:
        true = {('chr2', 136713559, 136713559, '+'): 1}
    else:
        true = {('chr2', 136713559, 136713559): 1}

    pdt.assert_dict_equal(test, true)


def test__combine_uniquely_multi(uniquely, multi, ignore_multimapping, empty,
                                 single_bam_combined_uniquely_multi_csv,
                                 stranded):
    from outrigger.io.bam import _combine_uniquely_multi

    u = uniquely
    m = multi
    if empty == 'uniquely':
        u = {}
    elif empty == 'multi':
        m = {}

    test = _combine_uniquely_multi(u, m, ignore_multimapping, stranded)
    true = pd.read_csv(single_bam_combined_uniquely_multi_csv)

    pdt.assert_frame_equal(test, true)


def test__get_junction_reads(bamfile, uniquely, multi, stranded):
    from outrigger.io.bam import _get_junction_reads

    test_uniquely, test_multi = _get_junction_reads(bamfile, stranded)

    true_uniquely = uniquely
    true_multi = multi

    pdt.assert_dict_equal(test_uniquely, true_uniquely)
    pdt.assert_dict_equal(test_multi, true_multi)


def test_bam_to_junction_reads_table(
        bamfile, single_bam_final_junction_reads_table_csv,
        ignore_multimapping, stranded):

    from outrigger.io.bam import bam_to_junction_reads_table

    test = bam_to_junction_reads_table(bamfile, ignore_multimapping, stranded)
    true = pd.read_csv(single_bam_final_junction_reads_table_csv)

    pdt.assert_frame_equal(test, true)


def test_read_multiple_bams(bam_filenames, multiple_bams_reads_table_csvs,
                            ignore_multimapping, stranded):
    from outrigger.io.bam import read_multiple_bams

    test = read_multiple_bams(bam_filenames, ignore_multimapping, stranded)

    dfs = [pd.read_csv(csv) for csv in multiple_bams_reads_table_csvs]
    true = pd.concat(dfs, ignore_index=True)

    # Sort and change the index because it's the contents not the order that
    # matters
    test = test.sort_values(test.columns.tolist())
    test.index = range(len(test.index))
    true = true.sort_values(true.columns.tolist())
    true.index = range(len(true.index))

    pdt.assert_frame_equal(test, true)
