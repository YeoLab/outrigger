import os
import glob

import pandas as pd
import pandas.util.testing as pdt
import pytest


@pytest.fixture
def bamfile(tasic2016_bam):
    return os.path.join(
        tasic2016_bam,
        'CAV_LP_Ipsi_tdTpos_cell_1-SRR2140356-GSM1840944_R1.polyATrim.'
        'adapterTrim.rmRep.sorted.rg.subset.sorted.bam')


@pytest.fixture
def junction_reads_table_csv(bamfile, tasic2016_intermediate_bam,
                             ignore_multimapping):
    """csv for bamfile (singular)"""
    basename = os.path.basename(bamfile)
    suffix = '.outrigger_junction_reads_ignore-multimapping{}.csv'.format(
        ignore_multimapping)
    basename2 = basename + suffix
    csv = os.path.join(tasic2016_intermediate_bam, basename2)
    return csv

@pytest.fixture
def bamfiles(tasic2016_bam):
    return glob.glob(os.path.join(tasic2016_bam, '*'))


@pytest.fixture
def junction_reads_table_csvs(tasic2016_intermediate_bam,
                             ignore_multimapping):
    """csvs for bamfiles (multiple)"""
    suffix = '*.outrigger_junction_reads_ignore-multimapping{}.csv'.format(
        ignore_multimapping)
    csvs = glob.glob(os.path.join(tasic2016_intermediate_bam, suffix))
    return csvs


def test__report_read_positions():
    pass


def test__choose_strand_and_sum():
    pass


def test__reads_dict_to_table():
    pass


def test__get_junction_reads():
    pass


def test_bam_to_junction_reads_table(bamfile, junction_reads_table_csv):
    from outrigger.io.bam import bam_to_junction_reads_table

    test = bam_to_junction_reads_table(bamfile)
    true = pd.read_csv(junction_reads_table_csv)

    pdt.assert_frame_equal(test, true)


def test_read_multiple_bams(bamfiles, junction_reads_table_csvs):
    from outrigger.io.bam import read_multiple_bams

    test = read_multiple_bams(bamfiles)

    dfs = [pd.read_csv(csv) for csv in junction_reads_table_csvs]
    true = pd.concat(dfs, ignore_index=True)

    pdt.assert_frame_equal(test, true)
