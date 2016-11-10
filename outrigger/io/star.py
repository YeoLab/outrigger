"""
Read splice junction output files from STAR aligner (SJ.out.tab)
"""
import os

import joblib
import numpy as np
import pandas as pd

from ..common import JUNCTION_ID, JUNCTION_START, JUNCTION_STOP, READS, \
    JUNCTION_MOTIF, EXON_START, EXON_STOP, CHROM, STRAND, ANNOTATED, \
    SAMPLE_ID, UNIQUE_READS, MULTIMAP_READS, MAX_OVERHANG

from .core import add_exons_and_junction_ids

COLUMN_NAMES = (CHROM, JUNCTION_START, JUNCTION_STOP, STRAND,
                JUNCTION_MOTIF, ANNOTATED, UNIQUE_READS, MULTIMAP_READS,
                MAX_OVERHANG)


def int_to_junction_motif(n):
    if n == 0:
        return 'non-canonical'
    if n == 1:
        return 'GT/AG'
    if n == 2:
        # Negative strand: CT/AC
        return 'GT/AG'
    if n == 3:
        return 'GC/AG'
    if n == 4:
        # Negative strand: CT/GC
        return 'GC/AG'
    if n == 5:
        return 'AT/AC'
    if n == 6:
        # Negative strand: GT/AT
        return 'AT/AC'


def read_sj_out_tab(filename):
    """Read an SJ.out.tab file as produced by the RNA-STAR aligner into a
    pandas Dataframe

    Parameters
    ----------
    filename : str of filename or file handle
        Filename of the SJ.out.tab file you want to read in

    Returns
    -------
    sj : pandas.DataFrame
        Dataframe of splice junctions with the columns,
        ('chrom', 'junction_start', 'junction_stop', 'strand',
        'junction_motif', 'exon_start', 'exon_stop', 'annotated',
        'unique_junction_reads', 'multimap_junction_reads', 'max_overhang')

    """
    sj = pd.read_table(filename, header=None, names=COLUMN_NAMES, sep='\s+')
    sj[JUNCTION_MOTIF] = sj[JUNCTION_MOTIF].map(int_to_junction_motif)

    # Convert unknown strands to explicitly say "undefined"
    sj[STRAND] = sj[STRAND].replace(0, 'undefined')

    # Convert integer strand to symbol
    # Use index-based replacement because it's 100x faster than map
    rows = sj.strand == 1
    sj.loc[rows, STRAND] = '+'
    rows = sj.strand == 2
    sj.loc[rows, STRAND] = '-'

    # Translate negative strand intron motifs
    # rows = sj.strand == '-'
    # sj.loc[rows, 'intron_motif'] = sj.intron_motif[rows].map(
    #     lambda x: NEG_STRAND_INTRON_MOTIF[x])
    sj.annotated = sj.annotated.astype(bool)

    sj = add_exons_and_junction_ids(sj)

    return sj


def _read_single_filename(filename, sample_id_func, ignore_multimapping=False):
    splice_junction = read_sj_out_tab(filename)
    sample_id = sample_id_func(filename)
    sample_id = sample_id.split('SJ.out.tab')[0].rstrip('.')
    splice_junction[SAMPLE_ID] = sample_id

    if not ignore_multimapping:
        splice_junction[READS] = splice_junction[UNIQUE_READS] \
                                  + splice_junction[MULTIMAP_READS]
    else:
        splice_junction[READS] = splice_junction[UNIQUE_READS]
    return splice_junction


def read_multiple_sj_out_tab(filenames, ignore_multimapping=False,
                             sample_id_func=os.path.basename, n_jobs=-1):
    """Read the splice junction files and return a tall, tidy dataframe

    Adds a column called "sample_id" based on the basename of the file, minus
    "SJ.out.tab"

    Parameters
    ----------
    filenames : iterator
        A list or other iterator of filenames to read
    multimapping : bool
        If True, include the multimapped reads in total read count
    sample_id_func : function
        A function to extract the sample id from the filenames

    Returns
    -------
    metadata : pandas.DataFrame
        A tidy dataframe, where each row has the observed reads for a sample
    """
    dfs = joblib.Parallel(n_jobs=n_jobs)(
        joblib.delayed(_read_single_filename)(
            filename, sample_id_func, ignore_multimapping)
        for filename in filenames)
    splice_junctions = pd.concat(dfs, ignore_index=True)

    splice_junctions[CHROM] = splice_junctions[CHROM].astype(str)
    splice_junctions = splice_junctions.sort_values(
        by=[SAMPLE_ID, CHROM, JUNCTION_START, JUNCTION_STOP])
    splice_junctions.index = np.arange(splice_junctions.shape[0])
    return splice_junctions


def make_metadata(spliced_reads, columns=(JUNCTION_ID, CHROM, JUNCTION_START,
                                          JUNCTION_STOP, STRAND, ANNOTATED,
                                          EXON_START, EXON_STOP)):
    """Get barebones junction chrom, start, stop, strand information

    Parameters
    ----------
    spliced_reads : pandas.DataFrame
        Concatenated SJ.out.tab files created by read_sj_out_tab
    columns : iterable
        Which columns to use to make the metadata

    Returns
    -------
    junctions : pandas.DataFrame
        A (n_junctions, 9) dataframe containing the columns:
         - junction_id
         - chrom
         - intron_start
         - intron_stop
         - exon_start
         - exon_stop
         - strand
         - intron_motif
         - annotated
    """
    columns = spliced_reads.columns.intersection(columns)
    metadata = spliced_reads[columns]
    metadata = metadata.drop_duplicates()
    metadata.index = np.arange(metadata.shape[0])

    # Force all chromosomes to be strings to be compatible with ENSEMBL's
    # integer chromosomes
    metadata[CHROM] = metadata[CHROM].astype(str)

    return metadata
