"""
Read splice junction output files from STAR aligner (SJ.out.tab)
"""
import os

import pandas as pd

COLUMN_NAMES = ('chrom', 'intron_start', 'intron_stop', 'strand',
                'intron_motif', 'annotated',
                'unique_junction_reads', 'multimap_junction_reads',
                'max_overhang')
NEG_STRAND_INTRON_MOTIF = {'CT/AC': 'GT/AG',
                           'CT/GC': 'GC/AG',
                           'GT/AT': 'AT/AC',
                           'non-canonical': 'non-canonical'}


def int_to_intron_motif(n):
    if n == 0:
        return 'non-canonical'
    if n == 1:
        return 'GT/AG'
    if n == 2:
        return 'CT/AC'
    if n == 3:
        return 'GC/AG'
    if n == 4:
        return 'CT/GC'
    if n == 5:
        return 'AT/AC'
    if n == 6:
        return 'GT/AT'


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
        ('chrom', 'intron_start', 'intron_stop', 'strand',
        'intron_motif', 'annotated', 'unique_junction_reads',
        'multimap_junction_reads', 'max_overhang')

    """
    sj = pd.read_table(filename, header=None, names=COLUMN_NAMES, sep='\s+')
    sj.intron_motif = sj.intron_motif.map(int_to_intron_motif)

    # Convert integer strand to symbol
    # Use index-based replacement because it's 100x faster than map
    rows = sj.strand == 1
    sj.loc[rows, 'strand'] = '+'
    rows = sj.strand == 2
    sj.loc[rows, 'strand'] = '-'

    # Translate negative strand intron motifs
    rows = sj.strand == '-'
    sj.loc[rows, 'intron_motif'] = sj.intron_motif[rows].map(
        lambda x: NEG_STRAND_INTRON_MOTIF[x])
    sj.annotated = sj.annotated.astype(bool)

    # Add intron location
    sj['junction_location'] = sj.chrom.astype(str) + ':' \
                              + sj.intron_start.astype(str) + '-' \
                              + sj.intron_stop.astype(str) + ':' \
                              + sj.strand.astype(str)

    return sj


def read_multiple_sj_out_tab(filenames, sample_id_func=os.path.basename):
    """Read the splice junction files and return a tall, tidy dataframe

    Adds a column called "sample_id" based on the basename of the file, minus
    "SJ.out.tab"

    Parameters
    ----------
    filenames : iterator
        A list or other iterator of filenames to read
    sample_id_func : function
        A function to extract the sample id from the filenames

    Returns
    -------
    splice_junctions : pandas.DataFrame
        A tidy dataframe, where each row has the observed reads for a sample
    """
    splice_junctions = []
    for filename in filenames:
        splice_junction = read_sj_out_tab(filename)
        sample_id = sample_id_func(filename)
        sample_id = sample_id.split('SJ.out.tab')[0].rstrip('.')
        splice_junction['sample_id'] = sample_id
    splice_junctions = pd.concat(splice_junctions, ignore_index=True)
    splice_junctions = splice_junctions.set_index('junction_location')
    return splice_junctions


def sj_count_to_metadata(splice_junctions):
    """Get just the junction metadata information"""
    metadata = splice_junctions.drop(['unique_junction_reads',
                                      'multimap_junction_reads',
                                      'max_overhang', 'sample_id'],
                                     axis=1)
    metadata = metadata.drop_duplicates()
    return metadata

