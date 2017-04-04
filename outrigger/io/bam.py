import collections
import os

import joblib
import numpy as np
import pandas as pd
import pysam

from ..common import UNIQUE_READS, MULTIMAP_READS, READS, CHROM, \
    JUNCTION_START, JUNCTION_STOP, STRAND
from .core import add_exons_and_junction_ids


def _report_read_positions(read, counter):
    chrom = read.reference_name
    strand = '-' if read.is_reverse else '+'

    last_read_pos = False
    for read_loc, genome_loc in read.get_aligned_pairs():
        if read_loc is None and last_read_pos:
            # Add one to be compatible with STAR output and show the
            # start of the intron (not the end of the exon)
            start = genome_loc + 1
        elif read_loc and last_read_pos is None:
            stop = genome_loc  # we are right exclusive ,so this is correct
            counter[(chrom, start, stop, strand)] += 1
            del start
            del stop
        last_read_pos = read_loc


def _choose_strand_and_sum(reads):
    """Use the strand with more counts and sum all reads with same junction

    STAR seems to take a simple majority to decide on strand when there are
    reads mapping to both, so we'll do the same

    Parameters
    ----------
    reads : pandas.Series
        A (chrom, start, stop, strand)-indexed series of read counts

    Returns
    -------
    reads_strand_chosen : pandas.Series
        A (chrom, start, stop, strand)-indexed series of read counts, with
        the majority strand as the "winner" and

    """
    if reads.empty:
        return pd.Series(name=reads.name)
    locations = reads.groupby(level=(0, 1, 2)).idxmax()
    counts = reads.groupby(level=(0, 1, 2)).sum()

    index = pd.MultiIndex.from_tuples(locations.values)

    return pd.Series(counts.values, index=index, name=reads.name)


def _combine_uniquely_multi(uniquely, multi, ignore_multimapping=False):
    """Combine uniquely and multi-mapped read counts into a single table

    Parameters
    ----------
    unqiuely, multi : dict
        A dictionary of {(chrom, start, end, strand) : n_reads} uniquely mapped
        and multi-mapped (reads that could map to multiple parts of the genome)
    ignore_multimapping : bool
        When summing all reads, whether or not to ignore the multimapping
        reads. Default is False.

    Returns
    -------
    reads : pandas.DataFrame
        A combined table of all uniquely and multi-mapped reads, with an
        additional column of "reads" which will ultimately be the reads used
        for creating an outrigger index and calculating percent spliced-in.
    """
    uniquely = pd.Series(uniquely, name=UNIQUE_READS)
    multi = pd.Series(multi, name=MULTIMAP_READS)

    uniquely = _choose_strand_and_sum(uniquely)
    multi = _choose_strand_and_sum(multi)

    # Join the data on the chromosome locations
    if multi.empty:
        reads = uniquely.to_frame()
        reads[MULTIMAP_READS] = np.nan
    elif uniquely.empty:
        reads = multi.to_frame()
        reads[UNIQUE_READS] = np.nan
    else:
        reads = uniquely.to_frame().join(multi)

    reads = reads.fillna(0)
    reads = reads.astype(int)

    if ignore_multimapping:
        reads[READS] = reads[UNIQUE_READS]
    else:
        reads[READS] = reads.sum(axis=1)
    reads = reads.reset_index()
    reads = reads.rename(columns={'level_0': CHROM, 'level_1': JUNCTION_START,
                                  'level_2': JUNCTION_STOP, 'level_3': STRAND})
    reads.index = np.arange(reads.shape[0])
    return reads


def _get_junction_reads(filename):
    """Read a sam file and extract unique and multi mapped junction reads"""
    samfile = pysam.AlignmentFile(filename, "rb")

    # Uniquely mapped reads
    uniquely = collections.Counter()

    # Multimapped reads
    multi = collections.Counter()

    for read in samfile.fetch():
        if "N" in read.cigarstring:
            if read.mapping_quality < 255:
                counter = multi
            else:
                counter = uniquely

            _report_read_positions(read, counter)
    samfile.close()
    return uniquely, multi


def bam_to_junction_reads_table(bam_filename, ignore_multimapping=False):
    """Create a table of reads for this bam file"""
    uniquely, multi = _get_junction_reads(bam_filename)
    reads = _combine_uniquely_multi(uniquely, multi, ignore_multimapping)

    # Remove "junctions" with same start and stop
    reads = reads.loc[reads[JUNCTION_START] != reads[JUNCTION_STOP]]
    reads.index = np.arange(reads.shape[0])

    reads['sample_id'] = os.path.basename(bam_filename)
    reads = add_exons_and_junction_ids(reads)
    return reads


def read_multiple_bams(bam_filenames, ignore_multimapping=False, n_jobs=-1):
    dfs = joblib.Parallel(n_jobs=n_jobs)(
        joblib.delayed(
            bam_to_junction_reads_table)(filename, ignore_multimapping)
        for filename in bam_filenames)
    reads = pd.concat(dfs, ignore_index=True)
    return reads
