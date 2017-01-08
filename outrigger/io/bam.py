import collections
import os

import joblib
import numpy as np
import pandas as pd
import pysam

from ..common import UNIQUE_READS, MULTIMAP_READS, READS, CHROM, \
    JUNCTION_START, JUNCTION_STOP, STRAND
from .core import add_exons_and_junction_ids


def _report_read_positions(read, counter, stranded=False):
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
            if stranded:
                location = (chrom, start, stop, strand)
            else:
                location = (chrom, start, stop)
            counter[location] += 1
            del start
            del stop
        last_read_pos = read_loc


def _combine_uniquely_multi(uniquely, multi, ignore_multimapping=False,
                            stranded=False):
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

    renamer = {'level_0': CHROM, 'level_1': JUNCTION_START,
               'level_2': JUNCTION_STOP}
    if stranded:
        renamer['level_3'] = STRAND

    reads = reads.rename(columns=renamer)
    reads.index = np.arange(reads.shape[0])
    return reads


def _get_junction_reads(filename, stranded):
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

            _report_read_positions(read, counter, stranded)
    samfile.close()
    return uniquely, multi


def bam_to_junction_reads_table(bam_filename, ignore_multimapping=False,
                                stranded=False):
    """Create a table of reads for this bam file"""
    uniquely, multi = _get_junction_reads(bam_filename, stranded)
    reads = _combine_uniquely_multi(uniquely, multi, ignore_multimapping,
                                    stranded)

    # Remove "junctions" with same start and stop
    reads = reads.loc[reads[JUNCTION_START] != reads[JUNCTION_STOP]]
    reads.index = np.arange(reads.shape[0])

    reads['sample_id'] = os.path.basename(bam_filename)
    reads = add_exons_and_junction_ids(reads, stranded)
    return reads


def read_multiple_bams(bam_filenames, ignore_multimapping=False,
                       stranded=False, n_jobs=-1):
    dfs = joblib.Parallel(n_jobs=n_jobs)(
        joblib.delayed(
            bam_to_junction_reads_table)(filename, ignore_multimapping,
                                         stranded)
        for filename in bam_filenames)
    reads = pd.concat(dfs, ignore_index=True)
    return reads
