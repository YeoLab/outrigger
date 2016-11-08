import collections
import os

import pandas as pd
import pysam

from ..common import UNIQUE_READS, MULTIMAP_READS, READS, CHROM, \
    JUNCTION_START, JUNCTION_STOP, STRAND


def report_read_position(read, counter):
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


def choose_strand_and_sum(reads):
    """Use the strand with more counts and sum all reads with same junction

    STAR seems to take a simple majority to decide on strand when there are
    reads mapping to both, so we'll do the same
    """
    locations = reads.groupby(level=(0, 1, 2)).idxmax()
    counts = reads.groupby(level=(0, 1, 2)).sum()

    index = pd.MultiIndex.from_tuples(locations.values)

    return pd.Series(counts.values, index=index, name=reads.name)


def reads_dict_to_table(uniquely, multi, ignore_multimapping=False):
    uniquely = pd.Series(uniquely, name=UNIQUE_READS)
    multi = pd.Series(multi, name=MULTIMAP_READS)

    uniquely = choose_strand_and_sum(uniquely)
    multi = choose_strand_and_sum(multi)

    # Join the data on the chromosome locations
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

    return reads


def get_junction_reads(filename):
    samfile = pysam.AlignmentFile(filename, "rb")

    # Uniquely mapped reads
    uniquely = collections.Counter()

    # Multimapped reads
    multi = collections.Counter()

    for read in samfile.fetch():
        if "N" in read.cigarstring:
            if read.is_secondary:
                counter = multi
            else:
                counter = uniquely

            report_read_position(read, counter)
    samfile.close()
    return uniquely, multi


def make_reads_table(filename, ignore_multimapping):
    uniquely, multi = get_junction_reads(filename)
    reads = reads_dict_to_table(uniquely, multi, ignore_multimapping)

    reads['sample_id'] = os.path.basename(filename)
    return reads
