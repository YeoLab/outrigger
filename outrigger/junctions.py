import argparse
import datetime
import os
import sys

import gffutils
import pandas as pd

UPSTREAM = 'upstream'
DOWNSTREAM = 'downstream'
DIRECTIONS = UPSTREAM, DOWNSTREAM


class JunctionAnnotator(object):
    """Annotate junctions with adjacent exons"""

    def __init__(self, splice_junctions, db, junction_id='junction_id',
                 exon_start='exon_start', exon_stop='exon_stop',
                 chrom='chrom', strand='strand'):
        """Initialize class to get upstream/downstream exons of junctions

        Parameters
        ----------
        splice_junctions : pandas.DataFrame
            A table of splice junctions with the columns indicated by the
            variables `junction_id`, `exon_start`, `exon_stop`, `chrom`,
            `strand`
        db : gffutils.FeatureDB
            Gffutils Database of gene, transcript, and exon features.
        junction_id, exon_start, exon_stop, chrom, strand : str
            Columns in `splice_junctions`
        """
        columns = self.junction_id, self.exon_start, self.exon_stop, \
                  self.chrom, self.strand

        for column in columns:
            if column not in splice_junctions:
                raise ValueError('The required column {} is not in the splice '
                                 'junction dataframe'.format(column))

        self.splice_junctions = splice_junctions.set_index(junction_id)
        self.splice_junctions = self.splice_junctions.sort_index()

        self.junction_id = junction_id
        self.exon_start = exon_start
        self.exon_stop = exon_stop
        self.chrom = chrom
        self.strand = strand

        self.db = db

    @staticmethod
    def _single_junction_exon_triple(direction_ind, direction, exon_id):
        """Create exon, direction, junction triples for an exon + its junctions

        Parameters
        ----------
        direction_ind : pandas.Series
            A boolean series of the indices of the junctions matching with the
            provided exon. The index of the series must be the junction ID
        direction : str
            The direction of the exon relative to the junctions, either
            "upstream" or "downstream"
        exon_id : str
            Unique identifier of the exon

        Returns
        -------


        """
        length = direction_ind.sum()

        exons = [exon_id] * length
        directions = [direction] * length
        junctions = direction_ind[direction_ind].index
        return pd.DataFrame(zip(exons, directions, junctions),
                            columns=['exon', 'direction', 'junction'])

    @staticmethod
    def genome_to_transcript_adjacency(adjacent_in_genome, strand):
        """If negative strand, swap the upstream/downstream adjacency"""
        if strand == '+':
            return {UPSTREAM: adjacent_in_genome[UPSTREAM],
                    DOWNSTREAM: adjacent_in_genome[DOWNSTREAM]}
        elif strand == '-':
            return {UPSTREAM: adjacent_in_genome[DOWNSTREAM],
                    DOWNSTREAM: adjacent_in_genome[UPSTREAM]}


    def genome_adjacent(self, exon, exon_start='exon_start',
                        exon_stop='exon_stop', chrom='chrom', strand='strand'):
        """Get indices of junctions next to an exon, in genome coordinates"""
        chrom_ind = self.splice_junctions[chrom] == exon.chrom

        strand_ind = self.splice_junctions[strand] == exon.strand

        upstream_in_genome = chrom_ind & strand_ind \
                             & (self.splice_junctions[exon_stop] == exon.stop)
        downstream_in_genome = chrom_ind & strand_ind \
                               & (self.splice_junctions[exon_start] == exon.start)
        return {UPSTREAM: upstream_in_genome, DOWNSTREAM: downstream_in_genome}


    def adjacent_junctions(self, exon, exon_start='exon_start',
                           exon_stop='exon_stop', chrom='chrom',
                           strand='strand'):
        """Get junctions adjacent to this exon"""
        dfs = []
        adjacent_in_genome = self.genome_adjacent(exon, self.splice_junctions,
                                                  exon_start, exon_stop, chrom,
                                                  strand)
        adjacent_in_transcriptome = self.genome_to_transcript_adjacency(
            adjacent_in_genome, exon.strand)

        exon_id = exon.id
        for direction, ind in adjacent_in_transcriptome.items():
            if ind.any():
                df = self._single_junction_exon_triple(ind, direction, exon_id)
                dfs.append(df)

        if len(dfs) > 0:
            return pd.concat(dfs, ignore_index=True)
        else:
            return pd.DataFrame()

    def get_adjacent_exons(self, exon_start='exon_start',
                           exon_stop='exon_stop', chrom='chrom'):
        """Get upstream and downstream exons of each junction

        Use junctions defined in ``sj_metadata`` and exons in ``db`` to create
        triples of (exon, direction, junction), which are read like
        (subject, object, verb) e.g. ('exon1', 'upstream', 'junction12'), for
        creation of a graph database.

        Parameters
        ----------
        sj_metadata : pandas.DataFrame
            A splice junction metadata dataframe with the junction id as the
            index, with  columns defined by variables ``exon_start`` and
            ``exon_stop``.
        db : gffutils.FeatureDB
            A database of gene annotations created by gffutils. Must have features
            of type "exon"
        exon_start : str, optional
            Name of the column in sj_metadata corresponding to the start of the
            exon
        exon_stop : str, optional
            Name of the column in sj_metadata corresponding to the end of the exon

        Returns
        -------
        junction_exon_triples : pandas.DataFrame
            A three-column dataframe describing the relationship of where an exon
            is relative to junctions
        """
        n_exons = sum(1 for _ in self.db.features_of_type('exon'))

        dfs = []

        sys.stdout.write('Starting annotation of all junctions with known '
                         'exons...\n')
        for i, exon in enumerate(self.db.features_of_type('exon')):
            if (i + 1) % 10000 == 0:
                sys.stdout.write('\t{}/{} exons completed\n'.format(i + 1,
                                                                    n_exons))
            df = self.adjacent_junctions(exon)
        junction_exon_triples = pd.concat(dfs, ignore_index=True)
        sys.stdout.write('Done.\n')
        return junction_exon_triples
