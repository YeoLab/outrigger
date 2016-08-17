"""
Find exons adjacent to junctions
"""

import pandas as pd

from ..io.common import JUNCTION_ID, EXON_START, EXON_STOP, CHROM, STRAND
from ..util import done, progress


UPSTREAM = 'upstream'
DOWNSTREAM = 'downstream'
DIRECTIONS = UPSTREAM, DOWNSTREAM


class ExonJunctionAdjacencies(object):
    """Annotate junctions with neighboring exon_cols (upstream or downstream)"""

    def __init__(self, metadata, db, junction_id=JUNCTION_ID,
                 exon_start=EXON_START, exon_stop=EXON_STOP,
                 chrom=CHROM, strand=STRAND):
        """Initialize class to get upstream/downstream exon_cols of junctions

        Parameters
        ----------
        metadata : pandas.DataFrame
            A table of splice junctions with the columns indicated by the
            variables `junction_id`, `exon_start`, `exon_stop`, `chrom`,
            `strand`
        db : gffutils.FeatureDB
            Gffutils Database of gene, transcript, and exon features.
        junction_id, exon_start, exon_stop, chrom, strand : str
            Columns in `metadata`
        """

        columns = junction_id, exon_start, exon_stop, chrom, strand

        for column in columns:
            if column not in metadata:
                raise ValueError('The required column {} is not in the splice '
                                 'junction dataframe'.format(column))

        self.metadata = metadata.set_index(junction_id)
        self.metadata = self.metadata.sort_index()

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
        triples : pandas.DataFrame
            A (n, 3) sized dataframe of an exon and its adjacent junctions
        """
        length = direction_ind.sum()

        exons = [exon_id] * length
        directions = [direction] * length
        junctions = direction_ind[direction_ind].index
        return pd.DataFrame(list(zip(exons, directions, junctions)),
                            columns=['exon', 'direction', 'junction'])

    @staticmethod
    def _to_stranded_transcript_adjacency(adjacent_in_genome, strand):
        """If negative strand, swap the upstream/downstream adjacency

        Parameters
        ----------
        adjacent_in_genome : dict
            dict of two keys, "upstream" and "downstream", mapping to a boolean
            series indicating whether the junction is upstream or downstream of
            a particular exon
        strand : "-" | "+"
            Positive or negative strand
        """
        if strand == '+':
            return {UPSTREAM: adjacent_in_genome[UPSTREAM],
                    DOWNSTREAM: adjacent_in_genome[DOWNSTREAM]}
        elif strand == '-':
            return {UPSTREAM: adjacent_in_genome[DOWNSTREAM],
                    DOWNSTREAM: adjacent_in_genome[UPSTREAM]}

    def _junctions_genome_adjacent_to_exon(self, exon):
        """Get indices of junctions next to an exon, in genome coordinates"""
        chrom_ind = self.metadata[self.chrom] == exon.chrom

        strand_ind = self.metadata[self.strand] == exon.strand

        upstream_in_genome = \
            chrom_ind & strand_ind \
            & (self.metadata[self.exon_stop] == exon.stop)
        downstream_in_genome = \
            chrom_ind & strand_ind & \
            (self.metadata[self.exon_start] == exon.start)
        return {UPSTREAM: upstream_in_genome, DOWNSTREAM: downstream_in_genome}

    def _adjacent_junctions_single_exon(self, exon):
        """Get junctions adjacent to this exon

        Parameters
        ----------
        exon : gffutils.Feature
            An item in

        """
        dfs = []
        adjacent_in_genome = self._junctions_genome_adjacent_to_exon(exon)
        adjacent_in_transcriptome = self._to_stranded_transcript_adjacency(
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

    def neighboring_exons(self):
        """Get upstream and downstream exon_cols of each junction

        The "upstream" and "downstream" is relative to the **junction**, e.g.

            exonA   upstream    junctionX
            exonB   downstream    junctionX

        should be read as "exonA is upstream of juction X" and "exonB is
        downstream of junctionX"

        Use junctions defined in ``sj_metadata`` and exon_cols in ``db`` to create
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
            A database of gene annotations created by gffutils. Must have
            features of type "exon"
        exon_start : str, optional
            Name of the column in sj_metadata corresponding to the start of the
            exon
        exon_stop : str, optional
            Name of the column in sj_metadata corresponding to the end of the
            exon

        Returns
        -------
        junction_exon_triples : pandas.DataFrame
            A three-column dataframe describing the relationship of where an
            exon is relative to junctions
        """
        n_exons = sum(1 for _ in self.db.features_of_type('exon'))

        dfs = []

        progress('Starting annotation of all junctions with known '
                 'neighboring exon_cols ...')
        for i, exon in enumerate(self.db.features_of_type('exon')):
            if (i + 1) % 10000 == 0:
                progress('\t{}/{} exon_cols completed'.format(i + 1, n_exons))
            df = self._adjacent_junctions_single_exon(exon)
            dfs.append(df)
        junction_exon_triples = pd.concat(dfs, ignore_index=True)
        done()
        return junction_exon_triples
