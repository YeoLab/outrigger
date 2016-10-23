"""
Find exons adjacent to junctions
"""
import sqlite3
import warnings

import gffutils
from gffutils.helpers import merge_attributes
import joblib

from ..common import JUNCTION_ID, EXON_START, EXON_STOP, CHROM, STRAND
from ..io.gtf import transform, maybe_analyze
from ..region import Region, STRANDS
from ..util import done, progress


with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import pandas as pd


UPSTREAM = 'upstream'
DOWNSTREAM = 'downstream'
DIRECTIONS = UPSTREAM, DOWNSTREAM
NOVEL_EXON = 'novel_exon'
OUTRIGGER_DE_NOVO = 'outrigger_de_novo'
MAX_DE_NOVO_EXON_LENGTH = 100


def _unify_strand(strand1, strand2):
    """If strands are equal, return the strand, otherwise return None"""
    if strand1 != strand2:
        strand = '.'
    else:
        strand = strand1
    return strand


def _exons_from_neighboring_junctions(junction, neighbors, side='left'):
    """Returns (chrom, start, stop, strand) of neighboring exons from junctions

    Parameters
    ----------
    junction : outrigger.Region
        A single junction whose neighboring exons you want
    neighbors : pandas.DataFrame
        All junctions which are within

    Returns
    -------
    exons : pandas.Series
        A column containing (chrom, start, stop, strand) for each detected exon
    """
    if neighbors.empty:
        return pd.Series()
    if side == 'left':
        return neighbors.apply(lambda x: (
            x.chrom, x.stop + 1, junction.start - 1,
            _unify_strand(x.strand, junction.strand)), axis=1)
    elif side == 'right':
        return neighbors.apply(lambda x: (
            x.chrom, junction.stop + 1, x.start - 1,
            _unify_strand(x.strand, junction.strand)), axis=1)


def _neighboring_exons(junction, df, side='left',
                       max_de_novo_exon_length=MAX_DE_NOVO_EXON_LENGTH):
    """Get either the left or right neighbors of a particular junction

    Used to find novel exons between junctions. Not part of the
    ExonJunctionAdjacencies object so it can be parallelized with joblib.
    Internal function

    Parameters
    ----------
    junction : outrigger.Region
        A Region object with .start and .stop
    df : pandas.DataFrame
        A data table with "start" and "stop" columns of all junctions
    side : 'left' | 'right'
        Specifies whether you want the left or right side neighbors of a
        junction (not strand-specific)

    Returns
    -------
    neighboring_exon_locations : pandas.Series
        A column containing (chrom, start, stop, strand) for each detected exon
    """
    if side == 'left':
        rows = junction.start - df.stop
    elif side == 'right':
        rows = df.start - junction.stop
    rows = rows[rows > 0]
    rows = rows[rows <= max_de_novo_exon_length]
    neighboring_junctions = df.loc[rows.index]

    neighboring_exon_locations = _exons_from_neighboring_junctions(
        junction, neighboring_junctions, side=side)

    return neighboring_exon_locations


def is_there_an_exon_here(self, junction1, junction2):
    """Check if there could be an exon between these two junctions

    Parameters
    ----------
    junction{1,2} : outrigger.Region
        Outrigger.Region objects

    Returns
    -------
    start, stop : (int, int) or (False, False)
        Start and stop of the new exon if it exists, else False, False
    """
    if junction1.overlaps(junction2):
        return False, False

    # These are junction start/stops, not exon start/stops
    # Move one nt upstream of starts for exon stops,
    # and one nt downstream of stops for exon starts.
    option1 = abs(junction1.stop - junction2.start) \
        < self.max_de_novo_exon_length
    option2 = abs(junction2.stop - junction1.start) \
        < self.max_de_novo_exon_length

    if option1:
        return junction1.stop + 1, junction2.start - 1
    elif option2:
        return junction2.stop + 1, junction1.start - 1
    return False, False


class ExonJunctionAdjacencies(object):
    """Annotate junctions with neighboring exons (upstream or downstream)"""

    exon_types = 'exon', NOVEL_EXON

    def __init__(self, metadata, db, junction_id=JUNCTION_ID,
                 exon_start=EXON_START, exon_stop=EXON_STOP,
                 chrom=CHROM, strand=STRAND,
                 max_de_novo_exon_length=MAX_DE_NOVO_EXON_LENGTH,
                 n_jobs=-1):
        """Initialize class to get upstream/downstream exons of junctions

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
        progress('\tLooking up which exons are already defined ...')
        self.existing_exons = set(
            i['id'] for i in self.db.execute(
                'select id from features where featuretype = "exon"'))
        done(n_tabs=3)
        self.max_de_novo_exon_length = max_de_novo_exon_length

        self.n_jobs = n_jobs

    def detect_exons_from_junctions(self):
        """Find exons based on gaps in junctions"""
        junctions = pd.Series(self.metadata.index.map(Region),
                              name='region', index=self.metadata.index)
        junctions = junctions.to_frame()
        junctions['chrom'] = junctions['region'].map(lambda x: x.chrom)
        junctions['start'] = junctions['region'].map(lambda x: x.start)
        junctions['stop'] = junctions['region'].map(lambda x: x.stop)
        junctions['strand'] = junctions['region'].map(lambda x: x.strand)

        for chrom, df in junctions.groupby('chrom'):
            # Only get left-adjacent novel exons since there has to be a
            # junction on both sides, and since we iterate over ALL junctions,
            # if we get all left and right exons for all junctions, we're
            # double-counting exons
            progress('\tFinding all exons on chromosome {chrom} '
                     '...'.format(chrom=chrom))
            exon_locations = pd.concat(joblib.Parallel(n_jobs=self.n_jobs)(
                joblib.delayed(_neighboring_exons)(junction, df, 'left')
                for junction in df.region))
            done(n_tabs=3)

            progress('\t\tFiltering for only novel exons on chromosome '
                     '{chrom} ...'.format(chrom=chrom))
            novel_exons = set(x for x in exon_locations if
                              'exon:{}:{}-{}:{}'.format(*x)
                              not in self.existing_exons)
            done(n_tabs=4)

            progress('\t\tCreating gffutils.Feature objects for each novel '
                     'exon, plus potentially its overlapping gene')
            exon_features = [self.exon_location_to_feature(*x)
                             for x in novel_exons]
            done(n_tabs=4)

            progress('\t\tUpdating gffutils database with {n} novel exons on '
                     'chromosome {chrom} ...'.format(chrom=chrom,
                                                     n=len(novel_exons)))
            try:
                self.db.update(exon_features,
                               make_backup=False,
                               id_spec={'novel_exon': 'location_id'},
                               transform=transform)
            except ValueError:
                progress('\tNo novel exons found on chromosome '
                         '{chrom}'.format(chrom=chrom))
            done(n_tabs=4)

        # For up to 1000x faster queries, re-Analyze the database now that it
        # has been updated
        maybe_analyze(self.db)

    def exon_location_to_feature(self, chrom, start, stop, strand):
        if strand not in STRANDS:
            strand = '.'
        overlapping_genes = self.db.region(seqid=chrom, start=start,
                                           end=stop, strand=strand,
                                           featuretype='gene')

        exon_id = 'exon:{chrom}:{start}-{stop}:{strand}'.format(
            chrom=chrom, start=start, stop=stop, strand=strand)

        attributes = {}
        for g in overlapping_genes:
            attributes = merge_attributes(attributes, g.attributes)

        exon = gffutils.Feature(chrom, source=OUTRIGGER_DE_NOVO,
                                featuretype=NOVEL_EXON, start=start,
                                end=stop, strand=strand, id=exon_id,
                                attributes=attributes)
        return exon

    def write_de_novo_exons(self, filename='novel_exons.gtf'):
        """Write all de novo exons to a gtf"""
        with open(filename, 'w') as f:
            for noveL_exon in self.db.features_of_type(NOVEL_EXON):
                f.write(str(noveL_exon) + '\n')

    def add_exon_to_db(self, chrom, start, stop, strand):
        if strand not in STRANDS:
            strand = None
        overlapping_genes = list(self.db.region(seqid=chrom, start=start,
                                                end=stop, strand=strand,
                                                featuretype='gene'))

        exon_id = 'exon:{chrom}:{start}-{stop}:{strand}'.format(
            chrom=chrom, start=start, stop=stop, strand=strand)

        if len(overlapping_genes) == 0:
            exon = gffutils.Feature(chrom, source=OUTRIGGER_DE_NOVO,
                                    featuretype=NOVEL_EXON, start=start,
                                    end=stop, strand=strand, id=exon_id)
            progress('\tAdded a novel exon ({}), located in an unannotated'
                     ' gene'.format(exon.id))
            self.db.update([exon], id_spec={'novel_exon': 'location_id'},
                           transform=transform)
            return

        de_novo_exons = [gffutils.Feature(
            chrom, source=OUTRIGGER_DE_NOVO, featuretype=NOVEL_EXON,
            start=start, end=stop, strand=g.strand, id=exon_id + g.strand,
            attributes=dict(g.attributes.items()))
                         for g in overlapping_genes]

        # Add all exons that aren't already there
        for exon in de_novo_exons:
            try:
                try:
                    gene_name = ','.join(exon.attributes['gene_name'])
                except KeyError:
                    try:
                        gene_name = ','.join(exon.attributes['gene_id'])
                    except KeyError:
                        gene_name = 'unknown_gene'
                try:
                    # Check that the non-novel exon doesn't exist already
                    self.db[exon_id + exon.strand]
                except gffutils.FeatureNotFoundError:
                    self.db.update([exon],
                                   id_spec={'novel_exon': 'location_id'},
                                   transform=transform)
                    progress(
                        '\tAdded a novel exon ({}) in the gene {} '
                        '({})'.format(
                            exon.id, ','.join(exon.attributes['gene_id']),
                            gene_name))
            except sqlite3.IntegrityError:
                continue

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
        else:
            # If strand is unknown, put both upstream and downstream for each
            # side
            adjacent = pd.concat(adjacent_in_genome.values())
            return {UPSTREAM: adjacent, DOWNSTREAM: adjacent}

    def _junctions_genome_adjacent_to_exon(self, exon):
        """Get indices of junctions next to an exon, in genome coordinates"""
        chrom_ind = self.metadata[self.chrom] == exon.chrom

        strand_ind = self.metadata[self.strand] == exon.strand
        common_ind = chrom_ind & strand_ind

        upstream_in_genome = \
            common_ind & (self.metadata[self.exon_stop] == exon.stop)
        downstream_in_genome = \
            common_ind & (self.metadata[self.exon_start] == exon.start)
        return {UPSTREAM: upstream_in_genome, DOWNSTREAM: downstream_in_genome}

    def junctions_adjacent_to_this_exon(self, exon):
        """Get junctions adjacent to this exon

        Parameters
        ----------
        exon : gffutils.Feature
            An item in a gffutils database

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

    def upstream_downstream_exons(self):
        """Get upstream and downstream exons of each junction

        The "upstream" and "downstream" is relative to the **junction**, e.g.

            exonA   upstream      junctionX
            exonB   downstream    junctionX

        should be read as "exonA is upstream of juction X" and "exonB is
        downstream of junctionX"

        Use junctions defined in ``sj_metadata`` and exons in ``db`` to
        create triples of (exon, direction, junction), which are read like
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
        n_exons = sum(1 for _ in self.db.features_of_type(self.exon_types))

        dfs = []

        progress('Starting annotation of all junctions with known '
                 'neighboring exons ...')
        for i, exon in enumerate(self.db.features_of_type(self.exon_types)):
            if (i + 1) % 10000 == 0:
                progress('\t{}/{} exons completed'.format(i + 1, n_exons))
            df = self.junctions_adjacent_to_this_exon(exon)
            dfs.append(df)
        junction_exon_triples = pd.concat(dfs, ignore_index=True)
        done()
        return junction_exon_triples
