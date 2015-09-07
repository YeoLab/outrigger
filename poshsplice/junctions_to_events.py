import itertools
import sys

import numpy as np
import pandas as pd
import graphlite
from graphlite import V

from .region import Region

_db_doc = """db : gffutils.FeatureDB
    Database of gene, transcript, and exon features. The exons must be
    accessible by the id provided on the exon_{5,3}p_col columns. If
    not provided, certain splice types which require information about
    the transcript (AFE, ALE) cannot be annotated."""

UPSTREAM = 'upstream'
DOWNSTREAM = 'downstream'
DIRECTIONS = [UPSTREAM, DOWNSTREAM]


def stringify_location(chrom, start, stop, strand, region=None):
    if region is not None:
        return '{}:{}:{}-{}:{}'.format(region, chrom, start, stop, strand)
    else:
        return '{}:{}-{}:{}'.format(chrom, start, stop, strand)


def opposite(direction):
    return UPSTREAM if direction == DOWNSTREAM else DOWNSTREAM


class JunctionAggregator(object):

    def __init__(self, junction_exon_triples, db=None, junction_col='junction',
                 exon_col='exon', debug=False):
        """Combine splice junctions into splicing events

        Parameters
        ----------
        junction_exon_triples : pandas.DataFrame

        db : gffutils.FeatureDB
            Gffutils Database of gene, transcript, and exon features. The exons
            must be accessible by the id provided on the `exon_col`
            columns. If not provided, certain splice types which require
            information about the transcript (AFE, ALE) cannot be annotated.

        Returns
        -------
        type
            Explanation of anonymous return value of type ``type``.
        describe : type
            Explanation of return value named `describe`.
        out : type
            Explanation of `out`.

        Other Parameters
        ----------------
        only_seldom_used_keywords : type
            Explanation
        common_parameters_listed_above : type
            Explanation

        Raises
        ------
        BadException
            Because you shouldn't have done that.

        See Also
        --------
        otherfunc : relationship (optional)
        newfunc : Relationship (optional), which could be fairly long, in which
                  case the line wraps here.
        thirdfunc, fourthfunc, fifthfunc

        Notes
        -----
        Notes about the implementation algorithm (if needed).

        This can have multiple paragraphs.

        You may include some math:

        .. math:: X(e^{j\omega } ) = x(n)e^{ - j\omega n}

        And even use a greek symbol like :math:`omega` inline.

        References
        ----------
        Cite the relevant literature, e.g. [1]_.  You may also cite these
        references in the notes section above.

        .. [1] O. McNoleg, "The integration of GIS, remote sensing,
           expert systems and adaptive co-kriging for environmental habitat
           modelling of the Highland Haggis using object-oriented, fuzzy-logic
           and neural-network techniques," Computers & Geosciences, vol. 22,
           pp. 585-588, 1996.
        """
        self.junction_exon_triples = junction_exon_triples
        self.db = db
        self.debug = debug

        self.graph = graphlite.connect(":memory:", graphs=DIRECTIONS)
        self.exons = junction_exon_triples[exon_col].unique()
        self.junctions = junction_exon_triples[junction_col].unique()

        self.items = np.concatenate([self.exons, self.junctions])
        self.int_to_item = pd.Series(self.items)
        self.item_to_int = pd.Series(
            dict((v, k) for k, v in self.int_to_item.iteritems()))
        self.item_to_region = pd.Series(map(Region, self.items),
                                        index=self.items)

        with self.graph.transaction() as tr:
            for i, row in self.junction_exon_triples.iterrows():
                junction = row[junction_col]
                exon = row[exon_col]

                junction_i = self.item_to_int[junction]
                exon_i = self.item_to_int[exon]

                if self.debug:
                    sys.stdout.write('\n{} is {} of {}\n'.format(
                        exon, row.direction, junction))
                    sys.stdout.write('{} is {} of {}\n'.format(
                        junction, opposite(row.direction), exon))

                tr.store(getattr(V(exon_i), row.direction)(junction_i))
                tr.store(getattr(V(junction_i),
                                 opposite(row.direction))(exon_i))

    @classmethod
    def from_sj(cls, sj_metadata, db=None):
        """Annotates junctions with nearby exons and initializes Annotator

        Parameters
        ----------
        var1 : array_like
            Array_like means all those objects -- lists, nested lists, etc. --
            that can be converted to an array.  We can also refer to
            variables like `var1`.
        db : gffutils.FeatureDB
            Database of gene, transcript, and exon features. The exons must be
            accessible by the id provided on the exon_{5,3}p_col columns. If
            not provided, certain splice types which require information about
            the transcript (AFE, ALE) cannot be annotated.

        Returns
        -------
        type
            Explanation of anonymous return value of type ``type``.
        describe : type
            Explanation of return value named `describe`.
        out : type
            Explanation of `out`.
        """

        sj_metadata[UPSTREAM] = ''
        sj_metadata[DOWNSTREAM] = ''

        n_exons = sum(1 for _ in db.features_of_type('exon'))

        sys.stdout.write('Starting annotation of all junctions with known '
                         'exons...\n')
        for i, exon in enumerate(db.features_of_type('exon')):
            if (i + 1) % 10000 == 0:
                sys.stdout.write('\t{}/{}\n'.format(i + 1, n_exons))
            chrom_ind = sj_metadata.chrom == exon.chrom
            strand_ind = sj_metadata.strand == exon.strand
            upstream_ind = chrom_ind & strand_ind & (
                sj_metadata.exon_stop == exon.stop)
            downstream_ind = chrom_ind & strand_ind & (
                sj_metadata.exon_start == exon.start)

            exon_id = exon.id
            if upstream_ind.any():
                if exon.strand == '+':
                    sj_metadata.loc[upstream_ind, UPSTREAM] = \
                        sj_metadata.loc[upstream_ind, UPSTREAM] + ',' \
                        + exon_id
                else:
                    sj_metadata.loc[upstream_ind, DOWNSTREAM] = \
                        sj_metadata.loc[upstream_ind, DOWNSTREAM] + ',' \
                        + exon_id

            if downstream_ind.any():
                if exon.strand == '+':
                    sj_metadata.loc[downstream_ind, DOWNSTREAM] = \
                        sj_metadata.loc[downstream_ind, DOWNSTREAM] + ',' \
                        + exon_id
                else:
                    sj_metadata.loc[downstream_ind, UPSTREAM] = \
                        sj_metadata.loc[downstream_ind, UPSTREAM] + ',' \
                        + exon_id
        sys.stdout.write('Done.\n')

        sj_metadata[UPSTREAM] = sj_metadata[UPSTREAM].map(
            lambda x: x.lstrip(',') if isinstance(x, str) else x)
        sj_metadata[DOWNSTREAM] = sj_metadata[DOWNSTREAM].map(
            lambda x: x.lstrip(',') if isinstance(x, str) else x)
        sj_metadata[DIRECTIONS] = sj_metadata[DIRECTIONS].replace('', np.nan)

        sj_metadata.loc[sj_metadata.index.map(
            lambda x: x.endswith('5p')), DOWNSTREAM] = np.nan
        sj_metadata.loc[sj_metadata.index.map(
            lambda x: x.endswith('3p')), UPSTREAM] = np.nan

    @classmethod
    def from_junction_to_exons(cls, junction_to_exons, db=None,
                               junction_col='junction', upstream_col=UPSTREAM,
                               downstream_col=DOWNSTREAM):
        """Initialize Annotator from table with junctions and nearby exons

        Parameters
        ----------
        sj_metadata : pandas.DataFrame
            A table with a column indicating "junction_location"
        db : gffutils.FeatureDB
            Database of gene, transcript, and exon features. The exons must be
            accessible by the id provided on the exon_{5,3}p_col columns. If
            not provided, certain splice types which require information about
            the transcript (AFE, ALE) cannot be annotated.
        junction_col : str
            Column name of the raw junction location (without |5p or |3p
            annotated)
        upstream_col : str
            Column name where exons upstream of the junction are stored
        downstream_col : str
            Column name where exons downstream of the junction are stored

        Returns
        -------
        junction_exon_triples
            A three-column table of junction_location, exon, and direction

        """
        junction_exons = cls.make_junction_exon_triples(
            junction_to_exons, junction_col=junction_col,
            upstream_col=upstream_col, downstream_col=downstream_col)
        return cls(junction_exons, db=db)

    @staticmethod
    def make_junction_exon_triples(junction_to_exons,
                                   junction_col='junction',
                                   upstream_col=UPSTREAM,
                                   downstream_col=DOWNSTREAM):
        """Create tidy table of exons upstream and downstream of a junction

        Parameters
        ----------
        sj_metadata : pandas.DataFrame
            A table with a column indicating "junction_location"
        junction_col : str
            Column name of the raw junction location (without |5p or |3p
            annotated)
        upstream_col : str
            Column name where exons upstream of the junction are stored
        downstream_col : str
            Column name where exons downstream of the junction are stored

        Returns
        -------
        junction_exon_triples
            A three-column table of junction_location, exon, and direction

        Examples
        --------
        >>> import pandas as pd
        >>> sj_metadata = pd.DataFrame(
        {'junction':['chr1:201-299:+', 'chr1:401:499:+'],
         'upstream': ['exon:chr1:100-200:+,exon:chr1:50-200:+',
                     'exon:chr1:300-400:+'],
         'downstream':['exon:chr1:300-400:+',
         'exon:chr1:500-600:+,exon:chr1:500-650:+']})
        >>> Annotator.make_junction_exon_triples(sj_metadata)

        """
        grouped = junction_to_exons.groupby(junction_col)
        direction_to_exon = {UPSTREAM: upstream_col,
                             DOWNSTREAM: downstream_col}
        dfs = []
        for direction, exon in direction_to_exon.items():
            df = grouped.apply(
                lambda x: x[exon].dropna().str.split(',').apply(pd.Series, 1))
            df = df.stack()
            df.index = df.index.droplevel((-2, -1))
            df = df.reset_index()
            df = df.rename(columns={0: 'exon'})
            df['direction'] = direction
            dfs.append(df)
        junction_exons = pd.concat(dfs, ignore_index=True)
        return junction_exons

    def skipped_exon(self):
        exons_to_junctions = {}
        n_exons = self.exons.shape[0]

        sys.stdout.write('Trying out {} exons'
                         '...\n'.format(n_exons))
        for i, exon3_str in enumerate(self.exons):
            if (i + 1) % 10000 == 0:
                sys.stdout.write('\t{}/{} '
                                 'exons tested'.format(i + 1, n_exons))

            exon3 = self.item_to_int[exon3_str]
            # Get upstream junctions
            upstream_junctions = self.graph.find(V().upstream(exon3))

            # Get upstream exons of this exon (which exons have these upstream
            # junctions, downstream of the exon)
            upstream_exons = upstream_junctions.traverse(V().downstream)

            # From those upstream exons, get both upstream and downstream
            # junctions
            upstream_exons_upstream_junctions = set(
                upstream_exons.traverse(V().upstream))
            upstream_exons_downstream_junctions = set(
                upstream_exons.traverse(V().downstream))

            # Of those junctions, which overlap --> SE event
            exon1exon2_junctions = upstream_exons_downstream_junctions \
                & upstream_exons_upstream_junctions

            upstream_junctions = set(upstream_junctions)

            for j in exon1exon2_junctions:
                exon1s = self.graph.find(V(j).downstream)
                exon2s = self.graph.find(V(j).upstream)

                exon1exon3_junctions = set(itertools.chain(
                    *[self.graph.find(V().downstream(x)) for x in exon1s])) \
                    & upstream_junctions
                exon2exon3_junctions = set(itertools.chain(
                    *[self.graph.find(V().downstream(x)) for x in exon2s])) \
                    & upstream_junctions

                exon1_confirmed = set(itertools.chain(
                    *[self.graph.find(V().upstream(x)) for x in
                      exon1exon3_junctions])) & set(exon1s)
                exon2_confirmed = set(itertools.chain(
                    *[self.graph.find(V().upstream(x)) for x in
                      exon2exon3_junctions])) & set(exon2s)

                if len(exon1_confirmed) == 1 and len(exon2_confirmed) == 1:
                    junctions = self.int_to_item[list(itertools.chain(
                        *[[j], exon2exon3_junctions, exon1exon3_junctions]))]
                    exons = self.int_to_item[
                        [exon1_confirmed.pop(), exon2_confirmed.pop(), exon3]]
                    exons_to_junctions[tuple(exons)] = tuple(junctions)

                if len(exon1_confirmed) == 1 and len(exon2_confirmed) > 1:
                    exon1 = exon1_confirmed.pop()

                    for exon2 in exon2_confirmed:
                        exon2exon3_junction = set(
                            self.graph.find(V().downstream(exon2))) \
                            & upstream_junctions
                        junctions = self.int_to_item[list(itertools.chain(
                            *[[j], exon2exon3_junction,
                              exon1exon3_junctions]))]
                        exons = self.int_to_item[[exon1, exon2, exon3]]

                        exons_to_junctions[tuple(exons)] = tuple(junctions)
                elif len(exon1_confirmed) > 1:
                    # If There's more than one exon1, use the shortest one.
                    # Since the 3' side of the exon is fixed, then the smallest
                    # size one will have the largest start number
                    exon1 = max(exon1_confirmed,
                                key=lambda x:
                                self.item_to_region[self.int_to_item[x]].start)

                    if len(exon2_confirmed) > 1:
                        for exon2 in exon2_confirmed:
                            exon2exon3_junction = set(
                                self.graph.find(V().downstream(exon2))) \
                                & upstream_junctions
                            junctions = self.int_to_item[list(itertools.chain(
                                *[[j], exon2exon3_junction,
                                  exon1exon3_junctions]))]
                            exons = self.int_to_item[[exon1, exon2, exon3]]

                            exons_to_junctions[tuple(exons)] = tuple(junctions)
                    else:
                        exon2 = exon2_confirmed.pop()
                        exon2exon3_junction = set(
                            self.graph.find(V().downstream(exon2))) \
                            & upstream_junctions
                        junctions = self.int_to_item[list(itertools.chain(
                            *[[j], exon2exon3_junction,
                              exon1exon3_junctions]))]
                        exons = self.int_to_item[[exon1, exon2, exon3]]
                        exons_to_junctions[tuple(exons)] = tuple(junctions)
        return exons_to_junctions

    def mutually_exclusive_exon(self):
        events_to_junctions = {}

        for exon1_name in self.exons:
            exon1 = self.item_to_int[exon1_name]

            downstream_junctions = set(self.graph.find(V().downstream(exon1)))
            downstream_junctions
            exon23s_from1 = list(
                self.graph.find(V().downstream(exon1)).traverse(V().upstream))
            exon4s = self.graph.find(V().downstream(exon1)).traverse(
                V().upstream).traverse(V().upstream).traverse(V().upstream)
            exon23s_from4 = exon4s.traverse(V().downstream).traverse(
                V().downstream)

            exon23s = self.int_to_item[set(exon23s_from4) & set(exon23s_from1)]

            exon23s = self.item_to_region[exon23s]

            for exon_a, exon_b in itertools.combinations(exon23s, 2):
                print exon_a.name, exon_b.name
                if not exon_a.overlaps(exon_b):
                    exon2 = min((exon_a, exon_b), key=lambda x: x.start)
                    exon3 = max((exon_a, exon_b), key=lambda x: x.start)
                    exon4_from2 = set(
                        self.graph.find(V(self.item_to_int[exon2.name]).upstream).traverse(
                            V().upstream))
                    exon4_from3 = set(
                        self.graph.find(V(self.item_to_int[exon3.name]).upstream).traverse(
                            V().upstream))
                    exon4 = exon4_from2 & exon4_from3
                    exon4_name = self.int_to_item[exon4]
                    if not exon4_name.empty:
                        exon4_name = exon4_name.values[0]
                        # Isoform 1 - corresponds to Psi=0. Inclusion of exon3
                        exon13_junction = self.graph.find(V(exon1).upstream) \
                            .intersection(V(self.item_to_int[exon3.name]).downstream)
                        exon34_junction = self.graph.find(
                            V(self.item_to_int[exon3.name]).upstream) \
                            .intersection(V(self.item_to_int[exon4_name]).downstream)

                        # Isoform 2 - corresponds to Psi=1. Inclusion of exon2
                        exon12_junction = self.graph.find(V(exon1).upstream) \
                            .intersection(V(self.item_to_int[exon2.name]).downstream)
                        exon24_junction = self.graph.find(
                            V(self.item_to_int[exon2.name]).upstream) \
                            .intersection(V(self.item_to_int[exon4_name]).downstream)

                        exon_tuple = exon1_name, exon2.name, exon3.name, \
                                     exon4_name
                        #             print exon12_junction.next()
                        junctions = list(
                            itertools.chain(*[exon13_junction, exon34_junction,
                                              exon12_junction,
                                              exon24_junction]))
                        junctions = self.int_to_item[junctions].tolist()

                        events_to_junctions[exon_tuple] = junctions
        return events_to_junctions

    def twin_cassette(self):
        pass

    def alt_5p_splice_site(self):
        pass

    def alt_3p_splice_site(self):
        pass

    def alt_first_exon(self):
        pass

    def alt_last_exon(self):
        pass
