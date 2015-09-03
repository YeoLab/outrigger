import itertools
import sys

import numpy as np
import pandas as pd
import graphlite
from graphlite import V

graph = graphlite.connect(":memory:", graphs=['upstream', 'downstream'])

_db_doc = """db : gffutils.FeatureDB
    Database of gene, transcript, and exon features. The exons must be
    accessible by the id provided on the exon_{5,3}p_col columns. If
    not provided, certain splice types which require information about
    the transcript (AFE, ALE) cannot be annotated."""

def stringify_location(chrom, start, stop, strand, region=None):
    if region is not None:
        return '{}:{}:{}-{}:{}'.format(region, chrom, start, stop, strand)
    else:
        return '{}:{}-{}:{}'.format(chrom, start, stop, strand)

class JunctionAggregator(object):

    def __init__(self, junction_exon_triples, db=None,
                 junction_col='junction', exon_col='exon'):
        """Combine splice junctions into splicing events

        A one-line summary that does not use variable names or the
        function name.

        Several sentences providing an extended description. Refer to
        variables using back-ticks, e.g. `var`.

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

        Examples
        --------
        These are written in doctest format, and should illustrate how to
        use the function.

        >>> a=[1,2,3]
        >>> print [x + 3 for x in a]
        [4, 5, 6]
        >>> print "a\n\nb"
        a
        b
        """
        self.junction_exon_triples = junction_exon_triples
        self.db = db

        self.graph = graphlite.connect(":memory:",
                                       graphs=['upstream', 'downstream'])
        self.all_exons = junction_exon_triples[exon_col].unique()
        self.all_junctions = junction_exon_triples[junction_col].unique()

        self.items = np.concatenate([self.all_exons, self.all_junctions])
        self.int_to_item = pd.Series(self.items)
        self.item_to_int = pd.Series(
            dict((v, k) for k, v in self.int_to_item.iteritems()))

        with graph.transaction() as tr:
            for i, row in self.junction_exon_triples.iterrows():
                #         print row
                junction = self.item_to_int[row[junction_col]]
                exon = self.item_to_int[row[exon_col]]
                opposite_direction = 'upstream' \
                    if row.direction == 'downstream' else 'downstream'

                eval1 = "tr.store(V({}).{}({}))".format(exon, row.direction,
                                                        junction)
                eval2 = "tr.store(V({}).{}({}))".format(junction,
                                                        opposite_direction,
                                                        exon)
                #         print '\n', eval1
                #         print eval2
                eval(eval1)
                eval(eval2)

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

        sj_metadata['upstream'] = ''
        sj_metadata['downstream'] = ''

        print 'Starting annotation of all junctions with known exons...'
        for i, exon in enumerate(db.features_of_type('exon')):
            if (i + 1) % 10000 == 0:
                print '\t{}/{}'.format(i + 1, n_exons)
            chrom_ind = sj_metadata.chrom == exon.chrom
            strand_ind = sj_metadata.strand == exon.strand
            upstream_ind = chrom_ind & strand_ind & (
            sj_metadata.exon_stop == exon.stop)
            downstream_ind = chrom_ind & strand_ind & (
            sj_metadata.exon_start == exon.start)

            exon_id = exon.id
            if upstream_ind.any():
                if exon.strand == '+':
                    sj_metadata.loc[upstream_ind, 'upstream'] = \
                        sj_metadata.loc[upstream_ind, 'upstream'] + ',' \
                        + exon_id
                else:
                    sj_metadata.loc[upstream_ind, 'downstream'] = \
                        sj_metadata.loc[upstream_ind, 'downstream'] + ',' \
                        + exon_id

            if downstream_ind.any():
                if exon.strand == '+':
                    sj_metadata.loc[downstream_ind, 'downstream'] = \
                    sj_metadata.loc[downstream_ind, 'downstream'] + ',' + exon_id
                else:
                    sj_metadata.loc[downstream_ind, 'upstream'] = \
                    sj_metadata.loc[downstream_ind, 'upstream'] + ',' + exon_id
        print 'Done.'

        sj_metadata['upstream'] = sj_metadata['upstream'].map(
            lambda x: x.lstrip(',') if isinstance(x, str) else x)
        sj_metadata['downstream'] = sj_metadata['downstream'].map(
            lambda x: x.lstrip(',') if isinstance(x, str) else x)
        sj_metadata[['upstream', 'downstream']] = sj_metadata[
            ['upstream', 'downstream']].replace('', np.nan)
        sj_metadata[['upstream', 'downstream']].head()

        sj_metadata.loc[sj_metadata.index.map(
            lambda x: x.endswith('5p')), 'downstream'] = np.nan
        sj_metadata.loc[sj_metadata.index.map(
            lambda x: x.endswith('3p')), 'upstream'] = np.nan

    @classmethod
    def from_junction_to_exons(cls, junction_to_exons, db=None,
                               junction_col='junction',
                               upstream_col='upstream',
                               downstream_col='downstream'):
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
                                   upstream_col='upstream',
                                   downstream_col='downstream'):
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
        direction_to_exon = {'upstream': upstream_col,
                             'downstream': downstream_col}
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

    def _skipped_exon(self):
        exons_to_junctions = {}
        n_exons = self.all_exons.shape[0]

        sys.stdout.write('Trying out {} exons'
                         '...\n'.format(n_exons))
        for i, exon3_str in enumerate(self.all_exons):
            if (i + 1) % 10000 == 0:
                sys.stdout.write('\t{}/{} '
                                 'exons tested'.format(i + 1, n_exons))

            exon3 = self.item_to_int[exon3_str]
            # Get upstream junctions
            upstream_junctions = graph.find(V().upstream(exon3))

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


            for j in exon1exon2_junctions:
                exon1s = graph.find(V(j).downstream)
                exon2s = graph.find(V(j).upstream)

                exon1exon3_junctions = set(itertools.chain(
                    *[graph.find(V().downstream(x)) for x in exon1s])) & set(
                    upstream_junctions)
                exon2exon3_junctions = set(itertools.chain(
                    *[graph.find(V().downstream(x)) for x in exon2s])) & set(
                    upstream_junctions)

                exon1_confirmed = set(itertools.chain(
                    *[graph.find(V().upstream(x)) for x in
                      exon1exon3_junctions])) & set(exon1s)
                exon2_confirmed = set(itertools.chain(
                    *[graph.find(V().upstream(x)) for x in
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
                            graph.find(V().downstream(exon2))) & set(
                            upstream_junctions)
                        junctions = self.int_to_item[list(itertools.chain(
                            *[[j], exon2exon3_junction,
                              exon1exon3_junctions]))]
                        exons = self.int_to_item[[exon1, exon2, exon3]]

                        exons_to_junctions[tuple(exons)] = tuple(junctions)
                elif len(exon1_confirmed) > 1:
                    # If There's more than one exon1, use the shortest one.
                    # Since the 3' side of the exon is fixed, then the smallest
                    # size one will have the smallest upstream side
                    exon1 = min(exon1_confirmed,
                                key=lambda x: len(self.db[self.int_to_item[x]]))

                    if len(exon2_confirmed) > 1:
                        for exon2 in exon2_confirmed:
                            exon2exon3_junction = set(
                                graph.find(V().downstream(exon2))) & set(
                                upstream_junctions)
                            junctions = self.int_to_item[list(itertools.chain(
                                *[[j], exon2exon3_junction,
                                  exon1exon3_junctions]))]
                            exons = self.int_to_item[[exon1, exon2, exon3]]

                            #                     print 'junctions', junctions
                            #                     print 'exons', exons
                            exons_to_junctions[tuple(exons)] = tuple(junctions)
                    else:
                        exon2 = exon2_confirmed.pop()
                        exon2exon3_junction = set(
                            graph.find(V().downstream(exon2))) & set(
                            upstream_junctions)
                        junctions = self.int_to_item[list(itertools.chain(
                            *[[j], exon2exon3_junction,
                              exon1exon3_junctions]))]
                        exons = self.int_to_item[[exon1, exon2, exon3]]
                        exons_to_junctions[tuple(exons)] = tuple(junctions)
        return exons_to_junctions

    def _mutually_exclusive_exon(self):
        pass

    def _twin_cassette(self):
        pass

    def _alt_5p_splice_site(self):
        pass

    def _alt_3p_splice_site(self):
        pass

    def _alt_first_exon(self):
        pass

    def _alt_last_exon(self):
        pass
