import itertools
import sys

import graphlite
from graphlite import V
import numpy as np
import pandas as pd

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
        return '{0}:{1}:{2}-{3}:{4}'.format(region, chrom, start, stop, strand)
    else:
        return '{0}:{1}-{2}:{3}'.format(chrom, start, stop, strand)


def opposite(direction):
    return UPSTREAM if direction == DOWNSTREAM else DOWNSTREAM


def make_junction_direction_df(direction_ind, direction, exon_id):
    return pd.DataFrame(zip(itertools.cycle((exon_id,)),
                     itertools.cycle((direction,)),
                     direction_ind[direction_ind].index),
                 columns=['exon', 'direction', 'junction'])

def make_junction_exon_triples(sj_metadata, db):
    """Get upstream and downstream exons in db, of introns in sj_metadata

    :param sj_metadata:
    :type sj_metadata:
    :param db:
    :type db:
    :return:
    :rtype:
    """
    if 'exon_start' not in sj_metadata:
        sj_metadata['exon_start'] = sj_metadata['intron_stop'] + 1
    if 'exon_stop' not in sj_metadata:
        sj_metadata['exon_stop'] = sj_metadata['intron_start'] - 1

    n_exons = sum(1 for _ in db.features_of_type('exon'))

    # sj_metadata['upstream_exon'] = ''
    # sj_metadata['downstream_exon'] = ''

    dfs = []

    sys.stdout.write('Starting annotation of all junctions with known '
                     'exons...\n')
    for i, exon in enumerate(db.features_of_type('exon')):
        if (i + 1) % 10000 == 0:
            sys.stdout.write('\t{}/{} exons completed\n'.format(i + 1, n_exons))
        chrom_ind = sj_metadata.chrom == exon.chrom
        strand_ind = sj_metadata.strand == exon.strand
        upstream_ind = chrom_ind & strand_ind & \
                       (sj_metadata.exon_stop == exon.stop)
        downstream_ind = chrom_ind & strand_ind & \
                         (sj_metadata.exon_start == exon.start)

        exon_id = exon.id
        if upstream_ind.any():
            upstream_df = make_junction_direction_df(upstream_ind, UPSTREAM,
                                                     exon_id)
            dfs.append(upstream_df)
            # if exon.strand == '+':
            #     sj_metadata.loc[upstream_ind, 'upstream_exon'] = \
            #     sj_metadata.loc[upstream_ind, 'upstream_exon'] + ',' + exon_id
            # else:
            #     sj_metadata.loc[upstream_ind, 'downstream_exon'] = \
            #     sj_metadata.loc[
            #         upstream_ind, 'downstream_exon'] + ',' + exon_id

        if downstream_ind.any():
            downstream_df = make_junction_direction_df(downstream_ind,
                                                       DOWNSTREAM, exon_id)
            dfs.append(downstream_df)
            # if exon.strand == '+':
            #     sj_metadata.loc[downstream_ind, 'downstream_exon'] = \
            #     sj_metadata.loc[
            #         downstream_ind, 'downstream_exon'] + ',' + exon_id
            # else:
            #     sj_metadata.loc[downstream_ind, 'upstream_exon'] = \
            #     sj_metadata.loc[
            #         downstream_ind, 'upstream_exon'] + ',' + exon_id
    junction_exon_triples = pd.concat(dfs, ignore_index=True)
    sys.stdout.write('Done.\n')
    return junction_exon_triples


class JunctionAggregator(object):

    def __init__(self, junction_exon_triples, db=None, junction_col='junction',
                 exon_col='exon', debug=False):
        """Combine splice junctions into splicing events

        Parameters
        ----------
        junction_exon_triples : pandas.DataFrame
            of "exon, direction, junction", e.g.:
            exon1, upstream, junction12

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

    def event_dict_to_df(self, events, exon_names, junction_names):
        columns = list(exon_names) + list(junction_names) + ['event_id']
        data = pd.DataFrame(index=np.arange(len(events)), columns=columns)
        for i, (exons, junctions) in enumerate(events.items()):
            event_id = '@'.join(exons)
            data.loc[i, exon_names] = list(exons)
            data.loc[i, junction_names] = list(junctions)
            data.loc[i, 'event_id'] = event_id
        return data

    def skipped_exon(self):
        events = {}
        n_exons = self.exons.shape[0]

        sys.stdout.write('Trying out {0} exons'
                         '...\n'.format(n_exons))
        for i, exon1 in enumerate(self.item_to_region.values):
            if (i + 1) % 10000 == 0:
                sys.stdout.write('\t{0}/{1} '
                                 'exons tested'.format(i + 1, n_exons))

            exon1_i = self.item_to_int[exon1.name]
            exon23s = list(
                self.graph.find(
                    V().downstream(exon1_i)).traverse(V().upstream))
            exon23s = self.item_to_region[self.int_to_item[exon23s]]

            for exon_a, exon_b in itertools.combinations(exon23s, 2):
                if not exon_a.overlaps(exon_b):
                    exon2 = min((exon_a, exon_b), key=lambda x: x.start)
                    exon3 = max((exon_a, exon_b), key=lambda x: x.start)

                    exon2_i = self.item_to_int[exon2.name]
                    exon3_i = self.item_to_int[exon3.name]

                    exon23_junction = self.graph.find(
                        V(exon2_i).upstream).intersection(
                        V().upstream(exon3_i))
                    exon23_junction = self.int_to_item[set(exon23_junction)]
                    if not exon23_junction.empty:
                        # Isoform 1 - corresponds to Psi=0. Exclusion of exon2
                        exon13_junction = self.graph.find(
                            V(exon1_i).upstream) \
                            .intersection(V(exon3_i).downstream)

                        # Isoform 2 - corresponds to Psi=1. Inclusion of exon2
                        exon12_junction = self.graph.find(
                            V(exon1_i).upstream) \
                            .intersection(V(exon2_i).downstream)
                        exon23_junction = self.graph.find(
                            V(exon2_i).upstream) \
                            .intersection(V(exon3_i).downstream)

                        junctions = list(itertools.chain(
                            *[exon12_junction, exon23_junction,
                              exon13_junction]))
                        junctions = self.int_to_item[junctions].tolist()
                        exons = exon1.name, exon2.name, exon3.name

                        events[exons] = junctions
        events = self.event_dict_to_df(events,
                                       exon_names=['exon1', 'exon2', 'exon3'],
                                       junction_names=['junction12',
                                                       'junction23',
                                                       'junction13'])
        return events

    def mutually_exclusive_exon(self):
        events = {}

        for exon1_name in self.exons:
            exon1_i = self.item_to_int[exon1_name]

            downstream_junctions = set(self.graph.find(
                V().downstream(exon1_i)))
            downstream_junctions
            exon23s_from1 = list(
                self.graph.find(V().downstream(
                    exon1_i)).traverse(V().upstream))
            exon4s = self.graph.find(V().downstream(exon1_i)).traverse(
                V().upstream).traverse(V().upstream).traverse(V().upstream)
            exon23s_from4 = exon4s.traverse(V().downstream).traverse(
                V().downstream)

            exon23s = self.int_to_item[set(exon23s_from4) & set(exon23s_from1)]

            exon23s = self.item_to_region[exon23s]

            for exon_a, exon_b in itertools.combinations(exon23s, 2):
                if not exon_a.overlaps(exon_b):
                    exon2 = min((exon_a, exon_b), key=lambda x: x.start)
                    exon3 = max((exon_a, exon_b), key=lambda x: x.start)

                    exon2_i = self.item_to_int[exon2.name]
                    exon3_i = self.item_to_int[exon3.name]

                    exon4_from2 = set(
                        self.graph.find(V(exon2_i).upstream).traverse(
                            V().upstream))
                    exon4_from3 = set(
                        self.graph.find(V(exon3_i).upstream).traverse(
                            V().upstream))
                    try:
                        exon4_i = (exon4_from2 & exon4_from3).pop()
                        exon4_name = self.int_to_item[exon4_i]
                        # Isoform 1 - corresponds to Psi=0. Inclusion of exon3
                        exon13_junction = self.graph.find(
                            V(exon1_i).upstream).intersection(
                            V(exon3_i).downstream)
                        exon34_junction = self.graph.find(
                            V(exon3_i).upstream) \
                            .intersection(V(exon4_i).downstream)

                        # Isoform 2 - corresponds to Psi=1. Inclusion of exon2
                        exon12_junction = self.graph.find(
                            V(exon1_i).upstream).intersection(
                            V(exon2_i).downstream)
                        exon24_junction = self.graph.find(
                            V(exon2_i).upstream) \
                            .intersection(V(exon4_i).downstream)

                        exon_tuple = exon1_name, exon2.name, exon3.name, \
                            exon4_name
                        #             print exon12_junction.next()
                        junctions = list(
                            itertools.chain(*[exon13_junction, exon34_junction,
                                              exon12_junction,
                                              exon24_junction]))
                        junctions = self.int_to_item[junctions].tolist()

                        events[exon_tuple] = junctions
                    except:
                        pass
        events = self.event_dict_to_df(events,
                                       exon_names=['exon1', 'exon2', 'exon3',
                                                   'exon4'],
                                       junction_names=['junction13',
                                                       'junction34',
                                                       'junction12',
                                                       'junction24'])
        return events

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

import itertools

BEST_TAGS = 'appris_principal', 'appris_candidate', 'CCDS', 'basic'

transcript_cols = ['isoform1_transcripts', 'isoform2_transcripts']

def get_attribute(features, attribute):
    try:
        for feature in features:
            try:
                yield feature[attribute]
            except KeyError:
                pass
    except TypeError:
        # The features aren't iterable
        pass

def get_feature_attribute_with_value(features, attribute, value):
    try:
        for feature in features:
            try:
                if value in feature[attribute]:
                    yield feature.id
            except KeyError:
                pass
    except TypeError:
        # The features aren't iterable
        pass

def get_feature_attribute_startswith_value(features, attribute, value):
    try:
        for feature in features:
            try:
                if any(map(lambda x: x.startswith(value),
                           feature[attribute])):
                    yield feature.id
            except KeyError:
                pass
    except TypeError:
        # The features aren't iterable
        pass

def consolidate_junction_events(df, db, event_col='event_id',
                                transcript_cols=['isoform1_transcripts',
                                                 'isoform2_transcripts']):
    if len(df) == 1:
        return 'only one', df[event_col].values[0]

    df_isoforms = df[transcript_cols].applymap(
        lambda x: np.nan if len(x) == 0 else map(lambda y: db[y], x))
    df_isoforms = df_isoforms.dropna(how='all')

    if df_isoforms.empty:
        return 'random,no gencode transcripts', df.loc[
            np.random.choice(df.index), event_col]

    if len(df_isoforms) == 1:
        return 'one event with gencode transcripts', df.loc[
            df_isoforms.index[0], event_col]

    df_tags = df_isoforms.applymap(
        lambda x: tuple(
            itertools.chain(*get_attribute(x, 'tag')))
        if not isinstance(x, float) else x)

    df_tags = df_tags.applymap(
        lambda x: x if not isinstance(x, list) or len(x) > 0 else np.nan)
    df_tags = df_tags.dropna(how='all')
    if df_tags.empty:
        return 'random df_isoforms', df_isoforms.loc[
            np.random.choice(df_isoforms.index)]

    for tag in BEST_TAGS:
        df_this_tag = df_tags.applymap(
            lambda x: map(lambda y: y.startswith(tag), x)
            if isinstance(x, tuple) else False)

        # Which isoform has at least one true
        df_this_tag = df_this_tag.any(axis=1)
        #         print df_this_tag
        if df_this_tag.any():
            best_index = np.random.choice(df_this_tag.index[df_this_tag])
            #             print '- best isoform:', tag, best_index
            #             print df.loc[best_index].event_id
            return 'best,{}'.format(tag), df.loc[best_index].event_id
    else:
        return 'random,no good tags', df.loc[np.random.choice(df.index),
                                             event_col]