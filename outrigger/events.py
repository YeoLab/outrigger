import itertools
import logging
import sys

import graphlite
from graphlite import V
import numpy as np
import pandas as pd

from .region import Region
from .junctions import UPSTREAM, DOWNSTREAM, DIRECTIONS


def stringify_location(chrom, start, stop, strand, region=None):
    """"""
    if region is not None:
        return '{0}:{1}:{2}-{3}:{4}'.format(region, chrom, start, stop,
                                            strand)
    else:
        return '{0}:{1}-{2}:{3}'.format(chrom, start, stop, strand)


def opposite(direction):
    return UPSTREAM if direction == DOWNSTREAM else DOWNSTREAM


class EventMaker(object):

    def __init__(self, junction_exon_triples, db=None, junction_col='junction',
                 exon_col='exon'):
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
        """
        self.log = logging.getLogger('EventMaker')
        self.junction_exon_triples = junction_exon_triples
        self.db = db

        self.graph = graphlite.connect(":memory:", graphs=DIRECTIONS)
        self.exons = junction_exon_triples[exon_col].unique()
        self.junctions = junction_exon_triples[junction_col].unique()

        self.items = tuple(np.concatenate([self.exons, self.junctions]))
        self.item_to_region = pd.Series(map(Region, self.items),
                                        index=self.items)

        with self.graph.transaction() as tr:
            for i, row in self.junction_exon_triples.iterrows():
                junction = row[junction_col]
                exon = row[exon_col]

                junction_i = self.items.index(junction)
                exon_i = self.items.index(exon)

                self.log.debug('\n{} is {} of {}\n'.format(
                    exon, row.direction, junction))
                self.log.debug('{} is {} of {}\n'.format(
                    junction, opposite(row.direction), exon))

                tr.store(getattr(V(exon_i), row.direction)(junction_i))
                tr.store(getattr(V(junction_i),
                                 opposite(row.direction))(exon_i))

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
        >>> EventMaker.get_adjacent_exons(sj_metadata)

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

    def _check_exon_in_se_event(self, exon1_name):
        event = {}

        exon1_i = self.items.index(exon1_name)
        exon23s = list(
            self.graph.find(
                V().downstream(exon1_i)).traverse(V().upstream))
        exon23s = self.item_to_region[[self.items[i] for i in exon23s]]

        for exon_a, exon_b in itertools.combinations(exon23s, 2):
            if not exon_a.overlaps(exon_b):
                exon2 = min((exon_a, exon_b), key=lambda x: x.start)
                exon3 = max((exon_a, exon_b), key=lambda x: x.start)

                exon2_i = self.items.index(exon2.name)
                exon3_i = self.items.index(exon3.name)

                exon23_junction_i = self.graph.find(
                    V(exon2_i).upstream).intersection(
                    V().upstream(exon3_i))
                exon23_junction = [self.items[i] for i in
                                   set(exon23_junction_i)]
                if len(exon23_junction) > 0:
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

                    junctions_i = list(itertools.chain(
                        *[exon12_junction, exon23_junction,
                          exon13_junction]))
                    junctions = [self.items[i] for i in junctions_i]
                    exons = exon1_name, exon2.name, exon3.name

                    event[exons] = junctions
        return event

    def skipped_exon(self):
        events = {}
        n_exons = self.exons.shape[0]

        sys.stdout.write('Trying out {0} exons'
                         '...\n'.format(n_exons))
        for i, exon1_name in enumerate(self.exons):
            if (i + 1) % 10000 == 0:
                sys.stdout.write('\t{0}/{1} '
                                 'exons tested'.format(i + 1, n_exons))

            exon1_i = self.items.index(exon1_name)
            exon23s = list(
                self.graph.find(
                    V().downstream(exon1_i)).traverse(V().upstream))
            exon23s = self.item_to_region[[self.items[i] for i in exon23s]]

            for exon_a, exon_b in itertools.combinations(exon23s, 2):
                if not exon_a.overlaps(exon_b):
                    exon2 = min((exon_a, exon_b), key=lambda x: x.start)
                    exon3 = max((exon_a, exon_b), key=lambda x: x.start)

                    exon2_i = self.items.index(exon2.name)
                    exon3_i = self.items.index(exon3.name)

                    exon23_junction_i = self.graph.find(
                        V(exon2_i).upstream).intersection(
                        V().upstream(exon3_i))
                    exon23_junction = [self.items[i] for i in
                                       set(exon23_junction_i)]
                    if len(exon23_junction) > 0:
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

                        junctions_i = list(itertools.chain(
                            *[exon12_junction, exon23_junction,
                              exon13_junction]))
                        junctions = [self.items[i] for i in junctions_i]
                        exons = exon1_name, exon2.name, exon3.name

                        events[exons] = junctions
        events = self.event_dict_to_df(events,
                                       exon_names=['exon1', 'exon2', 'exon3'],
                                       junction_names=['junction12',
                                                       'junction23',
                                                       'junction13'])
        return events

    def _check_exon_in_mxe_event(self, exon1_name):
        event = {}
        exon1_i = self.items.index(exon1_name)

        exon23s_from1 = list(
            self.graph.find(V().downstream(
                exon1_i)).traverse(V().upstream))
        exon4s = self.graph.find(V().downstream(exon1_i)).traverse(
            V().upstream).traverse(V().upstream).traverse(V().upstream)
        exon23s_from4 = exon4s.traverse(V().downstream).traverse(
            V().downstream)

        exon23s = set(exon23s_from4) & set(exon23s_from1)
        exon23s = [self.items[i] for i in exon23s]

        exon23s = self.item_to_region[exon23s]

        for exon_a, exon_b in itertools.combinations(exon23s, 2):
            if not exon_a.overlaps(exon_b):
                exon2 = min((exon_a, exon_b), key=lambda x: x.start)
                exon3 = max((exon_a, exon_b), key=lambda x: x.start)

                exon2_i = self.items.index(exon2.name)
                exon3_i = self.items.index(exon3.name)

                exon4_from2 = set(
                    self.graph.find(V(exon2_i).upstream).traverse(
                        V().upstream))
                exon4_from3 = set(
                    self.graph.find(V(exon3_i).upstream).traverse(
                        V().upstream))
                try:
                    exon4_i = (exon4_from2 & exon4_from3).pop()
                    exon4_name = self.items[exon4_i]
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
                    junctions = [self.items[i] for i in junctions]

                    event[exon_tuple] = junctions
                except:
                    pass
        return event

    def mutually_exclusive_exon(self):
        events = {}

        for exon1_name in self.exons:
            exon1_i = self.items.index(exon1_name)

            exon23s_from1 = list(
                self.graph.find(V().downstream(
                    exon1_i)).traverse(V().upstream))
            exon4s = self.graph.find(V().downstream(exon1_i)).traverse(
                V().upstream).traverse(V().upstream).traverse(V().upstream)
            exon23s_from4 = exon4s.traverse(V().downstream).traverse(
                V().downstream)

            exon23s = set(exon23s_from4) & set(exon23s_from1)
            exon23s = [self.items[i] for i in exon23s]

            exon23s = self.item_to_region[exon23s]

            for exon_a, exon_b in itertools.combinations(exon23s, 2):
                if not exon_a.overlaps(exon_b):
                    exon2 = min((exon_a, exon_b), key=lambda x: x.start)
                    exon3 = max((exon_a, exon_b), key=lambda x: x.start)

                    exon2_i = self.items.index(exon2.name)
                    exon3_i = self.items.index(exon3.name)

                    exon4_from2 = set(
                        self.graph.find(V(exon2_i).upstream).traverse(
                            V().upstream))
                    exon4_from3 = set(
                        self.graph.find(V(exon3_i).upstream).traverse(
                            V().upstream))
                    try:
                        exon4_i = (exon4_from2 & exon4_from3).pop()
                        exon4_name = self.items[exon4_i]
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
                            itertools.chain(*[exon13_junction,
                                              exon34_junction,
                                              exon12_junction,
                                              exon24_junction]))
                        junctions = [self.items[i] for i in junctions]

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


class EventConsolidator(object):

    def __init__(self, events, db, isoform1_exons, isoform2_exons, junctions,
                 best_tags):
        """Group events that share junctions and pick the best one

        For events that look identical based on their junctions, pick the best
        event based on biological information from the genome annotation
        """
        self.events = events
        self.db = db
        self.isoform1_exons = isoform1_exons
        self.isoform2_exons = isoform2_exons
        self.junctions = junctions
        self.best_tags = best_tags

        self.transcript_cols = ['isoform1_transcripts', 'isoform2_transcripts']

        # Add isoform columns
        self.events['isoform1_transcripts'] = self.events.apply(
            lambda row: map(
                lambda x: x.id, self._get_isoform_transcripts(
                    row, exons=self.isoform1_exons, exclude_exons='exon2',
                    db=self.db)), axis=1)
        self.events['isoform2_transcripts'] = self.events.apply(
            lambda row: map(
                lambda x: x.id, self._get_isoform_transcripts(
                    row, exons=self.isoform2_exons, db=self.db)), axis=1)

    @staticmethod
    def _get_isoform_transcripts(row, exons, db, exclude_exons=None,
                                 featuretype='transcript'):

        transcripts = map(
            lambda x: set(db.parents(db[row[x]],
                                     featuretype=featuretype)), exons)
        transcripts = set.intersection(*transcripts)

        if exclude_exons is not None:
            exclude_exons = [exclude_exons] if isinstance(exclude_exons, str) \
                else exclude_exons
            exclude_transcripts = map(
                lambda x: set(db.parents(db[row[x]], featuretype=featuretype)),
                exclude_exons)
            transcripts = transcripts.difference(
                set.intersection(*exclude_transcripts))
        return transcripts

    @staticmethod
    def _get_attribute(features, attribute):
        try:
            for feature in features:
                try:
                    yield feature[attribute]
                except KeyError:
                    pass
        except TypeError:
            # The features aren't iterable
            pass

    @staticmethod
    def _get_feature_attribute_with_value(features, attribute, value):
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

    @staticmethod
    def _get_feature_attribute_startswith_value(features, attribute, value):
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

    @staticmethod
    def _consolidated_series_to_dataframe(series):
        dataframe = series.reset_index()
        dataframe['criteria_full'] = dataframe[0].map(lambda x: x[0])
        dataframe['event_id'] = dataframe[0].map(lambda x: x[1])
        dataframe['criteria'] = dataframe[0].map(lambda x: x[0].split(',')[0])
        dataframe['criteria_additional'] = dataframe['criteria_full'].map(
            lambda x: x.split(',')[1] if len(x.split(',')) > 1 else np.nan)
        dataframe = dataframe.drop(0, axis=1)
        return dataframe

    def _consolidate_junction_events(self, df, db, event_col='event_id'):
        """Given events that share the same junctions, pick the best one"""
        if len(df) == 1:
            return 'only one', df[event_col].values[0]

        df_isoforms = df[self.transcript_cols].applymap(
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
                itertools.chain(*self._get_attribute(x, 'tag')))
            if not isinstance(x, float) else x)

        df_tags = df_tags.applymap(
            lambda x: x if not isinstance(x, list) or len(x) > 0 else np.nan)
        df_tags = df_tags.dropna(how='all')
        if df_tags.empty:
            return 'random,df_isoforms', df_isoforms.loc[
                np.random.choice(df_isoforms.index)]

        for tag in self.best_tags:
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

    def consolidate(self):
        """Of events that share the same junctions, pick the best exons"""
        consolidated = self.events.groupby(self.junctions).apply(
            lambda x: self._consolidate_junction_events(
                x, self.db, event_col='event_id'))
        consolidated_df = self._consolidated_series_to_dataframe(
            consolidated)
        return consolidated_df

