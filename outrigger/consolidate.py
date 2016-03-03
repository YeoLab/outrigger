"""
Consolidate alternative splicing events which use the same alternative exon,
but different flanking exons, based on biological annotation of the
transcripts.
"""

import itertools

import numpy as np


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
