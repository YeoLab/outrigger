
import glob
import os

import itertools
import pandas as pd

from .region import Region

SPLICE_TYPE_ISOFORM_EXONS = {'SE': {'isoform1': ('exon1', 'exon3'),
                                    'isoform2': ('exon1', 'exon2', 'exon3')},
                             'MXE': {'isoform1': ('exon1', 'exon3', 'exon4'),
                                     'isoform2': ('exon1', 'exon2', 'exon4')}
                             }


class SplicingAnnotator(object):
    """Annotates basic features of splicing events: gene ids and names"""

    def __init__(self, db, events, splice_type):
        self.db = db
        self.events = events
        self.splice_type = splice_type
        self.isoform_exons = SPLICE_TYPE_ISOFORM_EXONS[self.splice_type]

    # def get_gene_attributes(self):
    #     """Iterate over all events in all types and get gene ids"""
    #     for splice_type, df in self.event_dfs.items():
    #         exon_cols = SPLICE_TYPE_ISOFORM_EXONS[splice_type]
    #         sys.stdout.write('{timestamp}\tGetting gene attributes for '
    #                          '{splice_type} events ...\n'.format(
    #                     timestamp=timestamp(), splice_type=splice_type))
    #         annotated = self._get_gene_attributes(df, exon_cols)
    #         done()

    def _get_gene_attributes(self, df, exon_cols):
        """Get gene, transcript ids and names for one event dataframe"""
        df['gencode_id'] = df[exon_cols].apply(
            lambda x: ','.join(set(itertools.chain(
                *[self.db[i].attributes['gene_id'] for i in x]))), axis=1)
        df['transcript_id'] = df[exon_cols].apply(
            lambda x: ','.join(set(itertools.chain(
                *[self.db[i].attributes['transcript_id'] for i in x]))),
            axis=1)
        df['gene_name'] = df[exon_cols].apply(
            lambda x: ','.join(set(itertools.chain(
                *[self.db[i].attributes['gene_name'] for i in x]))), axis=1)
        df['ensembl_id'] = df[exon_cols].apply(
            lambda x: ','.join(set(itertools.chain(
                *[map(lambda y: y.split('.')[0],
                      self.db[i].attributes['gene_id']) for i in x]))),
            axis=1)

        for exon_col in exon_cols:
            df['{}_region'.format(exon_col)] = df[exon_col].map(Region)
            df['{}_length'.format(exon_col)] = df['{}_region'.format(
                exon_col)].map(len)
        df['strand'] = df[exon_cols[0]].str[-1]
        return df

    def _single_row_transcripts(self, row, isoform1_exons, isoform2_exons):

        isoform_to_exons = {
            'isoform1': map(lambda x: self.db[row[x]], isoform1_exons),
            'isoform2': map(lambda x: self.db[row[x]], isoform2_exons)}

        isoform_to_transcripts = dict((k, set.intersection(
            *map(lambda x: set(self.db.parents(x, featuretype='transcript')),
                 v)))
                                      for k, v in isoform_to_exons.items())

        return isoform_to_transcripts

    def _get_attributes(self):

        lines = []

        for row in self.events.iterrows():
            for isoform, exons in self.isoform_exons:
                exon1 = self.db[row[exons[0]]]
                other_exons = row[exons[1:]]
                attributes = {}
                for (k1, v1) in exon1.attributes.items():
                    v2 = set(itertools.chain(*[self.db[i].attributes[k1]
                                               for i in other_exons]))
                    value = set(v1) & v2
                    if len(value) > 0:
                        attributes[k1] = ','.join(sorted(list(value)))
                attributes = pd.Series(attributes, name=row['exons'])
                attributes.index = isoform + '_' + attributes.index
                lines.append(exons)
        return pd.concat(attributes)
