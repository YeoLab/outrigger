
import glob
import os

import itertools
import pandas as pd

from outrigger.region import Region

SPLICE_TYPE_ISOFORM_EXONS = {'SE': {'isoform1': ('exon1', 'exon3'),
                                    'isoform2': ('exon1', 'exon2', 'exon3')},
                             'MXE': {'isoform1': ('exon1', 'exon3', 'exon4'),
                                     'isoform2': ('exon1', 'exon2', 'exon4')}
                             }

class SplicingAnnotator(object):
    """Annotates basic features of splicing events: gene ids and names"""
    pass

    def __init__(self, outrigger_folder, events_folder, db):

        for filename in glob.iglob('{}/*.csv'.format(events_folder)):
            event_type = os.path.basename(filename).split('.csv')[0]

        mxe_df = pd.read_csv('{}/mxe.csv'.format(events_folder))
        mxe_df.head()
        mxe_exon_cols = ['exon{}'.format(i) for i in range(1, 5)]
        print(mxe_exon_cols)
        mxe_df_annotated = self.get_gene_attributes(mxe_df, mxe_exon_cols, db)

        se_exon_cols = ['exon{}'.format(i) for i in range(1, 4)]
        print(se_exon_cols)
        se_df_annotated = self.get_gene_attributes(se_df, se_exon_cols, db)

    @staticmethod
    def get_gene_attributes(df, exon_cols, db):
        df['gencode_id'] = df[exon_cols].apply(
            lambda x: ','.join(set(itertools.chain(
                *[db[i].attributes['gene_id'] for i in x]))), axis=1)
        df['transcript_id'] = df[exon_cols].apply(
            lambda x: ','.join(set(itertools.chain(
                *[db[i].attributes['transcript_id'] for i in x]))), axis=1)
        df['gene_name'] = df[exon_cols].apply(
            lambda x: ','.join(set(itertools.chain(*[db[i].attributes['gene_name'] for i in x]))), axis=1)
        df['ensembl_id'] = df[exon_cols].apply(
            lambda x: ','.join(set(itertools.chain(*[map(lambda y: y.split('.')[0],
                                                         db[i].attributes['gene_id']) for i in x]))),
            axis=1)

        for exon_col in exon_cols:
            df['{}_region'.format(exon_col)] = df[exon_col].map(Region)
            df['{}_length'.format(exon_col)] = df['{}_region'.format(exon_col)].map(len)
        df['strand'] = df[exon_cols[0]].str[-1]
        return df

    @staticmethod
    def get_transcripts(df, isoform1_exons, isoform2_exons, db):

        isoform_to_exons = {
            'isoform1': map(lambda x: db[row[x]], isoform1_exons),
            'isoform2': map(lambda x: db[row[x]], isoform2_exons)}

        isoform_to_transcripts = dict((k, set.intersection(
            *map(lambda x: set(db.parents(x, featuretype='transcript')),
                 v)))
                                      for k, v in isoform_to_exons.items())

