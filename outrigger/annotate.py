
import itertools
from outrigger.region import Region

SPLICE_TYPE_ISOFORM_EXONS = {'SE': {'isoform1': ('exon1', 'exon3'),
                                    'isoform2': ('exon1', 'exon2', 'exon3')},
                             'MXE': {'isoform1': ('exon1', 'exon3', 'exon4'),
                                     'isoform2': ('exon1', 'exon2', 'exon4')}
                             }

class SplicingAnnotator(object):
    """Annotates basic features of splicing events: gene ids and names"""
    pass

    @staticmethod
    def get_gene_attributes(df, exon_cols, db):
        df['gencode_id'] = df[exon_cols].apply(
            lambda x: ','.join(set(itertools.chain(*[db[i].attributes['gene_id'] for i in x]))), axis=1)
        df['transcript_id'] = df[exon_cols].apply(
            lambda x: ','.join(set(itertools.chain(*[db[i].attributes['transcript_id'] for i in x]))), axis=1)
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


mxe_df = pd.read_csv('{}/events/mxe.csv'.format(outrigger_folder))
print
mxe_df.shape
mxe_df.head()
mxe_exon_cols = ['exon{}'.format(i) for i in range(1, 5)]
print(mxe_exon_cols)
mxe_df_annotated = get_gene_attributes(mxe_df, mxe_exon_cols, v19db)

se_exon_cols = ['exon{}'.format(i) for i in range(1, 4)]
print(se_exon_cols)
se_df_annotated = get_gene_attributes(se_df, se_exon_cols, v19db)
