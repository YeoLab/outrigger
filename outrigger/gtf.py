import gffutils
import datetime

# Annotations from:
# ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz

gene_transcript = set(('gene', 'transcript'))


def transform(f):
    if f.featuretype in gene_transcript:
        return f
    else:
        exon_location = '{}:{}:{}-{}:{}'.format(f.featuretype, f.seqid, f.start, f.stop, f.strand)
        exon_id = exon_location
        if f.featuretype == 'CDS':
            exon_id += ':' + f.frame
        f.attributes['fancy_id'] = [exon_id]
        return f


def create_db(gtf_filename, db_filename=None):
    db_filename = ':memory:' if db_filename is None else db_filename

    return gffutils.create_db(
        gtf_filename,
        db_filename,
        merge_strategy='merge',
        id_spec={'gene': 'gene_id', 'transcript':
            'transcript_id', 'exon': 'fancy_id',
                 'CDS': 'fancy_id', 'start_codon': 'fancy_id',
                 'stop_codon': 'fancy_id', 'UTR': 'fancy_id'},
        transform=transform,
        force=True,
        verbose=True,
        infer_gene_extent=False,
        force_merge_fields=['source'])
