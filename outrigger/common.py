
# Column names and constants used across outrigger

CHROM = 'chrom'
STRAND = 'strand'

EXON_START = 'exon_start'
EXON_STOP = 'exon_stop'

JUNCTION_ID = 'junction_id'
JUNCTION_START = 'junction_start'
JUNCTION_STOP = 'junction_stop'
JUNCTION_MOTIF = 'junction_motif'

READS = 'reads'

ANNOTATED = 'annotated'

SAMPLE_ID = 'sample_id'

SPLICE_TYPES = (('skipped_exon', 'se'), ('mutually_exclusive_exon', 'mxe'))
ISOFORM_ORDER = 'isoform1', 'isoform2'
ISOFORM_COMPONENTS = {
    'se': {'isoform1': ('junction13',),
           'isoform2': ('junction12', 'exon2', 'junction23')},
    'mxe': {'isoform1': ('junction13', 'exon3', 'junction34'),
            'isoform2': ('junction12', 'exon2', 'junction24')}}
EVENT_ID_COLUMN = 'event_id'
ILLEGAL_JUNCTIONS = 'illegal_junctions'
SPLICE_TYPE_ISOFORM_EXONS = {'se': {'isoform1': ['exon1', 'exon3'],
                                    'isoform2': ['exon1', 'exon2', 'exon3']},
                             'mxe': {'isoform1': ['exon1', 'exon3', 'exon4'],
                                     'isoform2': ['exon1', 'exon2', 'exon4']}}
MIN_READS = 10
SE_ISOFORM1_JUNCTIONS = ['junction13']
SE_ISOFORM2_JUNCTIONS = ['junction12', 'junction23']
MXE_ISOFORM1_JUNCTIONS = ['junction13', 'junction34']
MXE_ISOFORM2_JUNCTIONS = ['junction12', 'junction24']
ISOFORM_JUNCTIONS = {
    'se': {'isoform1_junctions': SE_ISOFORM1_JUNCTIONS,
           'isoform2_junctions': SE_ISOFORM2_JUNCTIONS},
    'mxe': {'isoform1_junctions': MXE_ISOFORM1_JUNCTIONS,
            'isoform2_junctions': MXE_ISOFORM2_JUNCTIONS}}
