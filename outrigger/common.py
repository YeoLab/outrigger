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
EVENT_ID = 'event_id'
INCOMPATIBLE_JUNCTIONS = 'incompatible_junctions'
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


SPLICE_ABBREVS = 'se', 'mxe'
SPLICE_TYPE_ALL_EXONS = {'se': ['exon1', 'exon2', 'exon3'],
                         'mxe': ['exon1', 'exon2', 'exon3', 'exon4']}
SPLICE_TYPE_ALL_JUNCTIONS = {'se': ['junction13', 'junction12', 'junction23'],
                             'mxe': ['junction13',
                                     'junction34',
                                     'junction12',
                                     'junction24']}

# for gffutils to sort output features by something so that output files are
# comparable
ORDER_BY = ('seqid', 'start', 'end', 'frame', 'source', 'strand', 'attributes')
UNIQUE_READS = 'unique_junction_reads'
MULTIMAP_READS = 'multimap_junction_reads'
MAX_OVERHANG = 'max_overhang'

UNEVEN_COVERAGE_MULTIPLIER = 10


# --- Outrigger Psi --- #
NOTES = 'notes'
PSI = 'psi'
UPSTREAM = 'upstream'
DOWNSTREAM = 'downstream'
DIRECTIONS = UPSTREAM, DOWNSTREAM
NOVEL_EXON = 'novel_exon'
OUTRIGGER_DE_NOVO = 'outrigger_de_novo'
MAX_DE_NOVO_EXON_LENGTH = 100
