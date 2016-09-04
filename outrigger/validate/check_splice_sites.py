
from Bio import SeqIO
import pandas as pd

NT = 2


def read_splice_sites(bed, genome, fasta, direction='upstream'):
    """Read splice sites of an exon

    Parameters
    ----------
    bed : pybedtools.BedTool
        Exons whose splice sites you're interested in
    genome : str
        Location of a chromosome sizes file for the genome
    fasta : str
        Location of the genome fasta file
    direction : 'upstream' | 'downstream'
        Which direction of splice sites you want for the exon
    """

    if direction == 'upstream':
        left = NT
        right = 0
    elif direction == 'downstream':
        left = 0
        right =NT

    flanked = bed.flank(l=left, r=right, s=True, genome=genome)
    seqs = flanked.sequence(fi=fasta, s=True)

    with open(seqs.seqfn) as f:
        records = SeqIO.parse(f, 'fasta')
        records = pd.Series([str(x.seq) for x in records],
                            index=[x.name for x in bed])
    return records
