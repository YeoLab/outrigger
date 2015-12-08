import os
import subprocess
import tempfile

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
import pybedtools
import six

__author__ = 'olgabotvinnik'


VALID_SPLICE_SITES = 3, 5


def splice_site_sequences(exons, splice_site, genome_fasta, genome=None, g=None):
    """Get the sequence of the 5' or 3' splice site of the exons

    The 5' splice site for MaxEntScan is defined as 9 bases:
        [3 bases in exon][6 bases in intron]
    The 3' splice site for MaxEntScan is defined as 23 bases:
        [20 bases in the intron][3 base in the exon]

    Parameters
    ----------
    exons : pybedtools.BedTool | str
        Exons for which you want to find the splice site sequences, either a
        string to the full path of the bed file, or a pre-created pybedtool
    splice_site : 5 | 3
        Either the 5' or 3' splice site
    genome : str | OrderedDict
        Name of the genome to use, e.g. "hg19", which looks up the chromosome
        sizes of that genome. If a non-standard genome, must instead provide
        an OrderedDict of chromosome name to sizes, e.g.
        OrderedDict([('chr1', (0, 1000))])
    genome_fasta : str
        Location of the (indexed) genome fasta file

    Returns
    -------
    filename : str
        Location of the fasta file of the splice site sequences
    """
    if not isinstance(exons, pybedtools.BedTool):
        exons = pybedtools.BedTool(exons)
    names = [x.name for x in exons]

    if splice_site == 5:
        left = 3
        right = 6
    elif splice_site == 3:
        left = 21
        right = 2
    else:
        raise ValueError('{0} is not a valid splice site. Only 5 and 3 are '
                         'acceptable'.format(splice_site))

    if g is None and genome is None:
        raise ValueError('Either `g` or `genome` must be specified. `g` is '
                         'a chromosome sizes file, while `genome` is a string '
                         'like "hg19" that looks up the chromosome sizes from'
                         ' a database')
    kwargs = dict(g=g) if g is not None else dict(genome=genome)

    negative = exons.filter(lambda x: x.strand == '-')
    positive = exons.filter(lambda x: x.strand == '+')

    seqs = []
    beds = {'-': negative, '+': positive}
    for strand, bed in beds.items():
        if strand == '-':
            right += 1
            left -= 1

        if splice_site == 5:
            flanked = bed.flank(l=0, r=right, s=True, **kwargs)
            slopped = flanked.slop(r=0, l=left, s=True, **kwargs)
        else:
            flanked = bed.flank(l=left, r=0, s=True, **kwargs)
            slopped = flanked.slop(l=0, r=right, s=True, **kwargs)
        seq = slopped.sequence(fi=genome_fasta, s=True)
        with open(seq.seqfn) as f:
            records = SeqIO.parse(f, 'fasta')
            records = pd.Series([str(x.seq) for x in records],
                                index=[x.name for x in slopped])
        seqs.append(records)
    seqs = pd.concat(seqs)

    # Reorder the sequences into the original order
    seqs = seqs[names]

    reordered = [SeqRecord(Seq(s), id=name, description='')
                 for name, s in seqs.iteritems()]

    return reordered


def score_splice_fasta(fasta, splice_site, filename=None):
    """Get the Maximum Entropy scores of a splice site

    Parameters
    ----------
    fasta : str
        Location of the fasta files you want to test
    splice_site : 5 | 3
        Either the 5' or 3' splice site
    filename : str, optional
        Where to output the splice site scores to

    Returns
    -------
    scores : str


    """
    if splice_site not in VALID_SPLICE_SITES:
        raise ValueError('{0} is not a valid splice site. Only 5 and 3 are '
                         'acceptable'.format(splice_site))

    fasta = os.path.abspath(os.path.expanduser(fasta))
    maxentscan_dir = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                                  os.path.join(*['external', 'maxentscan']))
    currdir = os.getcwd()
    os.chdir(maxentscan_dir)

    # Check that the input fasta files are the right length
    length = 9 if splice_site == 5 else 23
    try:
        parsed = SeqIO.parse(open(fasta), 'fasta')
    except IOError:
        parsed = fasta
        tempfasta = tempfile.NamedTemporaryFile()
        with open(tempfasta.name, 'w') as f:
            SeqIO.write(parsed, f, 'fasta')
        fasta = tempfasta.name

    if sum(len(x.seq) != length for x in parsed) > 0:
        raise ValueError("Not all the sequences in the fasta file are {0}nt "
                         "long. For {1}' splice sites, the required"
                         "input is fasta files with sequence "
                         "length {0}nt.".format(length, splice_site))

    program = 'score{0}.pl'.format(splice_site)
    pipe = subprocess.Popen(["perl", program, fasta], stdout=subprocess.PIPE)
    output = pipe.stdout.read()

    if filename is not None:
        with open(filename, 'w') as f:
            f.write(output)
    os.chdir(currdir)
    if six.PY3:
        output = output.decode('utf-8')
    return output


def read_splice_scores(scores):
    """Read splice site scores and return a pandas series

    Parameters
    ----------
    scores : str
        Either the output from :py:func:`score_splice_fasta` or a filename of
        scores from the original MaxEntScan ``score{5,3}.pl`` functions

    Returns
    -------
    scores_series : pandas.Series
        A pandas series of the splicing scores, in exactly the order they
        appeared
    """
    if not os.path.exists(scores):
        filename = six.StringIO(scores)
    else:
        filename = scores
    return pd.read_table(filename, squeeze=True, header=None, index_col=0,
                         sep='\s+')


def score_exons(exons, genome, genome_fasta):
    """One-stop shop to score 5' and 3' ends of exons

    Get the splice site strength of the 5' and 3' end of the exons.

    Parameters
    ----------
    exons : str
        Full path to the exons bed file
    genome : str
        Genome build, e.g. 'hg19'
    genome_fasta : str
        Full path to the genome fasta file

    Returns
    -------
    scores : pandas.DataFrame
        A (n_exons, 4) dataframe of the 3' and 5' splice site scores and the
        sequences from which they are derived. The index (row names) of the
        dataframe is taken from the "name" field of the exon bed file, and the
        scores are in the exact same order.
    """
    if not isinstance(exons, pybedtools.BedTool):
        bed = pybedtools.BedTool(exons)
    else:
        bed = exons

    df = pd.DataFrame(index=[x.name for x in bed])

    for splice_site in VALID_SPLICE_SITES:
        ss_seqs = splice_site_sequences(bed, genome_fasta, splice_site, genome)
        scores = score_splice_fasta(ss_seqs, splice_site)
        scores = read_splice_scores(scores)
        df['splice_site_{}p_score'.format(splice_site)] = scores.values
        df['splice_site_{}p_seq'.format(splice_site)] = scores.index
    return df
