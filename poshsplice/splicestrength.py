import os
import subprocess

from Bio import SeqIO
import pandas as pd
import pybedtools
import six

__author__ = 'olgabotvinnik'


VALID_SPLICE_SITES = (3, 5)


def get_splice_site(exons, genome, splice_site):
    """Get the 5' or 3' splice site of a set of exons

    The 5' splice site for MaxEntScan is defined as 9 bases:
        [3 bases in exon][6 bases in intron]
    The 3' splice site for MaxEntScan is defined as 23 bases:
        [20 bases in the intron][3 base in the exon]

    Parameters
    ----------
    exons : pybedtools.BedTool
        Exons for which you want to find the splice site locations
    genome : str
        Name of the genome to use, e.g. "hg19"
    splice_site : 5 | 3
        Either the 5' or 3' splice site

    Returns
    -------
    five_prime : pybedtools.BedTool
        5' Splice site of the exons
    """
    if splice_site not in VALID_SPLICE_SITES:
        raise ValueError('{0} is not a valid splice site. Only 5 and 3 are '
                         'acceptable'.format(splice_site))
    if splice_site == 5:
        left = 3
        right = 6
    elif splice_site == 3:
        left = 20
        right = 3

    return exons.flank(genome=genome, l=0, r=right,
                       s=True).slop(genome=genome, r=0, l=left, s=True)


def get_ss_sequence(exons, genome, splice_site, genome_fasta, filename=None):
    """Get the sequence of the 5' or 3' splice site

    Parameters
    ----------
    exons : pybedtools.BedTool
        Exons for which you want to find the splice site sequences
    genome : str
        Name of the genome to use, e.g. "hg19"
    splice_site : 5 | 3
        Either the 5' or 3' splice site
    genome_fasta : str
        Location of the (indexed) genome fasta file
    filename : str, optional
        Where to write the splice site sequences to, in fasta format

    Returns
    -------
    filename : str
        Location of the fasta file of the splice site sequences
    """
    splice_sites = get_splice_site(exons, genome, splice_site)
    splice_sites = splice_sites.sequence(fi=genome_fasta, s=True)
    if filename is None:
        return splice_sites.seqfn
    with open(filename, 'w') as f:
        f.write(open(splice_sites.seqfn).read())
    return filename


def score_splice_fasta(ss_fasta, splice_site, filename=None):
    """Get the Maximum Entropy scores of a splice site

    Parameters
    ----------
    ss_fasta : str
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

    ss_fasta = os.path.abspath(os.path.expanduser(ss_fasta))
    maxentscan_dir = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                                  os.path.join(*['external', 'maxentscan']))
    currdir = os.getcwd()
    os.chdir(maxentscan_dir)

    # Check that the input fasta files are the right length
    length = 9 if splice_site == 5 else 23
    parsed = SeqIO.parse(open(ss_fasta), 'fasta')
    if sum(len(x.seq) != length for x in parsed) > 0:
        raise ValueError("Not all the sequences in the fasta file are {0}nt "
                         "long. For {1}' splice sites, the required"
                         "input is fasta files with sequence "
                         "length {0}nt.".format(length, splice_site))

    program = 'score{0}.pl'.format(splice_site)
    pipe = subprocess.Popen(["perl", program, ss_fasta],
                            stdout=subprocess.PIPE)
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
    scores = []
    bed = pybedtools.BedTool(exons)
    df = pd.DataFrame(index=[x.name for x in bed])
    for splice_site in VALID_SPLICE_SITES:
        ss_seqs = get_ss_sequence(exons, genome, splice_site, genome_fasta)
        score = score_splice_fasta(ss_seqs, splice_site)
        score = read_splice_scores(scores)
        df['splice_site_score_{}p'.format(splice_site)] = score.values
    return df