import os

from Bio import SeqIO
from Bio.Alphabet import Alphabet

import pandas as pd
import pandas.util.testing as pdt
import pybedtools
import pytest
import six

__author__ = 'olgabotvinnik'


@pytest.fixture
def example_three_prime_fasta(tmpdir):
    s = '''> dummy1
ttccaaacgaacttttgtAGgga
> dummy2
tgtctttttctgtgtggcAGtgg
> dummy3
ttctctcttcagacttatAGcaa'''
    filename = '{0}/three_prime.fasta'.format(tmpdir)
    with open(filename, 'w') as f:
        f.write(s)
    return filename


@pytest.fixture
def example_five_prime_fasta(tmpdir):
    s = '''> dummy1
cagGTAAGT
> dummy2
gagGTAAGT
> dummy3
taaATAAGT'''
    filename = '{0}/five_prime.fasta'.format(tmpdir)
    with open(filename, 'w') as f:
        f.write(s)
    return filename


@pytest.fixture(params=[5, 3])
def splice_site_combo(request, example_five_prime_fasta,
                      example_three_prime_fasta):
    if request.param == 5:
        # From actually running the perl program
        true = six.u('''cagGTAAGT	10.86
gagGTAAGT	11.08
taaATAAGT	-0.12
''')
        return example_five_prime_fasta, request.param, true
    if request.param == 3:
        # From actually running the perl program
        true = six.u('''ttccaaacgaacttttgtAGgga	2.89
tgtctttttctgtgtggcAGtgg	8.19
ttctctcttcagacttatAGcaa	-0.08
''')
        return example_three_prime_fasta, request.param, true


@pytest.fixture
def dirname():
    return os.path.dirname(__file__)


@pytest.fixture
def bed_filename(dirname):
    return '{}/test.bed'.format(dirname)


@pytest.fixture(params=['BedTool', 'filename'])
def exons(request, bed_filename):
    if request.param == 'filename':
        return bed_filename
    # elif request == 'BedTool':
    else:
        return pybedtools.BedTool(bed_filename)


@pytest.fixture
def genome(dirname):
    return '{}/test.chromsizes'.format(dirname)


@pytest.fixture
def genome_fasta(dirname):
    return '{}/test.fasta'.format(dirname)


@pytest.fixture(params=[5, 3])
def splice_site(request):
    return request.param


def test_get_ss_sequence(exons, genome, splice_site, genome_fasta):
    from outrigger.splicestrength import get_ss_sequence

    seqs = get_ss_sequence(exons, splice_site, genome_fasta, g=genome)

    string = six.StringIO()
    SeqIO.write(seqs, string, 'fasta')
    test = string.getvalue()

    if splice_site == 5:
        true = """>exon4_positive
AAGGCCATC
>exon4_negative
TTCGACCAC
>exon1alt_positive
GACATCCAC
>exon1alt_negative
AGACCCGTT
>exon2a5ss_positive
GCATATCAA
>exon2a5ss_negative
CTCTACTCC
>exon1_positive
TATTGTGTC
>exon1_negative
CGTCGCATA
>exon3_positive
TCTAGTAAT
>exon3_negative
CGCCTTCTG
>exon2_positive
GCGTGTACC
>exon2_negative
CTCTACTCC
>exon2a3ss_positive
GCGTGTACC
>exon2a3ss_negative
TCAGTGTAG
>exon4alt_positive
AATCATTAC
>exon4alt_negative
TCCAGTCAA
"""
    elif splice_site == 3:
        true = """>exon4_positive
ACTATAATGCGTAAGTGGTCGAA
>exon4_negative
CATTGCCTTGCTCAGATGGCCTT
>exon1alt_positive
TCCCTACTGTCCCCAACGGGTCT
>exon1alt_negative
TAGCCGAGTCCTAAGTGGATGTC
>exon2a5ss_positive
AGTGCATGTCGTTCGGAGTAGAG
>exon2a5ss_negative
TGACGTTCGAACGATTGATATGC
>exon1_positive
ACTTAGGACTCGGCTATGCGACG
>exon1_negative
AGGGGGAACAGGTTGACACAATA
>exon3_positive
AATCGTTCGAACGTCAGAAGGCG
>exon3_negative
TAGAAACTATTATCATTACTAGA
>exon2_positive
AGTGCATGTCGTTCGGAGTAGAG
>exon2_negative
TGGTATGGACATGAGGTACACGC
>exon2a3ss_positive
TCAACCTGTTCCCCCTACACTGA
>exon2a3ss_negative
TGGTATGGACATGAGGTACACGC
>exon4alt_positive
AGGCTTGACGCATTTTGACTGGA
>exon4alt_negative
TGTCGAATACTACTGTAATGATT
"""

    assert test == true


def test_score_splice_fasta(splice_site_combo):
    from outrigger.splicestrength import score_splice_fasta

    fasta, splice_site, true = splice_site_combo
    test = score_splice_fasta(fasta, splice_site)
    pdt.assert_equal(test, true)


@pytest.fixture(params=['string', 'filename'])
def splice_scores(request, splice_site_combo, tmpdir):
    from outrigger.splicestrength import score_splice_fasta
    fasta, splice_site, true = splice_site_combo
    scores = score_splice_fasta(fasta, splice_site)
    if request.param == 'string':
        return scores
    elif request.param == 'filename':
        filename = '{0}/splice_scores.txt'.format(tmpdir)
        with open(filename, 'w') as f:
            f.write(scores)
        return filename


def test_read_splice_scores(splice_scores):
    from outrigger.splicestrength import read_splice_scores

    test = read_splice_scores(splice_scores)
    if not os.path.exists(splice_scores):
        filename = six.StringIO(splice_scores)
    else:
        filename = splice_scores
    true = pd.read_table(filename, squeeze=True, header=None,
                         index_col=0, sep='\s+')
    pdt.assert_series_equal(test, true)
