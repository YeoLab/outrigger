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
    from outrigger.splicestrength import splice_site_sequences

    seqs = splice_site_sequences(exons, splice_site, genome_fasta, g=genome)

    string = six.StringIO()
    SeqIO.write(seqs, string, 'fasta')
    test = string.getvalue()

    if splice_site == 5:
        true = """>exon4_positive
GAGAATGGC
>exon4_negative
GGTTCAATT
>exon1alt_positive
TGGCTTGTA
>exon1alt_negative
ATCTATGGC
>exon2a5ss_positive
ACTGTCACT
>exon2a5ss_negative
TGAATGAGA
>exon1_positive
AGAGAACCC
>exon1_negative
TCAGCAAGC
>exon3_positive
TGTAACAAC
>exon3_negative
CATAGGCGA
>exon2_positive
GTTGTCGAT
>exon2_negative
TGAATGAGA
>exon2a3ss_positive
GTTGTCGAT
>exon2a3ss_negative
AATTTATTC
>exon4alt_positive
TAATCTCTA
>exon4alt_negative
TCCCACCAA
"""
    elif splice_site == 3:
        true = """>exon4_positive
AGGCGTACCATAGTAATTGAACC
>exon4_negative
ACAAACTGTAGGGTGCCATTCTC
>exon1alt_positive
ATAATTACACAATAGCCATAGAT
>exon1alt_negative
GCCGTTTGTGCCTATACAAGCCA
>exon2a5ss_positive
TTGGGGAAGTAGTATCTCATTCA
>exon2a5ss_negative
GAAAGAAGCAGCCTAGTGACAGT
>exon1_positive
TATAGGCACAAACGGCTTGCTGA
>exon1_negative
TCCCTCCACGTTATGGGTTCTCT
>exon3_positive
CTAGGCTGCTTCTTTCGCCTATG
>exon3_negative
TATGGAGCGATGACGTTGTTACA
>exon2_positive
TTGGGGAAGTAGTATCTCATTCA
>exon2_negative
TGAGAGCTTTGACTATCGACAAC
>exon2a3ss_positive
CCATAACGTGGAGGGAATAAATT
>exon2a3ss_negative
TGAGAGCTTTGACTATCGACAAC
>exon4alt_positive
AGGCCGGGCAGTGATTGGTGGGA
>exon4alt_negative
TGGGTCAATAAATTTAGAGATTA
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
