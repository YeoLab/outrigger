__author__ = 'olgabotvinnik'

import pandas.util.testing as pdt
import pytest

@pytest.fixture
def example_three_prime_fasta(tmpdir):
    s = '''> dummy1
ttccaaacgaacttttgtAGgga
> dummy2
tgtctttttctgtgtggcAGtgg
> dummy3
ttctctcttcagacttatAGcaa'''
    filename = '{}/three_prime.fasta'.format(tmpdir)
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
    filename = '{}/five_prime.fasta'.format(tmpdir)
    with open(filename, 'w') as f:
        f.write(s)
    return filename

@pytest.fixture(params=[5, 3])
def splice_site_combo(request, example_five_prime_fasta, example_three_prime_fasta):
    if request.param == 5:
        # From actually running the perl program
        true = '''cagGTAAGT	10.86
gagGTAAGT	11.08
taaATAAGT	-0.12
'''
        return example_five_prime_fasta, request.param, true
    if request.param == 3:
        # From actually running the perl program
        true = '''ttccaaacgaacttttgtAGgga	2.89
tgtctttttctgtgtggcAGtgg	8.19
ttctctcttcagacttatAGcaa	-0.08
'''
        return example_three_prime_fasta, request.param, true


def test_score_splice_fasta(splice_site_combo):
    from poshsplice.splicestrength import score_splice_fasta

    fasta, splice_site, true = splice_site_combo
    test = score_splice_fasta(fasta, splice_site)
    pdt.assert_equal(test, true)