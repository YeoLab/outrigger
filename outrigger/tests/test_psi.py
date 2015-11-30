import numpy as np
import pytest
import pandas as pd
import pandas.util.testing as pdt
import six

@pytest.fixture
def exons_to_junctions():
    exons_to_junctions = {
        ('exon:chr1:150-175:+',  # Exon 1
         'exon:chr1:225-250:+',  # Exon 2
         'exon:chr1:300-350:+'):  # Exon 3
            ('junction:chr1:176-224:+',
             'junction:chr1:251-299:+',
             'junction:chr1:176-299:+'),
        ('exon:chr1:150-175:+',  # Exon 1
         'exon:chr1:225-275:+',  # Exon 2, alt 5' splice site
         'exon:chr1:300-350:+'):  # Exon 3
            ('junction:chr1:176-224:+',
             'junction:chr1:276-299:+',
             'junction:chr1:176-299:+'),
        ('exon:chr1:150-175:+',  # Exon 1
         'exon:chr1:300-350:+',  # Exon 3
         'exon:chr1:400-425:+'):  # Exon 4
            ('junction:chr1:176-299:+',
             'junction:chr1:351-399:+',
             'junction:chr1:176-399:+'),
        ('exon:chr1:150-175:+',  # Exon 1
         'exon:chr1:225-250:+',  # Exon 2
         'exon:chr1:400-425:+'):  # Exon 4
            ('junction:chr1:176-224:+',
             'junction:chr1:251-399:+',
             'junction:chr1:176-399:+'),
        ('exon:chr1:225-250:+',  # Exon 2
         'exon:chr1:300-350:+',  # Exon 3
         'exon:chr1:400-425:+'):  # Exon 4
            ('junction:chr1:251-299:+',
             'junction:chr1:351-399:+',
             'junction:chr1:251-399:+'),
        ('exon:chr1:150-175:+',  # Exon 1
         'exon:chr1:200-250:+',  # Exon 2, alt 3' splice site
         'exon:chr1:300-350:+'):  # Exon 4
            ('junction:chr1:176-199:+',
             'junction:chr1:251-299:+',
             'junction:chr1:176-299:+')}


    n_exons_se = 3
    print len(exons_to_junctions)
    exons_to_junctions = pd.DataFrame(exons_to_junctions).T.reset_index()
    exons_to_junctions = exons_to_junctions.rename(
        columns=dict(('level_{}'.format(i), 'exon{}'.format(i + 1)) for i in
                     range(n_exons_se)))
    exons_to_junctions = exons_to_junctions.rename(
        columns={0: 'junction12', 1: 'junction23', 2: 'junction13'})
    exons_to_junctions['event_id'] = exons_to_junctions.exon1 + '@' \
                                       + exons_to_junctions.exon2 + '@' \
                                       + exons_to_junctions.exon3
    return exons_to_junctions


@pytest.fixture
def junction12():
    """Junction between exons 1 and 2"""
    return "junction:chr1:176-224:+"


@pytest.fixture
def junction23():
    """Junction between exons 2 and 3"""
    return 'junction:chr1:251-299:+'


@pytest.fixture
def junction13():
    """Junction between exons 1 and 3"""
    return 'junction:chr1:176-299:+'


@pytest.fixture
def event_name():
    return 'exon:chr1:150-175:+@exon:chr1:225-250:+@exon:chr1:300-350:+'


def test_psi1(junction12, junction23, junction13, exons_to_junctions):
    from outrigger.psi import calculate_psi

    s = """sample_id,junction,reads
sample1,{},100
sample1,{},200
""".format(junction12, junction23)
    splice_junction_reads = pd.read_csv(six.StringIO(s), comment='#')
    splice_junction_reads = splice_junction_reads.set_index(
        ['junction', 'sample_id'])
    splice_junction_reads = splice_junction_reads.sort_index()

    test = calculate_psi(exons_to_junctions, splice_junction_reads)
    true = pd.read_csv(six.StringIO("""
sample_id,exon:chr1:150-175:+@exon:chr1:225-250:+@exon:chr1:300-350:+
sample1,1.0"""), index_col=0)

    pdt.assert_frame_equal(test, true)

def test_psi1_junction13_not_enough_reads(junction12, junction23, junction13,
                                          exons_to_junctions):
    from outrigger.psi import calculate_psi

    s = """sample_id,junction,reads
sample1,{},100
sample1,{},200
sample1,{},2
""".format(junction12, junction23, junction13)
    splice_junction_reads = pd.read_csv(six.StringIO(s), comment='#')
    splice_junction_reads = splice_junction_reads.set_index(
        ['junction', 'sample_id'])
    splice_junction_reads = splice_junction_reads.sort_index()

    test = calculate_psi(exons_to_junctions, splice_junction_reads)
    true = pd.read_csv(six.StringIO("""
sample_id,exon:chr1:150-175:+@exon:chr1:225-250:+@exon:chr1:300-350:+
sample1,1.0"""), index_col=0)

    pdt.assert_frame_equal(test, true)

@pytest.fixture(params=[100, 2, np.nan],
                ids=['enough reads', 'not enough reads', 'not there'])
def junction12_reads(request):
    return request.param

@pytest.fixture(params=[100, 2, np.nan],
                ids=['enough reads', 'not enough reads', 'not there'])
def junction23_reads(request):
    return request.param

@pytest.fixture(params=[100, 2, np.nan],
                ids=['enough reads', 'not enough reads', 'not there'])
def junction13_reads(request):
    return request.param


def test_psi_se(junction12, junction12_reads,
                junction23, junction23_reads,
                junction13, junction13_reads, exons_to_junctions):
    from outrigger.psi import calculate_psi, MIN_READS

    s = """sample_id,junction,reads
sample1,{0},{1}
sample1,{2},{3}
sample1,{4},{5}""".format(junction12, junction12_reads,
                          junction23, junction23_reads,
                          junction13, junction13_reads)
    splice_junction_reads = pd.read_csv(six.StringIO(s), comment='#')
    splice_junction_reads = splice_junction_reads.dropna()
    splice_junction_reads = splice_junction_reads.set_index(
        ['junction', 'sample_id'])
    splice_junction_reads = splice_junction_reads.sort_index()

    reads12 = junction12_reads if junction12_reads >= MIN_READS else 0
    reads23 = junction23_reads if junction23_reads >= MIN_READS else 0
    reads13 = junction13_reads if junction13_reads >= MIN_READS else 0

    if reads12 == 0 or reads23 == 0:
        isoform2_reads = 0
    else:
        isoform2_reads = reads12 + reads23
    isoform1_reads = reads13

    # This tests whether both are greater than zero
    if isoform1_reads or isoform2_reads:
        true_psi = isoform2_reads/(isoform2_reads + 2.*isoform1_reads)
    else:
        true_psi = np.nan

    other_isoform1_psi = 0. if isoform1_reads > 0 else np.nan

    test = calculate_psi(exons_to_junctions, splice_junction_reads,
                         isoform1_junctions=['junction13'],
                         isoform2_junctions=['junction12', 'junction23'])
    true = pd.read_csv(six.StringIO("""
sample_id,exon:chr1:150-175:+@exon:chr1:200-250:+@exon:chr1:300-350:+,exon:chr1:150-175:+@exon:chr1:225-250:+@exon:chr1:300-350:+,exon:chr1:150-175:+@exon:chr1:225-275:+@exon:chr1:300-350:+
sample1,{1},{0},{1}""".format(true_psi, other_isoform1_psi)), index_col=0)
    true = true.dropna(axis=1)

    pdt.assert_frame_equal(test, true)

