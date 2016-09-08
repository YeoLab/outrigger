import collections
import os

import pandas as pd
import pandas.util.testing as pdt
import pytest

@pytest.fixture
def negative_control_folder(data_folder):
    return os.path.join(data_folder, 'simulated', 'validate_negative_control')


@pytest.fixture
def exon2_bed(negative_control_folder):
    return os.path.join(negative_control_folder, 'outrigger_output', 'index',
                        'se', 'exon2.bed')


@pytest.fixture
def do_not_have_mysql():
    """Checks if mysql exists on this platform"""
    for path in os.environ["PATH"].split(os.pathsep):
        if os.path.exists(os.path.join(path, 'msql')):
                return True
    return False


@pytest.fixture
def genome_name():
    return 'mm10'


@pytest.fixture(params=['upstream', 'downstream'])
def direction(request):
    return request.param

@pytest.fixture
def simulated_chromsizes(negative_control_folder):
    return os.path.join(negative_control_folder, 'chromsizes')


@pytest.fixture
def simulated_fasta(negative_control_folder):
    return os.path.join(negative_control_folder, 'genome.fasta')


@pytest.fixture(params=[
    pytest.mark.skipif('do_not_have_mysql')('genome_name'),
    'filename'])
def maybe_chromsizes(request, simulated_chromsizes, genome_name):
    if request.param == 'genome_name':
        return genome_name
    elif request.param == 'filename':
        return simulated_chromsizes


def test_splice_site_str_to_tuple():
    from outrigger.validate.check_splice_sites import splice_site_str_to_tuple

    test = splice_site_str_to_tuple('GT/AG,AT/AC')
    true = 'GT/AG', 'AT/AC'

    assert test == true


def test_maybe_read_chromsizes(maybe_chromsizes, genome_name):
    from outrigger.validate.check_splice_sites import maybe_read_chromsizes
    test = maybe_read_chromsizes(maybe_chromsizes)

    if maybe_chromsizes != genome_name:
        true = collections.OrderedDict({'simulated': (0, 1000)})
    else:
        chromsizes = '''chr1	195471971
chr2	182113224
chrX	171031299
chr3	160039680
chr4	156508116
chr5	151834684
chr6	149736546
chr7	145441459
chr10	130694993
chr8	129401213
chr14	124902244
chr9	124595110
chr11	122082543
chr13	120421639
chr12	120129022
chr15	104043685
chr16	98207768
chr17	94987271
chrY	91744698
chr18	90702639
chr19	61431566
chr5_JH584299_random	953012
chrX_GL456233_random	336933
chrY_JH584301_random	259875
chr1_GL456211_random	241735
chr4_GL456350_random	227966
chr4_JH584293_random	207968
chr1_GL456221_random	206961
chr5_JH584297_random	205776
chr5_JH584296_random	199368
chr5_GL456354_random	195993
chr4_JH584294_random	191905
chr5_JH584298_random	184189
chrY_JH584300_random	182347
chr7_GL456219_random	175968
chr1_GL456210_random	169725
chrY_JH584303_random	158099
chrY_JH584302_random	155838
chr1_GL456212_random	153618
chrUn_JH584304	114452
chrUn_GL456379	72385
chr4_GL456216_random	66673
chrUn_GL456393	55711
chrUn_GL456366	47073
chrUn_GL456367	42057
chrUn_GL456239	40056
chr1_GL456213_random	39340
chrUn_GL456383	38659
chrUn_GL456385	35240
chrUn_GL456360	31704
chrUn_GL456378	31602
chrUn_GL456389	28772
chrUn_GL456372	28664
chrUn_GL456370	26764
chrUn_GL456381	25871
chrUn_GL456387	24685
chrUn_GL456390	24668
chrUn_GL456394	24323
chrUn_GL456392	23629
chrUn_GL456382	23158
chrUn_GL456359	22974
chrUn_GL456396	21240
chrUn_GL456368	20208
chrM	16299
chr4_JH584292_random	14945
chr4_JH584295_random	1976'''.split('\n')
        chromsizes = [x.split() for x in chromsizes]
        true = collections.OrderedDict(
            [(chrom, (0, int(length))) for chrom, length in chromsizes])
    pdt.assert_dict_equal(test, true)


def test_read_splice_sites(exon2_bed, direction, simulated_fasta,
                           simulated_chromsizes, negative_control_folder):
    from outrigger.validate.check_splice_sites import read_splice_sites

    test = read_splice_sites(exon2_bed, simulated_chromsizes, simulated_fasta,
                             direction)

    csv = os.path.join(negative_control_folder,
                       'exon2_{}_splice_sites.csv'.format(direction))
    true = pd.read_csv(csv, index_col=0, squeeze=True, header=None, names=None)
    true.name = None
    true.index.name = None

    pdt.assert_series_equal(test, true)
