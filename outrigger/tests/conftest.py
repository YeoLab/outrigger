import os

import gffutils
import pytest


@pytest.fixture
def data_folder():
    """Absolute path to where test data is stored"""
    return os.path.join(os.path.abspath(os.path.dirname(__file__)), 'data')


@pytest.fixture
def treutlein_folder(data_folder):
    return '{}/treutlein2014'.format(data_folder)


@pytest.fixture
def positive_strand_gtf_filename(data_folder):
    """Gene annotation file for positive strand tests"""
    return '{}/gencode.v19.rps24.positive.strand.gtf'.format(data_folder)


@pytest.fixture
def negative_strand_gtf_filename(data_folder):
    """Gene annotation file for negative strand tests"""
    return '{}/gencode.v19.pkm.negative.strand.gtf'


@pytest.fixture(params=['positive', 'negative'])
def strand(request):
    if request.param == 'positive':
        return '+'
    else:
        return '-'


@pytest.fixture
def gtf_filename(treutlein_folder):
    return '{}/gencode.vM2.annotation.fgfr2.gtf'.format(treutlein_folder)


@pytest.fixture
def db_filename(treutlein_folder):
    return '{}/gencode.vM2.annotation.fgfr2.gtf.db'.format(treutlein_folder)


@pytest.fixture
def db(db_filename):

    return gffutils.FeatureDB(db_filename)


@pytest.fixture
def splice_junctions(sj_filenames):
    from outrigger.io import star

    return star.read_multiple_sj_out_tab(sj_filenames)


@pytest.fixture
def metadata(splice_junctions):
    from outrigger.io import star

    return star.sj_count_to_metadata(splice_junctions)
