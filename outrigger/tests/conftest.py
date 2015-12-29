import glob
import os

import pytest
import sj2psi


@pytest.fixture
def data_folder():
    """Absolute path to where test data is stored"""
    return os.path.join(os.path.abspath(os.path.dirname(__file__)), 'data')


@pytest.fixture
def positive_strand_gtf_filename(data_folder):
    """Gene annotation file for positive strand tests"""
    return '{}/gencode.v19.rps24.positive.strand.gtf'


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
def gtf_filename(strand, positive_strand_gtf_filename,
        negative_strand_gtf_filename):
    if strand == '+':
        return positive_strand_gtf_filename
    elif strand == '-':
        return negative_strand_gtf_filename


@pytest.fixture
def db(gtf_filename):
    from outrigger.gtf import create_db
    return create_db(gtf_filename)


@pytest.fixture
def splice_junctions(sj_filenames):
    return sj2psi.read_multiple_sj_out_tab(sj_filenames)


@pytest.fixture
def metadata(splice_junctions):
    return sj2psi.sj_count_to_metadata(splice_junctions)
