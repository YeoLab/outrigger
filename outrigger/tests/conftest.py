import os

import pytest

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

@pytest.fixture(params=['negative', 'positive'])
def strand(request):
    return request.param

@pytest.fixture
def gtf_filename(strand, positive_strand_gtf_filename,
        negative_strand_gtf_filename):
    if strand == 'positive':
        return positive_strand_gtf_filename
    elif strand == 'negative':
        return negative_strand_gtf_filename

@pytest.fixture
def db(gtf_filename):
    from outrigger.gtf import create_db
    return create_db(gtf_filename)
