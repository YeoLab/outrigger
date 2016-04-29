import glob
import os

import gffutils
import pandas as pd
import pytest


@pytest.fixture
def data_folder():
    """Absolute path to where test data is stored"""
    return os.path.join(os.path.abspath(os.path.dirname(__file__)), 'data')


@pytest.fixture
def treutlein_folder(data_folder):
    return os.path.join(data_folder, 'treutlein2014')


@pytest.fixture
def sj_filenames(treutlein_folder):
    globber = os.path.join(treutlein_folder, 'sj_out_tab', '*SJ.out.tab')
    return glob.glob(globber)


@pytest.fixture
def splice_junctions(treutlein_folder):
    filename = os.path.join(treutlein_folder,
                            'splice_junctions_multimappingFalse.csv')
    return pd.read_csv(filename)


@pytest.fixture
def metadata(treutlein_folder):
    filename = os.path.join(treutlein_folder, 'junction_metadata.csv')
    return pd.read_csv(filename)


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
