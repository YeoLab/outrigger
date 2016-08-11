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
def treutlein():
    """Suffix for data created from Treutlein et al, 2014 (Nature) paper"""
    return 'treutlein2014'


@pytest.fixture(params=['positive', 'negative'])
def strand(request):
    if request.param == 'positive':
        return '+'
    else:
        return '-'


# - Input/Output (IO) test data folders - #
@pytest.fixture
def io_folder(data_folder):
    return os.path.join(data_folder, 'io')


# -- IO: STAR SJ.out tab folders -- #
@pytest.fixture
def star_folder(io_folder):
    return os.path.join(io_folder, 'star')


@pytest.fixture
def treutlein_star(star_folder, treutlein):
    return os.path.join(star_folder, treutlein)


@pytest.fixture
def sj_filenames(treutlein_star):
    globber = os.path.join(treutlein_star, 'sj_out_tab', '*SJ.out.tab')
    return glob.glob(globber)


@pytest.fixture
def splice_junctions(treutlein_star):
    filename = os.path.join(treutlein_star,
                            'splice_junctions_multimappingFalse.csv')
    return pd.read_csv(filename)


@pytest.fixture
def metadata(treutlein_star):
    filename = os.path.join(treutlein_star, 'junction_metadata.csv')
    return pd.read_csv(filename)


# -- IO: GTF input folders -- #
@pytest.fixture
def gtf_folder(io_folder):
    return os.path.join(io_folder, 'gtf')


@pytest.fixture
def treutlein_gtf(io_folder, treutlein):
    return os.path.join(io_folder, 'gtf', treutlein)


@pytest.fixture
def gtf_filename(treutlein_gtf):
    return os.path.join(treutlein_gtf, 'gencode.vM2.annotation.fgfr2.gtf')


@pytest.fixture
def db_filename(treutlein_gtf):
    return os.path.join(treutlein_gtf, 'gencode.vM2.annotation.fgfr2.gtf.db')


@pytest.fixture
def db(db_filename):
    return gffutils.FeatureDB(db_filename)


# - Outrigger index test data folders - #
@pytest.fixture
def index_folder(data_folder):
    return os.path.join(data_folder, 'index')


@pytest.fixture
def adjacencies_folder(index_folder):
    return os.path.join(index_folder, 'adjacencies')


@pytest.fixture
def treutlein_adjacencies(adjacencies_folder, treutlein):
    return os.path.join(adjacencies_folder, treutlein)


# - Outrigger psi test data folders - #
@pytest.fixture
def psi_folder(data_folder):
    return os.path.join(data_folder, 'psi')
