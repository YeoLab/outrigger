import glob
import os

import gffutils
import pandas as pd
import pytest


@pytest.fixture
def data_folder():
    """Absolute path to where test data is stored"""
    return os.path.join(os.path.abspath(os.path.dirname(__file__)),
                        '../test_data')


@pytest.fixture
def tasic2016(data_folder):
    """Suffix for data created from Tasic et al, Nat Neurosci (2016)"""
    return os.path.join(data_folder, 'tasic2016')


@pytest.fixture
def tasic2016_unprocessed(tasic2016):
    """Suffix for raw data created from Tasic et al, Nat Neurosci (2016)"""
    return os.path.join(tasic2016, 'unprocessed')

@pytest.fixture
def tasic2016_intermediate(tasic2016):
    """Suffix for intermediate files from Tasic et al Nat Neurosci (2016)"""
    return os.path.join(tasic2016, 'intermediate')


@pytest.fixture
def tasic2016_outrigger_output(tasic2016):
    """Suffix for outrigger_output files from Tasic et al Nat Neurosci (2016)"""
    return os.path.join(tasic2016, 'outrigger_output')


@pytest.fixture(params=['positive', 'negative'])
def strand(request):
    if request.param == 'positive':
        return '+'
    else:
        return '-'


@pytest.fixture
def sj_filenames(tasic2016_unprocessed):
    globber = os.path.join(tasic2016_unprocessed, 'sj_out_tab', '*SJ.out.tab')
    return glob.glob(globber)


@pytest.fixture
def splice_junctions(tasic2016_outrigger_output):
    filename = os.path.join(tasic2016_outrigger_output, 'junction_reads.csv')
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


@pytest.fixture
def snap25_exon():
    """Exon1 of MXE event in SNAP25. Should be in the gffutils db"""
    return 'exon:chr2:136763573-136763621:+'
