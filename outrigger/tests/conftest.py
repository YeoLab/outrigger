
import glob
import os

import gffutils
import pandas as pd
import pytest


@pytest.fixture
def outrigger_folder():
    """Absolute path root folder of outrigger package"""
    return os.path.join(os.path.abspath(os.path.dirname(__file__)), '../../')


@pytest.fixture
def data_folder():
    """Absolute path to where test data is stored"""
    return os.path.join(os.path.abspath(os.path.dirname(__file__)),
                        './data')


@pytest.fixture
def tasic2016(data_folder):
    """Suffix for data created from Tasic et al, Nat Neurosci (2016)"""
    return os.path.join(data_folder, 'tasic2016')


@pytest.fixture
def tasic2016_unprocessed(tasic2016):
    """Suffix for raw data created from Tasic et al, Nat Neurosci (2016)"""
    return os.path.join(tasic2016, 'unprocessed')


@pytest.fixture
def simulated(data_folder):
    """Suffix for simulated data"""
    return os.path.join(data_folder, 'simulated')


@pytest.fixture
def simulated_unprocessed(simulated):
    """Suffix for raw simulated data"""
    return os.path.join(simulated, 'unprocessed')


@pytest.fixture
def simulated_outrigger_output(simulated):
    """Suffix for simulated data outrigger output"""
    return os.path.join(simulated, 'outrigger_output')


@pytest.fixture
def simulated_outrigger_index(simulated_outrigger_output):
    """Suffix for simulated data splicing events"""
    return os.path.join(simulated_outrigger_output, 'index')


@pytest.fixture
def simulated_outrigger_se(simulated_outrigger_index):
    """Suffix for simulated data splicing events"""
    return os.path.join(simulated_outrigger_index, 'se')


@pytest.fixture
def simulated_outrigger_mxe(simulated_outrigger_index):
    """Suffix for simulated data splicing events"""
    return os.path.join(simulated_outrigger_index, 'mxe')


@pytest.fixture
def simulated_outrigger_psi(simulated_outrigger_output):
    """Suffix for simulated data splicing events"""
    return os.path.join(simulated_outrigger_output, 'psi')


@pytest.fixture
def negative_control_folder(data_folder):
    return os.path.join(data_folder, 'simulated', 'validate_negative_control')


@pytest.fixture
def tasic2016_sj_out_tab(tasic2016_unprocessed):
    return os.path.join(tasic2016_unprocessed, 'sj_out_tab')


@pytest.fixture
def tasic2016_gtf(tasic2016_unprocessed):
    return os.path.join(tasic2016_unprocessed, 'gtf')


@pytest.fixture
def tasic2016_intermediate(tasic2016):
    """Suffix for intermediate files from Tasic et al Nat Neurosci (2016)"""
    return os.path.join(tasic2016, 'intermediate')


@pytest.fixture
def tasic2016_outrigger_output(tasic2016):
    """Suffix for outrigger_output files from Tasic Nat Neurosci (2016)"""
    return os.path.join(tasic2016, 'outrigger_output')


@pytest.fixture
def tasic2016_outrigger_junctions(tasic2016_outrigger_output):
    return os.path.join(tasic2016_outrigger_output, 'junctions')


@pytest.fixture(params=['se', 'mxe'])
def splice_type(request):
    return request.param


@pytest.fixture(params=['positive', 'negative'])
def strand(request):
    if request.param == 'positive':
        return '+'
    else:
        return '-'


@pytest.fixture
def sj_filenames(tasic2016_sj_out_tab):
    globber = os.path.join(tasic2016_sj_out_tab, '*SJ.out.tab')
    return glob.glob(globber)


@pytest.fixture
def junction_reads(tasic2016_outrigger_junctions):
    filename = os.path.join(tasic2016_outrigger_junctions, 'reads.csv')
    return pd.read_csv(filename)


@pytest.fixture
def junction_metadata(tasic2016_outrigger_junctions):
    filename = os.path.join(tasic2016_outrigger_junctions, 'metadata.csv')
    return pd.read_csv(filename)


# -- IO: GTF input folders -- #
@pytest.fixture
def gtf_filename(tasic2016_gtf):
    return os.path.join(tasic2016_gtf,
                        'gencode.vM10.annotation.subset.gtf')


@pytest.fixture
def db_filename(gtf_filename):
    return gtf_filename + '.db'


@pytest.fixture
def db(db_filename):
    return gffutils.FeatureDB(db_filename)


@pytest.fixture
def snap25_exon_id():
    """Exon1 of MXE event in SNAP25. Should be in the gffutils db"""
    return 'exon:chr2:136763573-136763621:+'


@pytest.fixture
def myl6_novel_exon_junction():
    """Junction downstream of a novel exon in MYL6 test set"""
    return 'junction:chr10:128491033-128491719:-'


def pytest_addoption(parser):
    parser.addoption("--skip-slow", action="store_true",
                     help="Don't run the 'slow' (~15m) tests")

