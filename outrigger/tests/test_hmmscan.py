import os

import pandas as pd
import pandas.util.testing as pdt
import pytest

try:
    # For Python 2
    from StringIO import StringIO
except ImportError:
    # For Python 3
    from io import StringIO


@pytest.fixture
def example_hmmscan_out():
    directory = os.path.dirname(__file__)
    filename = '{0}/hmmscan_out.txt'.format(directory)
    return filename


def test_read_hmmscan_out(example_hmmscan_out):
    from outrigger.hmmscan import read_hmmscan

    test = read_hmmscan(example_hmmscan_out)
    true_s = """,target_name,target_accession,target_length,query_name,query_accession,query_length,sequence_e_value,sequence_score,sequence_bias,domain_number,domain_total,domain_conditional_e_value,domain_independent_e_value,domain_score,domain_bias,target_start,target_stop,query_start,query_stop,query_domain_envelope_start,query_domain_envelope_stop,mean_posterior_probability,target_description# noqa
0,Fox-1_C,PF12414.3,93,sp|O43251|RFOX2_HUMAN,-,390,3.2e-39,133.2,29.5,1,2,0.23,670.0,0.7,0.0,14,48,177,213,166,243,0.66,Calcitonin gene-related peptide regulator C terminal# noqa
1,Fox-1_C,PF12414.3,93,sp|O43251|RFOX2_HUMAN,-,390,3.2e-39,133.2,29.5,2,2,8.900000000000001e-42,2.6e-38,130.2,27.3,2,93,265,362,264,362,0.97,Calcitonin gene-related peptide regulator C terminal# noqa
2,RRM_1,PF00076.17,70,sp|O43251|RFOX2_HUMAN,-,390,8e-19,67.0,0.1,1,1,5.9e-22,1.7000000000000002e-18,65.9,0.1,2,70,124,191,123,191,0.97,"RNA recognition motif. (a.k.a. RRM, RBD, or RNP domain)"
3,RRM_6,PF14259.1,70,sp|O43251|RFOX2_HUMAN,-,390,2.4e-15,56.2,0.1,1,1,1.3999999999999999e-18,4.3e-15,55.4,0.1,1,70,123,191,123,191,0.95,"RNA recognition motif (a.k.a. RRM, RBD, or RNP domain)"
4,RRM_5,PF13893.1,56,sp|O43251|RFOX2_HUMAN,-,390,8.099999999999999e-11,41.6,0.1,1,1,5.9e-14,1.8000000000000002e-10,40.5,0.1,1,54,137,193,137,195,0.9,"RNA recognition motif. (a.k.a. RRM, RBD, or RNP domain)"
5,RRM_3,PF08777.6,105,sp|O43251|RFOX2_HUMAN,-,390,0.084,12.7,0.0,1,1,6.7e-05,0.2,11.5,0.0,17,79,136,202,127,206,0.83,RNA binding motif# noqa
"""
    true = pd.read_csv(StringIO(true_s), index_col=0, comment='#')
    pdt.assert_frame_equal(test, true)
