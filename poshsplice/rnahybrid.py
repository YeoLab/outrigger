__author__ = 'olgabotvinnik'

from collections import Counter

import pandas as pd


def read_rnahybrid_output(rnahybrid_output):
    """Read output from RNAHybrid algorithm

    :param rnahybrid_output:
    :return:
    """
    return pd.read_csv(rnahybrid_output,
                       sep=':', header=None,
                       names=['chrom', 'start-stop', 'exon_length', 'mirna', 'mirna_length',
                              'minimum_free_energy', 'p_value', 'target_bound_start', 'mirna_unbound',
                              'mirna_bound', 'exon_bound', 'exon_bound'])
