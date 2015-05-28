import pandas as pd

__author__ = 'olgabotvinnik'

RNAHYBRID_COLUMNS = ['chrom', 'start-stop', 'exon_length', 'mirna',
                     'mirna_length', 'minimum_free_energy', 'p_value',
                     'target_bound_start', 'mirna_unbound',
                     'mirna_bound', 'exon_bound', 'exon_bound']


def read_rnahybrid_output(rnahybrid_output):
    """Read output from RNAHybrid algorithm

    Parameters
    ----------
    rnahybrid_output : str
        Filename of the rnahybrid output

    Returns
    -------
    rnahybrid_df : pandas.DataFrame
        A pandas dataframe of the RNAHybrid output
    """
    return pd.read_csv(rnahybrid_output, sep=':', header=None,
                       names=RNAHYBRID_COLUMNS)
