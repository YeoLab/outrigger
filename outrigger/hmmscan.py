import pandas as pd


def read_hmmscan(hmmscan_out):
    """Read output from hmmscan

    Parameters
    ----------
    hmmscan_out : str
        Filename of hmmscan output

    Returns
    -------
    hmmscan_df : pandas.DataFrame
        Parsed string of the hmmscan output
    """
    entries = []
    with open(hmmscan_out) as f:
        for line in f.readlines():
            if line.startswith('#'):
                continue
            split = line.split()
            beginning = split[:22]
            end = ' '.join(split[22:])
            entries.append(beginning + [end])
    columns = ['target_name', 'target_accession', 'target_length',
               'query_name', 'query_accession', 'query_length',
               'sequence_e_value', 'sequence_score', 'sequence_bias',
               'domain_number', 'domain_total', 'domain_conditional_e_value',
               'domain_independent_e_value', 'domain_score', 'domain_bias',
               'target_start', 'target_stop', 'query_start', 'query_stop',
               'query_domain_envelope_start', 'query_domain_envelope_stop',
               'mean_posterior_probability', 'target_description']
    df = pd.DataFrame.from_records(entries, columns=columns)
    df = df.convert_objects(convert_numeric=True)
    return df
