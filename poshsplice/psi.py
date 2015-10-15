import six
import pandas as pd

idx = pd.IndexSlice

MIN_READS = 10

def filter_and_sum(reads, min_reads, junctions):
    """Require minimum reads and sum junctions from the same sample

    Remove all samples that don't have enough reads

    """
    if reads.empty:
        return reads
    reads = reads.groupby(level=1).filter(lambda x: all(x >= min_reads) and
                                                    len(x) == len(junctions))
    return reads.groupby(level=1).sum().dropna()

def maybe_get_isoform_reads(splice_junction_reads, row, junctions, reads_col):
    if splice_junction_reads.index.levels[0].isin(row[junctions]).sum() > 0:
        return splice_junction_reads.loc[idx[row[junctions], :], reads_col]
    else:
        return pd.Series()



def calculate_psi(exons_to_junctions, splice_junction_reads,
                  isoform1_junctions, isoform2_junctions,
                  illegal_junctions=None, reads_col='reads',
                  min_reads=10, event_col='event_id', debug=False):
    """Compute percent-spliced-in of events based on junction reads

    Parameters
    ----------
    exons_to_junctions : pandas.DataFrame
        A table where each row represents a single splicing event. The required
        columns are the ones specified in `isoform1_junctions`,
        `isoform2_junctions`, and `event_col`.
    splice_junction_reads : pandas.DataFrame
        A table of the number of junction reads observed in each sample.
        *Important:* This must be multi-indexed by the junction and sample id
    reads_col : str, optional (default "reads")
        The name of the column in `splice_junction_reads` which represents the
        number of reads observed at a splice junction of a particular sample.
    min_reads : int, optional (default 10)
        The minimum number of reads that need to be observed at each junction
        for an event to be counted.
    isoform1_junctions : list
        Columns in `exons_to_junctions` which represent junctions that
        correspond to isoform1, the Psi=0 isoform
    isoform2_junctions : list
        Columns in `exons_to_junctions` which represent junctions that
        correspond to isoform2, the Psi=1 isoform
    event_col : str
        Column in `exons_to_junctions` which is a unique identifier for each
        row, e.g.

    Returns
    -------
    psi : pandas.DataFrame
        An (samples, events) dataframe of the percent spliced-in values
    """
    psis = []

    for i, row in exons_to_junctions.iterrows():
        isoform1 = maybe_get_isoform_reads(splice_junction_reads, row,
                                           isoform2_junctions, reads_col)
        isoform2 = maybe_get_isoform_reads(splice_junction_reads, row,
                                           isoform1_junctions, reads_col)

        if debug:
            six.print_('\n\n', row.junction12, row.junction23, row.junction13)
            six.print_('--- isoform1 ---\n', isoform1)
            six.print_('--- isoform2 ---\n', isoform2)

        isoform1 = filter_and_sum(isoform1, min_reads, isoform1_junctions)
        isoform2 = filter_and_sum(isoform2, min_reads, isoform2_junctions)

        if isoform1.empty and isoform2.empty:
            # If both are empty after filtering this event --> don't calculate
            continue

        if debug:
            six.print_('\n- After filter and sum -')
            six.print_('--- isoform1 ---\n', isoform1)
            six.print_('--- isoform2 ---\n', isoform2)


        isoform1, isoform2 = isoform1.align(isoform2, 'outer')

        isoform1 = isoform1.fillna(0)
        isoform2 = isoform2.fillna(0)

        multiplier = float(len(isoform1_junctions))/len(isoform2_junctions)
        psi = (isoform2)/(isoform2 + multiplier * isoform1)
        if debug:
            six.print_('--- Psi ---\n', psi)
        psi.name = row[event_col]
        psis.append(psi)
    if len(psis) > 0:
        psi_df = pd.concat(psis, axis=1)
    else:
        psi_df = pd.DataFrame(index=splice_junction_reads.index.levels[1])
    return psi_df
