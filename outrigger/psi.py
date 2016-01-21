import logging

import pandas as pd

idx = pd.IndexSlice
MIN_READS = 10

SE_ISOFORM1_JUNCTIONS = ['junction13']
SE_ISOFORM2_JUNCTIONS = ['junction12', 'junction23']

MXE_ISOFORM1_JUNCTIONS = ['junction13', 'junction34']
MXE_ISOFORM2_JUNCTIONS = ['junction12', 'junction24']

ISOFORM_JUNCTIONS = {'se': {'isoform1_junctions': SE_ISOFORM1_JUNCTIONS,
                       'isoform2_junctions': SE_ISOFORM2_JUNCTIONS},
                'mxe': {'isoform1_junctions': MXE_ISOFORM1_JUNCTIONS,
                        'isoform2_junctions': MXE_ISOFORM2_JUNCTIONS}}


def filter_and_sum(reads, min_reads, junctions):
    """Require minimum reads and sum junctions from the same sample

    Remove all samples that don't have enough reads

    """
    if reads.empty:
        return reads

    # Remove all samples that don't have enough reads in all required junctions
    reads = reads.groupby(level=1).filter(
        lambda x: all(x >= min_reads) and len(x) == len(junctions))

    # Sum all reads from junctions in the same samples (level=1), remove NAs
    reads = reads.groupby(level=1).sum().dropna()

    return reads


def maybe_get_isoform_reads(splice_junction_reads, junction_locations,
                            isoform_junctions, reads_col):
    """If there are junction reads at a junction_locations, return the reads

    Parameters
    ----------
    splice_junction_reads : pandas.DataFrame
        A table of the number of junction reads observed in each sample.
        *Important:* This must be multi-indexed by the junction and sample id
    junction_locations : pandas.Series
        Mapping of junction names (from var:`isoform_junctions`) to their
        chromosomal locations
    isoform_junctions : list of str
        Names of the isoform_junctions in question, e.g.
        ['junction12', 'junction13']
    reads_col : str
        Name of the column containing the number of reads in
        `splice_junction_reads`

    Returns
    -------
    reads : pandas.Series
        Number of reads at each junction for this isoform
    """
    junctions = junction_locations[isoform_junctions]
    if splice_junction_reads.index.levels[0].isin(junctions).sum() > 0:
        return splice_junction_reads.loc[idx[junctions, :], reads_col]
    else:
        return pd.Series()


def calculate_psi(event_annotation, splice_junction_reads,
                  isoform1_junctions, isoform2_junctions, reads_col='reads',
                  min_reads=10, event_col='event_id', debug=False):
    """Compute percent-spliced-in of events based on junction reads

    Parameters
    ----------
    event_annotation : pandas.DataFrame
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
        Columns in `event_annotation` which represent junctions that
        correspond to isoform1, the Psi=0 isoform, e.g. ['junction13']
    isoform2_junctions : list
        Columns in `event_annotation` which represent junctions that
        correspond to isoform2, the Psi=1 isoform, e.g.
        ['junction12', 'junction23']
    event_col : str
        Column in `event_annotation` which is a unique identifier for each
        row, e.g.

    Returns
    -------
    psi : pandas.DataFrame
        An (samples, events) dataframe of the percent spliced-in values
    """
    log = logging.getLogger('calculate_psi')

    psis = []

    junction_cols = isoform1_junctions + isoform2_junctions

    for i, row in event_annotation.iterrows():
        isoform1 = maybe_get_isoform_reads(splice_junction_reads, row,
                                           isoform1_junctions, reads_col)
        isoform2 = maybe_get_isoform_reads(splice_junction_reads, row,
                                           isoform2_junctions, reads_col)

        log.debug('\n\n%s\t%s\t%s', row[junction_cols])
        log.debug('--- isoform1 ---\n%s', repr(isoform1))
        log.debug('--- isoform2 ---\n%s', repr(isoform2))

        isoform1 = filter_and_sum(isoform1, min_reads, isoform1_junctions)
        isoform2 = filter_and_sum(isoform2, min_reads, isoform2_junctions)

        if isoform1.empty and isoform2.empty:
            # If both are empty after filtering this event --> don't calculate
            continue

        log.debug('\n- After filter and sum -')
        log.debug('--- isoform1 ---\n%s', repr(isoform1))
        log.debug('--- isoform2 ---\n%s', repr(isoform2))

        isoform1, isoform2 = isoform1.align(isoform2, 'outer')

        isoform1 = isoform1.fillna(0)
        isoform2 = isoform2.fillna(0)

        multiplier = float(len(isoform2_junctions))/len(isoform1_junctions)
        psi = (isoform2)/(isoform2 + multiplier * isoform1)
        log.debug('--- Psi ---\n%.2f', psi)
        psi.name = row[event_col]
        psis.append(psi)

    if len(psis) > 0:
        psi_df = pd.concat(psis, axis=1)
    else:
        psi_df = pd.DataFrame(index=splice_junction_reads.index.levels[1])
    return psi_df
