import logging
import sys
import warnings

from ..common import ILLEGAL_JUNCTIONS
from ..util import timestamp


with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import pandas as pd

logging.basicConfig()

idx = pd.IndexSlice


def filter_and_sum(reads, min_reads, junctions, debug=False):
    """Require minimum reads and sum junctions from the same sample

    Remove all samples that don't have enough reads

    """
    logger = logging.getLogger('outrigger.psi.filter_and_sum')
    if debug:
        logger.setLevel(10)

    if reads.empty:
        return reads

    # Remove all samples that don't have enough reads in all required junctions
    reads = reads.groupby(level=1).filter(
        lambda x: (x >= min_reads).all() and len(x) == len(junctions))
    logger.debug('filtered reads:\n' + repr(reads.head()))

    # Sum all reads from junctions in the same samples (level=1), remove NAs
    reads = reads.groupby(level=1).sum().dropna()
    logger.debug('summed reads:\n' + repr(reads.head()))

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
    junction_names = splice_junction_reads.index.levels[0]

    # Get junctions that can't exist for this to be a valid event
    illegal_junctions = junction_locations[ILLEGAL_JUNCTIONS]
    illegal_samples = []
    if isinstance(illegal_junctions, str):
        illegal_junctions = illegal_junctions.split('|')

        # Get samples that have illegal junctions
        if junction_names.isin(illegal_junctions).sum() > 0:
            illegal_reads = splice_junction_reads.loc[
                idx[illegal_junctions, :], reads_col]
            illegal_samples = illegal_reads.index.get_level_values('sample_id')

    if junction_names.isin(junctions).sum() > 0:
        reads_subset = splice_junction_reads.loc[idx[junctions, :], reads_col]
        legal_samples = reads_subset.index.get_level_values('sample_id')
        legal_samples = legal_samples.difference(illegal_samples)
        return reads_subset.loc[idx[:, legal_samples]]
    else:
        return pd.Series()


def calculate_psi(event_annotation, splice_junction_reads,
                  isoform1_junctions, isoform2_junctions, reads_col='reads',
                  min_reads=10, debug=False):
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
        correspond to isoform1, the Psi=0 isoform, e.g. ['junction13'] for SE
        (junctions between exon1 and exon3)
    isoform2_junctions : list
        Columns in `event_annotation` which represent junctions that
        correspond to isoform2, the Psi=1 isoform, e.g.
        ['junction12', 'junction23'] (junctions between exon1, exon2, and
        junction between exon2 and exon3)
    event_col : str
        Column in `event_annotation` which is a unique identifier for each
        row, e.g.

    Returns
    -------
    psi : pandas.DataFrame
        An (samples, events) dataframe of the percent spliced-in values
    """
    log = logging.getLogger('outrigger.psi.calculate_psi')

    if debug:
        log.setLevel(10)

    junction_cols = isoform1_junctions + isoform2_junctions

    # There are multiple rows with the same event id because the junctions
    # are the same, but the flanking exons may be a little wider or shorter,
    # but ultimately the event Psi is calculated only on the junctions so the
    # flanking exons don't matter for this. But, all the exons are in
    # exon\d.bed in the index! And you, the lovely user, can decide what you
    # want to do with them!
    grouped = event_annotation.groupby(level=0, axis=0)

    n_events = len(grouped.size())

    sys.stdout.write('{}\t\tIterating over {} events ...\n'.format(
        timestamp(), n_events))

    psi_df = pd.DataFrame(index=splice_junction_reads.index.levels[1],
                          columns=sorted(grouped.groups.keys()))

    for i, (event_id, event_df) in enumerate(grouped):
        if (i+1) % 1000 == 0:
            sys.stdout.write('{}\t\t\t{} events completed\n'.format(
                timestamp(), i))
        junction_locations = event_df.iloc[0]

        isoform1 = maybe_get_isoform_reads(splice_junction_reads,
                                           junction_locations,
                                           isoform1_junctions, reads_col)
        isoform2 = maybe_get_isoform_reads(splice_junction_reads,
                                           junction_locations,
                                           isoform2_junctions, reads_col)

        log.debug('--- junction columns of event ---\n%s',
                  repr(junction_locations[junction_cols]))
        log.debug('--- isoform1 ---\n%s', repr(isoform1))
        log.debug('--- isoform2 ---\n%s', repr(isoform2))

        isoform1 = filter_and_sum(isoform1, min_reads, isoform1_junctions,
                                  debug=debug)
        isoform2 = filter_and_sum(isoform2, min_reads, isoform2_junctions,
                                  debug=debug)

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
        psi = isoform2/(isoform2 + multiplier * isoform1)
        log.debug('--- Psi ---\n%s', repr(psi))
        if not psi.empty:
            psi.name = event_id
            psi_df[event_id] = psi
    sys.stdout.write('{}\t\t\tDone.\n'.format(timestamp()))

    return psi_df
