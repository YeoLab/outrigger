import logging

import joblib
import pandas as pd

from ..common import ILLEGAL_JUNCTIONS, MIN_READS, READS
from ..util import progress, done


logging.basicConfig()

idx = pd.IndexSlice


def _scale(x, n_junctions, method='mean', min_reads=MIN_READS):
    multiplier = -1 if (x < min_reads).any() else 1
    if method == 'mean':
        return multiplier * x.sum()/float(n_junctions)
    elif method == 'min':
        return multiplier * x.min()


def _filter_and_scale(reads, n_junctions, debug=False, min_reads=MIN_READS,
                      method='mean'):
    """Remove samples without reads on all junctions and sum across junctions

    If any junction on the isoform has fewer than the minimum reads, flag it
    as -1

    Parameters
    ----------

    method : 'mean' | 'min'
        If there are more than 2 junctions for this isoform, then they have to
        be consolidated somehow.
        - "mean": Sum the reads and divide by the number of junctions
        - "min": Use the minimum number of reads
    """
    logger = logging.getLogger('outrigger.psi.filter_and_scale')
    if debug:
        logger.setLevel(10)

    if reads.empty:
        return reads

    if n_junctions > 1:
        # Remove all samples that don't have all required junctions
        reads = reads.groupby(level=1).filter(
            lambda x: len(x) == n_junctions)
        logger.debug('filtered reads:\n' + repr(reads.head()))

    reads = reads.groupby(level=1).apply(
        lambda x: _scale(x, n_junctions, method, min_reads))

    # Sum all reads from junctions in the same samples (level=1), remove NAs
    # reads = reads.groupby(level=1).sum().dropna()
    logger.debug('summed reads:\n' + repr(reads.head()))

    return reads


def _maybe_get_isoform_reads(splice_junction_reads, junction_locations,
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


def _remove_insufficient_reads(isoform1, isoform2):
    """Exclude samples whose junctions have insufficient reads (marked with -1)

    In filter_and_sum, if any junction had fewer than expected reads, set it
     to -1

    ONLY Use junctions with individually sufficient reads
    """
    isoform1, isoform2 = isoform1.align(isoform2, 'outer')

    isoform1 = isoform1.fillna(0)
    isoform2 = isoform2.fillna(0)

    invalid = (isoform1 < 0) | (isoform2 < 0)

    isoform1 = isoform1[~invalid]
    isoform2 = isoform2[~invalid]
    return isoform1, isoform2


def _maybe_sufficient_reads(isoform1, isoform2, n_junctions, min_reads,
                            case=1):
    """Check if the sum of reads is enough compared to number of junctions"""
    if (isoform1.sum() + isoform2.sum()) >= (min_reads * n_junctions):
        # Case 5a: There are sufficient junction reads
        return isoform1, isoform2, 'Case {case}a: There are sufficient ' \
                                   'junction reads'.format(case=case)
    else:
        # Case 5b: There are insufficient junction reads
        return None, None, 'Case {case}b: There are insufficient ' \
                           'junction reads'.format(case=case)


def _maybe_reject_events(isoform1, isoform2, n_junctions, min_reads=MIN_READS):
    """Given the junction reads of isoform1 and isoform2, remove them if they are bad"""

    if (isoform1 >= min_reads).all() and (isoform2 == 0).all():
        # Case 1: Perfect exclusion
        return isoform1, isoform2, 'Case 1: Perfect exclusion'
    elif (isoform1 == 0).all() and (isoform2 >= min_reads).all():
        # Case 2: Perfect inclusion
        return isoform1, isoform2, 'Case 2: Perfect inclusion'
    elif (isoform1 >= min_reads).all() and (isoform2 >= min_reads).all():
        # Case 3: coverage on both isoforms
        return isoform1, isoform2, 'Case 3: coverage on both isoforms'
    elif (isoform1 == 0).any() or (isoform2 == 0).any():
        # Case 4: Any observed junction is zero and it's not all of one isoform
        return None, None, "Case 4: Any observed junction is zero and it's not all of one isoform"
    elif (isoform1 >= min_reads).all() and (isoform2 < min_reads).all():
        # Case 5: isoform1 totally covered and isoform2 not
        if (isoform1.sum() + isoform2.sum()) >= (min_reads * n_junctions):
            # Case 5a: There are sufficient junction reads
            return isoform1, isoform2, 'Case 5a: There are sufficient junction reads'
        else:
            # Case 5b: There are insufficient junction reads
            return None, None, 'Case 5b: There are insufficient junction reads'
    elif (isoform1 < min_reads).all() and (isoform2 >= min_reads).all():
        # Case 6: Isoform2 is totally covered and isoform1 is not
        if (isoform1.sum() + isoform2.sum()) >= (min_reads * n_junctions):
            # Case 6a: There are sufficient junction reads
            return isoform1, isoform2, 'Case 6a: There are sufficient junction reads'
        else:
            # Case 6b: There are insufficient junction reads
            return None, None, 'Case 6b: There are insufficient junction reads'
    elif (isoform1 >= min_reads).all() and (isoform2 < min_reads).any():
        # Case 7: Isoform 1 is fully covered and isoform2 is questionable
        if (isoform1.sum() + isoform2.sum()) >= (min_reads * n_junctions):
            # Case 7a: There are sufficient junction reads
            return isoform1, isoform2, 'Case 7a: There are sufficient junction reads'
        else:
            # Case 7b: There are insufficient junction reads
            return None, None, 'Case 7b: There are insufficient junction reads'
    elif (isoform1 < min_reads).any() and (isoform2 >= min_reads).all():
        # Case 8: Isoform 1 is fully covered and isoform2 is questionable
        if (isoform1.sum() + isoform2.sum()) >= (min_reads * n_junctions):
            # Case 8a: There are sufficient junction reads
            return isoform1, isoform2, '8a: There are sufficient junction reads'
        else:
            # Case 8b: There are insufficient junction reads
            return None, None, 'Case 8b: There are insufficient junction reads'
    if (isoform1 < min_reads).any() or (isoform2 < min_reads).any():
        # Case 9: insufficient reads somehow
        if (isoform1 < min_reads).all() and (isoform2 < min_reads).any():
            # Case 9a: 3 junctions have less than minimum reads (2 on iso1 and 1 on iso2)
            return None, None, 'Case 9a: 3 junctions have less than minimum reads (2 on iso1 and 1 on iso2)'
        if (isoform1 < min_reads).any() and (isoform2 < min_reads).all():
            # Case 9b: 3 junctions have less than minimum reads (2 on iso2 and one on iso1)
            return None, None, 'Case 9b: 3 junctions have less than minimum reads (2 on iso2 and one on iso1)'

        if (isoform1.sum() + isoform2.sum()) >= (min_reads * n_junctions):
            # Case 9c: There are sufficient junction reads
            return isoform1, isoform2, 'Case 9c: There are sufficient junction reads'
        else:
            # Case 9d: There are insufficient junction reads
            return None, None, 'Case 9d: There are insufficient junction reads'

    elif (isoform1 < min_reads).any() or (isoform2 < min_reads).any():
        # Case 9: isoform1 and isoform2 don't have sufficient reads
        return None, None, "Case 9: isoform1 and isoform2 don't have sufficient reads"

    # If none of these is true, then there's some uncaught case
    return '???', '???', "Case ???"


def _single_event_psi(event_id, event_df, splice_junction_reads,
                      isoform1_junctions, isoform2_junctions, reads_col=READS,
                      min_reads=MIN_READS, method='mean', debug=False,
                      log=None):
    """Calculate percent spliced in for a single event across all samples

    Returns
    -------
    psi : pandas.Series or None
        If unable to calculate psi for this event due to insufficient junctions
        then None is returned.
    """
    junction_locations = event_df.iloc[0]

    n_junctions1 = len(isoform1_junctions)
    n_junctions2 = len(isoform2_junctions)

    isoform1 = _maybe_get_isoform_reads(splice_junction_reads,
                                        junction_locations,
                                        isoform1_junctions, reads_col)
    isoform2 = _maybe_get_isoform_reads(splice_junction_reads,
                                        junction_locations,
                                        isoform2_junctions, reads_col)

    if debug and log is not None:
        junction_cols = isoform1_junctions + isoform2_junctions
        log.debug('--- junction columns of event ---\n%s',
                  repr(junction_locations[junction_cols]))
        log.debug('--- isoform1 ---\n%s', repr(isoform1))
        log.debug('--- isoform2 ---\n%s', repr(isoform2))

    isoform1 = _filter_and_scale(isoform1, n_junctions1, min_reads, method)
    isoform2 = _filter_and_scale(isoform2, n_junctions2, min_reads, method)

    if isoform1.empty and isoform2.empty:
        # If both are empty after filtering this event --> don't calculate
        return

    if debug and log is not None:
        log.debug('\n- After filter and sum -')
        log.debug('--- isoform1 ---\n%s', repr(isoform1))
        log.debug('--- isoform2 ---\n%s', repr(isoform2))

    isoform1, isoform2 = _remove_insufficient_reads(isoform1, isoform2)

    if debug and log is not None:
        log.debug('\n- After removing insufficient reads -')
        log.debug('--- isoform1 ---\n%s', repr(isoform1))
        log.debug('--- isoform2 ---\n%s', repr(isoform2))

    psi = isoform2 / (isoform2 + isoform1)
    psi.name = event_id

    if debug and log is not None:
        log.debug('--- Psi ---\n%s', repr(psi))

    if not psi.empty:
        return psi, isoform1, isoform2


def _maybe_parallelize_psi(event_annotation, splice_junction_reads,
                           isoform1_junctions, isoform2_junctions,
                           reads_col=READS, min_reads=MIN_READS, method='mean',
                           n_jobs=-1, debug=False, log=None):
    # There are multiple rows with the same event id because the junctions
    # are the same, but the flanking exons may be a little wider or shorter,
    # but ultimately the event Psi is calculated only on the junctions so the
    # flanking exons don't matter for this. But, all the exons are in
    # exon\d.bed in the index! And you, the lovely user, can decide what you
    # want to do with them!
    grouped = event_annotation.groupby(level=0, axis=0)

    n_events = len(grouped.size())

    if n_jobs == 1:
        progress('\tIterating over {} events ...\n'.format(n_events))
        psis = []
        isoform1s = []
        isoform2s = []
        for event_id, event_df in grouped:
            psi, isoform1, isoform2 = _single_event_psi(
                event_id, event_df, splice_junction_reads,
                isoform1_junctions, isoform2_junctions,
                reads_col, min_reads, method, debug, log)
            psis.append(psi)
            isoform1s.append(isoform1)
            isoform2s.append(isoform2)
    else:
        processors = n_jobs if n_jobs > 0 else joblib.cpu_count()
        progress("\tParallelizing {} events' Psi calculation across {} "
                 "CPUs ...\n".format(n_events, processors))
        outputs = joblib.Parallel(n_jobs=n_jobs)(
            joblib.delayed(_single_event_psi)(
                event_id, event_df, splice_junction_reads,
                isoform1_junctions, isoform2_junctions, reads_col,
                min_reads, method) for event_id, event_df in grouped)
        import pdb; pdb.set_trace()

    return psis


def calculate_psi(event_annotation, splice_junction_reads,
                  isoform1_junctions, isoform2_junctions, reads_col=READS,
                  min_reads=MIN_READS, method='mean', n_jobs=-1,
                  debug=False):
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

    psis = _maybe_parallelize_psi(event_annotation, splice_junction_reads,
                                  isoform1_junctions, isoform2_junctions,
                                  reads_col, min_reads, method, n_jobs, debug,
                                  log)

    # use only non-empty psi outputs
    psis = filter(lambda x: x is not None, psis)
    try:
        psi_df = pd.concat(psis, axis=1)
    except ValueError:
        psi_df = pd.DataFrame(index=splice_junction_reads.index.levels[1])
    done(n_tabs=3)
    psi_df.index.name = splice_junction_reads.index.names[1]
    return psi_df
