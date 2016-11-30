import logging

import joblib
import pandas as pd

from ..common import ILLEGAL_JUNCTIONS, MIN_READS, READS, \
    UNEVEN_COVERAGE_MULTIPLIER, SAMPLE_ID, EVENT_ID, NOTES, PSI, \
    JUNCTION_ID
from ..util import progress, done


logging.basicConfig()

idx = pd.IndexSlice


def _scale(x, n_junctions, method='mean'):
    if method == 'mean':
        return x.sum()/float(n_junctions)
    elif method == 'min':
        return x.min()


def _single_sample_maybe_sufficient_reads(isoform1, isoform2, n_junctions,
                                          min_reads, case, letters='ab'):
    """Check if the sum of reads is enough compared to number of junctions

    Parameters
    ----------
    isoform1, isoform2 : pandas.Series
        Number of reads found for a single sample's isoform1 and isoform2
        junctions
    n_junctions : int
        Total number of junctions, also len(isoform1) + len(isoform2)
    min_reads : int
        Minimum number of reads per junction
    case : str
        English explanation with case number for why an event was rejected or
        not
    letters : str, optional
        Cases could have multiple clauses, by default the case is "option a" or
         "option b" but other letters could be used. (default='ab')

    Returns
    -------
    isoform1, isoform2 : pandas.Series or None
        If the sum of reads on all isoforms is equal to or greater than the
        minimum number of reads times the number of junctions, then return the
        event. Otherwise, there are insufficient reads and the event is
        rejected.
    """
    if (isoform1.sum() + isoform2.sum()) >= (min_reads * n_junctions):
        # Case 5a: There are sufficient junction reads
        return isoform1, isoform2, '{case}, option {letter}: There are ' \
                                   'sufficient junction ' \
                                   'reads'.format(case=case, letter=letters[0])
    else:
        # Case 5b: There are insufficient junction reads
        return None, None, '{case}, option {letter}: There are insufficient ' \
                           'junction reads'.format(case=case,
                                                   letter=letters[1])


def _single_sample_check_unequal_read_coverage(
    isoform, uneven_coverage_multiplier=UNEVEN_COVERAGE_MULTIPLIER):
    """If one junction of an isoform is more heavily covered, reject it

    If the difference in read depth between two junctions of an isoform is
    higher than a multiplicative amount, reject the isoform.

    If the isoform only has one junction, don't reject it.

    Possible fallacy: assumes there can be at most two junctions per isoform...

    Parameters
    ----------
    isoform : pandas.Series
        Number of exon-exon junction-spanning reads found for an isoform
    uneven_coverage_multiplier : int
        Scale factor for the maximum amount bigger one side of a junction can
        be before rejecting the event, e.g. for an SE event with two junctions,
        junction12 and junction23, junction12=40 but junction23=500, then this
        event would be rejected because 500 > 40*10

    Returns
    -------
    isoform : pandas.Series or None
        If not rejected, the original is returned, otherwise return None
    """

    if len(isoform) == 1:
        return isoform

    multiplied = isoform * uneven_coverage_multiplier

    junction0 = isoform.iloc[0]
    junction1 = isoform.iloc[1]

    if junction0 > junction1 and junction0 > multiplied.iloc[1]:
        return None
    elif junction1 > junction0 and junction1 > multiplied.iloc[0]:
        return None
    else:
        return isoform


def _maybe_reject(reads, isoform1_ids, isoform2_ids, illegal_ids,
                  n_junctions, min_reads=MIN_READS,
                  uneven_coverage_multiplier=UNEVEN_COVERAGE_MULTIPLIER):
    """Remove samples with reads that are incompatible with event definition

    Parameters
    ----------
    reads : pandas.DataFrame
        A (n_samples, n_junctions) table of the number of reads found in each
        samples' exon-exon junctions
    isoform1_ids : list of str
        Column names in ``reads`` tha correspond to the junction ids that are
        contained within isoform 1, e.g. ['junction1:chr1:100-400:+']
    isoform2_ids :
        Column names in ``reads`` tha correspond to the junction ids that are
        contained within isoform 2, e.g. ['junction:chr1:100-200:+',
        'junction:chr1:300-400:+']
    illegal_ids : list of str
        Column names in ``reads`` tha correspond to the junction ids that are
        contained within junctions that are not compatible with the event
        definition
    n_junctions : int
        Total number of legal junctions,
        i.e. len(isoform1_ids) + len(isoform2_ids)
    min_reads : int
        Minimum number of reads for a junction to be viable. The rules
        governing compatibility of events are complex, and it is recommended to
        read the documentation for ``outrigger psi``
    uneven_coverage_multiplier : int
        Scale factor for the maximum amount bigger one side of a junction can
        be before rejecting the event, e.g. for an SE event with two junctions,
        junction12 and junction23, junction12=40 but junction23=500, then this
        event would be rejected because 500 > 40*10

    Returns
    -------

    """
    if not isinstance(illegal_ids, float):
        junctions_with_illegal_coverage = reads[illegal_ids] >= min_reads
        samples_with_illegal_coverage = junctions_with_illegal_coverage.any(
            axis=1)
        reads = reads.loc[~samples_with_illegal_coverage]

    maybe_rejected = reads.apply(
        lambda sample: _single_maybe_reject(
            sample, isoform1_ids, isoform2_ids,
            n_junctions=n_junctions, min_reads=min_reads,
            uneven_coverage_multiplier=uneven_coverage_multiplier), axis=1)
    return maybe_rejected


def _single_maybe_reject(
    sample, isoform1_ids, isoform2_ids, n_junctions, min_reads=MIN_READS,
    uneven_coverage_multiplier=UNEVEN_COVERAGE_MULTIPLIER):
    """Given a row of junction reads, return a filtered row of reads

    For a single sample's junction reads of an isoform, check if they should be
    rejected, and if they are, return a row with all NAs for the reads. Always
    include the case by which the reads were or were not rejected

    Parameters
    ----------
    sample : pandas.Series
        A single sample's junction reads across all isoforms of an alternative
        event
    isoform1_ids : list
        List of strings that correspond to indicies in ``sample``
    isoform2_ids : list
        List of strings that correspond to indicies in ``sample``
    n_junctions : int
        Total number of junctions expected in the splicing event
    min_reads : int, optional
        Minimum number of junction reads for an event to be valid. See
        documentation for much more detailed information regarding when events
        are rejected or retained. (default=10)
    uneven_coverage_multiplier : int, optional
        When checking for uneven coverage between two sides of a junction, one
        side must be this amount bigger than the other side to be rejected.
        For example, if one side has 10x (default) more read coverage than the
        other, then reject the event. (default=10)

    Returns
    -------
    single_maybe_rejected : pandas.Series
        Unrejected reads of a single samples' splicing event. If the event was
        rejected, the reads are replaced with NAs. This series has one more
        field than the input "sample", with the field of "notes" that explains
        why this event was or was not rejected.
    """
    isoform1, isoform2, case = _single_isoform_maybe_reject(
            sample[isoform1_ids], sample[isoform2_ids],
            n_junctions=n_junctions, min_reads=min_reads,
            uneven_coverage_multiplier=uneven_coverage_multiplier)
    if isoform1 is None:
        single_maybe_rejected = pd.Series(None, index=sample.index)
    else:
        single_maybe_rejected = sample
    single_maybe_rejected[NOTES] = case
    return single_maybe_rejected



def _single_isoform_maybe_reject(
    isoform1, isoform2, n_junctions, min_reads=MIN_READS,
    uneven_coverage_multiplier=UNEVEN_COVERAGE_MULTIPLIER):
    """Given junction reads of isoform1 and isoform2, remove if they are bad

    Parameters
    ----------
    isoform1, isoform2 : pandas.Series
        Number of reads found on exon-exon junctions for isoform1 and isoform2
    n_junctions : int
        Total number of junctions. Could also be found by
        len(isoform1) + len(isoform2)
    min_reads : int, optional
        Minimum number of reads for a junction to be counted, though the full
        explanation is a little more complicated, please see the documentation
        for more details. (default=10)
    uneven_coverage_multiplier : int, optional
        Scale factor for the maximum amount bigger one side of a junction can
        be before rejecting the event, e.g. for an SE event with two junctions,
        junction12 and junction23, junction12=40 but junction23=500, then this
        event would be rejected because 500 > 40*10

    Returns
    -------
    isoform1, isoform2 : pandas.Series or None
        If the event was not rejected, return the original event, otherwise
        return None
    case : str
        Reason for rejecting or retaining the event
    """

    isoform1 = _single_sample_check_unequal_read_coverage(
        isoform1, uneven_coverage_multiplier)
    isoform2 = _single_sample_check_unequal_read_coverage(
        isoform2, uneven_coverage_multiplier)

    if isoform1 is None or isoform2 is None:
        # Case 1: Unbalanced number of reads between two sides of an isoform
        return None, None, "Case 1: Unequal read coverage"
    elif (isoform1 == 0).all() and (isoform2 == 0).all():
        # Case 2: All reads are zero
        return None, None, 'Case 2: No observed reads'
    elif (isoform1 >= min_reads).all() and (isoform2 == 0).all():
        # Case 3: Perfect exclusion
        return isoform1, isoform2, 'Case 3: Perfect exclusion'
    elif (isoform1 == 0).all() and (isoform2 >= min_reads).all():
        # Case 4: Perfect inclusion
        return isoform1, isoform2, 'Case 4: Perfect inclusion'
    elif (isoform1 >= min_reads).all() and (isoform2 >= min_reads).all():
        # Case 5: Sufficient coverage on both isoforms
        return isoform1, isoform2, 'Case 5: Sufficient coverage on both ' \
                                   'isoforms'
    elif (isoform1 == 0).any() or (isoform2 == 0).any():
        # Case 6: Any observed junction is zero and it's not all of one isoform
        return None, None, "Case 6: Any observed junction is zero and it's " \
                           "not all of one isoform"
    elif (isoform1 >= min_reads).all() and (isoform2 < min_reads).all():
        # Case 7: Isoform1 totally covered and isoform2 not
        return _single_sample_maybe_sufficient_reads(
            isoform1, isoform2, n_junctions, min_reads,
            'Case 7: Isoform1 totally covered and isoform2 not')
    elif (isoform1 < min_reads).all() and (isoform2 >= min_reads).all():
        # Case 8: Isoform2 is totally covered and isoform1 is not
        return _single_sample_maybe_sufficient_reads(
            isoform1, isoform2, n_junctions, min_reads,
            'Case 8: Isoform2 is totally covered and isoform1 is not')
    elif (isoform1 >= min_reads).all() and (isoform2 < min_reads).any():
        # Case 9: Isoform 1 is fully covered and isoform2 is questionable
        return _single_sample_maybe_sufficient_reads(
            isoform1, isoform2, n_junctions, min_reads,
            'Case 9: Isoform 1 is fully covered and isoform2 is questionable')
    elif (isoform1 < min_reads).any() and (isoform2 >= min_reads).all():
        # Case 10: Isoform 1 is fully covered and isoform2 is questionable
        return _single_sample_maybe_sufficient_reads(
            isoform1, isoform2, n_junctions, min_reads,
            'Case 10: Isoform 1 is fully covered and isoform2 is questionable')
    elif (isoform1 < min_reads).any() or (isoform2 < min_reads).any():
        # Case 11: insufficient reads somehow
        if (isoform1 < min_reads).all() and (isoform2 < min_reads).any():
            # Case 11a: 3 junctions have less than minimum reads (2 on iso1
            # and 1 on iso2)
            return None, None, 'Case 11a: 3 junctions have less than minimum ' \
                               'reads (2 on iso1 and 1 on iso2)'
        if (isoform1 < min_reads).any() and (isoform2 < min_reads).all():
            # Case 11b: 3 junctions have less than minimum reads (2 on iso2
            # and one on iso1)
            return None, None, 'Case 11b: 3 junctions have less than minimum ' \
                               'reads (2 on iso2 and one on iso1)'

        return _single_sample_maybe_sufficient_reads(
            isoform1, isoform2, n_junctions, min_reads,
            case='Case 11: Insufficient reads somehow', letters='cd')
    elif (isoform1 < min_reads).any() or (isoform2 < min_reads).any():
        # Case 12: isoform1 and isoform2 don't have sufficient reads
        return None, None, "Case 12: isoform1 and isoform2 don't have " \
                           "sufficient reads"

    # If none of these is true, then there's some uncaught case
    return '???', '???', "Case ???"


def _summarize_event(event_id, reads, maybe_rejected, psi,
                     isoform1_junction_ids, isoform2_junction_ids,
                     isoform1_junction_numbers, isoform2_junction_numbers):
    """Make table summarizing junction reads, psi, and notes for an event

    Parameters
    ----------
    event_id : str
        Uniquely identifying string for a splicing event
    reads : pandas.DataFrame
        A (n_samples, n_junctions) table of the number of reads found in each
        samples' exon-exon junctions for this event
    maybe_rejected : pandas.DataFrame
        A (n_samples, n_junctions + 1) table that is nearly a copy of ``reads``
        except has rejected events' junction reads replaced with NAs, and for
        all events, also has a column called "notes" which has the reason the
        event was or was not rejected
    psi : pandas.Series
        A (n_samples,) sized column of the percent spliced-in values for each
        sample, for this event
    isoform1_junction_ids : list of str
        Column names in ``reads`` tha correspond to the junction ids that are
        contained within isoform 1, e.g. ['junction1:chr1:100-400:+']
    isoform2_junction_ids : list of str
        Column names in ``reads`` tha correspond to the junction ids that are
        contained within isoform 2, e.g. ['junction:chr1:100-200:+',
        'junction:chr1:300-400:+']
    isoform1_junction_numbers : list of str
         Junction numbers corresponding to isoform 1, e.g. ['junction13']
    isoform2_junction_numbers : list of str
        Junction numbers corresponding to isoform 2, e.g. ['junction12',
        'junction23']

    Returns
    -------
    summary : pandas.DataFrame
        A (n_samples, 7) shaped table with the sample id, junction reads,
        percent spliced-in (Psi), and notes on each event in each sample, that
        explain why or why not Psi was calculated
    """
    isoform1_numbers = ['isoform1_' + x for x in isoform1_junction_numbers]
    isoform2_numbers = ['isoform2_' + x for x in isoform2_junction_numbers]

    column_renamer = dict(zip(isoform1_junction_ids,
                              isoform1_numbers))
    column_renamer.update(dict(zip(isoform2_junction_ids,
                                   isoform2_numbers)))

    summary = reads.rename(columns=column_renamer)
    summary[NOTES] = maybe_rejected[NOTES]
    summary[PSI] = psi
    summary = summary.reset_index()
    summary[EVENT_ID] = event_id
    summary.columns.name = None

    column_order = [SAMPLE_ID, EVENT_ID] + isoform1_numbers \
                   + isoform2_numbers + [PSI, NOTES]
    summary = summary[column_order]
    return summary


def _single_event_psi(event_id, event_df, junction_reads_2d,
                      isoform1_junction_numbers, isoform2_junction_numbers,
                      min_reads=MIN_READS, method='mean',
                      uneven_coverage_multiplier=UNEVEN_COVERAGE_MULTIPLIER):
    """Calculate percent spliced in for a single event across all samples

    Parameters
    ----------
    event_id : str
        Uniquely identifying string for a splicing event
    event_df : pandas.DataFrame
        A table with the event id as the index (row names) and the junction
        locations for the different isoforms. This may have multiple rows (or
        not) depending on the different widths of the flanking exons
    junction_reads_2d : pandas.DataFrame
        A (n_samples, n_total_junctions) table of the number of reads found in
        all samples' exon-exon, all junctions. Very very large, e.g.
        1000 samples x 50,000 junctions = 50 million elements
    isoform1_junction_numbers : list of str
         Junction numbers corresponding to the isoform,
         e.g. ``['junction13']``. Must be columns in ``event_df``
    isoform2_junction_numbers : list of str
        Junction numbers corresponding to the isoform, e.g. ``['junction12',
        'junction23']``. Must be columns in ``event_df``
    min_reads : int
        Minimum number of reads for a junction to be viable. The rules
        governing compatibility of events are complex, and it is recommended to
        read the documentation for ``outrigger psi``
    method : "mean" | "min"
        Denotes the method by which to aggregate junctions from the same
        isoform - either use the mean (default) or the minimum
    uneven_coverage_multiplier : int
        Scale factor for the maximum amount bigger one side of a junction can
        be before rejecting the event, e.g. for an SE event with two junctions,
        junction12 and junction23, junction12=40 but junction23=500, then this
        event would be rejected because 500 > 40*10

    Returns
    -------
    summary : pandas.DataFrame
        A (n_samples, 7) shaped table with the sample id, junction reads,
        percent spliced-in (Psi), and notes on each event in each sample, that
        explain why or why not Psi was calculated

    >>> # Example Summary output
    >>> summary.head()
                                               sample_id  junction13  \
    0  CAV_LP_Ipsi_tdTpos_cell_1-SRR2140356-GSM184094...          26
    1  CAV_LP_Ipsi_tdTpos_cell_10-SRR2140365-GSM18409...           0
    2  CAV_LP_Ipsi_tdTpos_cell_11-SRR2140366-GSM18409...          52
    3  CAV_LP_Ipsi_tdTpos_cell_12-SRR2140367-GSM18409...           0
    4  CAV_LP_Ipsi_tdTpos_cell_13-SRR2140368-GSM18409...          31

       junction12  junction23                      notes  psi  \
    0           0           0  Case 3: Perfect exclusion  0.0
    1           0           0  Case 2: No observed reads  NaN
    2           0           0  Case 3: Perfect exclusion  0.0
    3           0           0  Case 2: No observed reads  NaN
    4           0           0  Case 3: Perfect exclusion  0.0

                                                event_id
    0  isoform1=junction:chr10:128491034-128491719:-|...
    1  isoform1=junction:chr10:128491034-128491719:-|...
    2  isoform1=junction:chr10:128491034-128491719:-|...
    3  isoform1=junction:chr10:128491034-128491719:-|...
    4  isoform1=junction:chr10:128491034-128491719:-|...

    """
    junction_locations = event_df.iloc[0]

    n_junctions1 = len(isoform1_junction_numbers)
    n_junctions2 = len(isoform2_junction_numbers)
    n_junctions = n_junctions1 + n_junctions2

    isoform1_junction_ids = junction_locations[
        isoform1_junction_numbers].tolist()
    isoform2_junction_ids = junction_locations[
        isoform2_junction_numbers].tolist()
    illegal_junction_ids = junction_locations[ILLEGAL_JUNCTIONS]

    junction_cols = isoform1_junction_ids + isoform2_junction_ids

    if not isinstance(illegal_junction_ids, float):
        illegal_junction_ids = illegal_junction_ids.split('|')
        illegal_junction_ids = junction_reads_2d.columns.intersection(
            illegal_junction_ids)
        junction_cols += illegal_junction_ids

    reads = junction_reads_2d[junction_cols]

    maybe_rejected = _maybe_reject(
        reads, isoform1_junction_ids, isoform2_junction_ids,
        illegal_junction_ids, n_junctions, min_reads=min_reads,
        uneven_coverage_multiplier=uneven_coverage_multiplier)

    isoform1 = maybe_rejected[isoform1_junction_ids].apply(
        _scale, n_junctions=n_junctions1, method=method, axis=1)
    isoform2 = maybe_rejected[isoform2_junction_ids].apply(
        _scale, n_junctions=n_junctions2, method=method, axis=1)

    psi = isoform2 / (isoform2 + isoform1)

    summary = _summarize_event(event_id, reads, maybe_rejected, psi,
                               isoform1_junction_ids, isoform2_junction_ids,
                               isoform1_junction_numbers,
                               isoform2_junction_numbers)
    return summary

def _maybe_parallelize_psi(
    event_annotation, junction_reads_2d, isoform1_junctions,
    isoform2_junctions, min_reads=MIN_READS, method='mean',
    uneven_coverage_multiplier=UNEVEN_COVERAGE_MULTIPLIER, n_jobs=-1):
    """If n_jobs!=1, run the parallelized version of psi

    Parameters
    ----------
    event_annotation : pandas.DataFrame
        A table of all possible events, with event ids as the index (row names)
        and all junctions described, and contains the columns described by
        ``isoform1_junctions`` and ``isoform_junctions``
    junction_reads_2d : pandas.DataFrame
        A (n_samples, n_total_junctions) table of the number of reads found in
        all samples' exon-exon, all junctions. Very very large, e.g.
        1000 samples x 50,000 junctions = 50 million elements
    isoform1_junctions : list of str
        Junction numbers corresponding to isoform 1, e.g. ['junction13']
    isoform2_junctions : list of str
        Junction numbers corresponding to isoform 2, e.g. ['junction12',
        'junction23']
    min_reads : int, optional
        Minimum number of reads for a junction to be viable. The rules
        governing compatibility of events are complex, and it is recommended to
        read the documentation for ``outrigger psi`` (default=10)
    method : "mean" | "min", optional
        Denotes the method by which to aggregate junctions from the same
        isoform - either use the mean (default) or the minimum.
        (default="mean")
    uneven_coverage_multiplier : int, optional
        Scale factor for the maximum amount bigger one side of a junction can
        be before rejecting the event, e.g. for an SE event with two junctions,
        junction12 and junction23, junction12=40 but junction23=500, then this
        event would be rejected because 500 > 40*10 (default=10)
    n_jobs : int, optional
        Number of subprocesses to create. Default is -1, which is to use as
        many processes/cores as possible

    Returns
    -------
    summary : pandas.DataFrame
        A (n_samples * n_events, 7) shaped table with the sample id, junction
        reads, percent spliced-in (Psi), and notes on each event in each
        sample, that explains why or why not Psi was calculated
    """
    # There are multiple rows with the same event id because the junctions
    # are the same, but the flanking exons may be a little wider or shorter,
    # but ultimately the event Psi is calculated only on the junctions so the
    # flanking exons don't matter for this. But, all the exons are in
    # exon\d.bed in the index! And you, the lovely user, can decide what you
    # want to do with them!
    grouped = event_annotation.groupby(level=0, axis=0)

    n_events = len(grouped.size())

    if n_jobs == 1:
        # Do a separate branch because joblib doesn't do a good job of managing
        # the python debugger so use --n-jobs=1 (n_jobs=1) when debugging
        progress('\tIterating over {} events ...\n'.format(n_events))
        summaries = []
        for event_id, event_df in grouped:
            summary = _single_event_psi(
                event_id, event_df, junction_reads_2d,
                isoform1_junctions, isoform2_junctions,
                min_reads=min_reads,
                uneven_coverage_multiplier=uneven_coverage_multiplier,
                method=method)
            summaries.append(summary)
    else:
        processors = n_jobs if n_jobs > 0 else joblib.cpu_count()
        progress("\tParallelizing {} events' Psi calculation across {} "
                 "CPUs ...\n".format(n_events, processors))
        summaries = joblib.Parallel(n_jobs=n_jobs)(
            joblib.delayed(_single_event_psi)(
                event_id, event_df, junction_reads_2d,
                isoform1_junctions, isoform2_junctions,
                min_reads=min_reads,
                uneven_coverage_multiplier=uneven_coverage_multiplier,
                method=method)
            for event_id, event_df in grouped)

    return summaries


def calculate_psi(event_annotation, junction_reads_2d,
                  isoform1_junctions, isoform2_junctions,
                  min_reads=MIN_READS, method='mean',
                  uneven_coverage_multiplier=UNEVEN_COVERAGE_MULTIPLIER,
                  n_jobs=-1):
    """Compute percent-spliced-in of events based on junction reads

    Parameters
    ----------
    event_annotation : pandas.DataFrame
        A table where each row represents a single splicing event. The required
       columns are the ones specified in `isoform1_junctions`,
        `isoform2_junctions`, and `event_col`.
    junction_reads_2d : pandas.DataFrame
        A (n_samples, n_total_junctions) table of the number of reads found in
        all samples' exon-exon, all junctions. Very very large, e.g.
        1000 samples x 50,000 junctions = 50 million elements
        number of reads observed at a splice junction of a particular sample.
    isoform1_junctions : list
        Columns in `event_annotation` which represent junctions that
        correspond to isoform1, the Psi=0 isoform, e.g. ['junction13'] for SE
        (junctions between exon1 and exon3)
    isoform2_junctions : list
        Columns in `event_annotation` which represent junctions that
        correspond to isoform2, the Psi=1 isoform, e.g.
        ['junction12', 'junction23'] (junctions between exon1, exon2, and
        junction between exon2 and exon3)
    min_reads : int, optional
        Minimum number of reads for a junction to be viable. The rules
        governing compatibility of events are complex, and it is recommended to
        read the documentation for ``outrigger psi`` (default=10)
    method : "mean" | "min", optional
        Denotes the method by which to aggregate junctions from the same
        isoform - either use the mean (default) or the minimum.
        (default="mean")
    uneven_coverage_multiplier : int, optional
        Scale factor for the maximum amount bigger one side of a junction can
        be before rejecting the event, e.g. for an SE event with two junctions,
        junction12 and junction23, junction12=40 but junction23=500, then this
        event would be rejected because 500 > 40*10 (default=10)
    n_jobs : int, optional
        Number of subprocesses to create. Default is -1, which is to use as
        many processes/cores as possible

    Returns
    -------
    psi : pandas.DataFrame
        An (samples, events) dataframe of the percent spliced-in values
    summary : pandas.DataFrame
        A (n_samples * n_events, 7) shaped table with the sample id, junction
        reads, percent spliced-in (Psi), and notes on each event in each
        sample, that explains why or why not Psi was calculated
    """
    summaries = _maybe_parallelize_psi(event_annotation, junction_reads_2d,
                                  isoform1_junctions, isoform2_junctions,
                                  min_reads, method,
                                  uneven_coverage_multiplier, n_jobs)
    summary = pd.concat(summaries, ignore_index=True)

    psi = summary.pivot(index=SAMPLE_ID, columns=EVENT_ID, values=PSI)
    return psi, summary
