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
    return reads.groupby(level=1).sum()

def maybe_get_isoform_reads(splice_junction_reads, row, junctions, reads_col):
    if splice_junction_reads.index.levels[0].isin(row[junctions]).sum() > 0:
        return splice_junction_reads.loc[idx[row[junctions], :], reads_col]
    else:
        return pd.Series()



def calculate_psi(exons_to_junctions, splice_junction_reads,
                  reads_col='reads', min_reads=10,
                  isoform1_junctions=['junction12', 'junction23'],
                  isoform2_junctions=['junction13'], event_col='event_id',
                  debug=False):
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

        isoform1 = filter_and_sum(isoform1, min_reads, isoform2_junctions)
        isoform2 = filter_and_sum(isoform2, min_reads, isoform1_junctions)

        if isoform1.empty and isoform2.empty:
            # If both are empty after looking at this event
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
    psi_df = pd.concat(psis, axis=1)
    return psi_df
