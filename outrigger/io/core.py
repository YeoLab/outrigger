from ..common import EXON_START, EXON_STOP, JUNCTION_START, JUNCTION_STOP, \
    JUNCTION_ID, CHROM, STRAND


def add_exons_and_junction_ids(junction_reads):
    """Given junction locations, add exon locations and junction ids

    This assumes that the junction starts are one nucleotide after an exon
    (i.e. where the intron starts not where the exon ends) and similarly for
    the junction stops, that the junctions are one nucleotide before the exon.

    Parameters
    ----------
    junction_reads : pandas.DataFrame
        A tidy table of junction reads, with columns {chrom}, {start}, {stop}
        and {strand}.

    Returns
    -------
    reads_with_exons : pandas.DataFrame
        The same table, with the columns {exon_start}, {exon_stop} and
        {junction_id} added.
    """.format(chrom=CHROM, start=JUNCTION_START, stop=JUNCTION_STOP,
               strand=STRAND, exon_start=EXON_START, exon_stop=EXON_STOP,
               junction_id=JUNCTION_ID)

    # From STAR, exon_cols start one base pair down from the end of the intron
    junction_reads[EXON_START] = junction_reads[JUNCTION_STOP] + 1

    # From STAR, exon_cols stop one base pair up from the start of the intron
    junction_reads[EXON_STOP] = junction_reads[JUNCTION_START] - 1

    junction_reads[JUNCTION_ID] = 'junction:' + \
                                  junction_reads[CHROM].astype(str) + ':' + \
                                  junction_reads[JUNCTION_START].astype(str) \
                                  + '-' \
                                  + junction_reads[JUNCTION_STOP].astype(str) \
                                  + ':' \
                                  + junction_reads[STRAND].astype(str)
    return junction_reads
