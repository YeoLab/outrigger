__author__ = 'olgabotvinnik'

from collections import defaultdict

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

import pandas as pd

isoform1_seqs = []
isoform2_seqs = []

isoform1_translations = defaultdict(list)
isoform2_translations = defaultdict(list)

exon1_filename = '/projects/ps-yeolab/obotvinnik/miso_helpers/hg19/se_exon1.fasta'
exon2_filename = '/projects/ps-yeolab/obotvinnik/miso_helpers/hg19/se_exon2.fasta'
exon3_filename = '/projects/ps-yeolab/obotvinnik/miso_helpers/hg19/se_exon3.fasta'

def seq_name_to_exon_id(seq_name):
    """Convert a region location created from a BedTool.seq() call

    Parameters
    ----------
    seq_name : str
        A sequence name, as generated from a pybedtools.BedTool.seq() function

    Returns
    -------
    exon_id : str
        The base-1-ify'd version of this region into an exon id

    >>> seq_name_to_exon_id('chr1:100-200(+)')
    'exon:chr1:101-200:+
    """
    chr_start_stop, strand = seq_name.split('(')
    chrom, startstop = chr_start_stop.split(':')
    start, stop = startstop.split('-')

    # TODO: Double check this plus-oneing if it's necessary for base-0-ifying
    start = int(start)+1
    chr_start_stop = '{}:{}-{}'.format(chrom, start, stop)
    strand = strand.rstrip(')')
    exon = 'exon:{}:{}'.format(chr_start_stop, strand)
    return exon

def splice_type_exons(splice_type, exons):
    """Get exons corresponding to a particular isoform of a splice type

    For SE:
        isoform1: exon1, exon2
        isoform2: exon1, exon2, exon3
    for MXE:
        isoform1: exon1, exon3, exon4
        isoform2: exon1, exon2, exon4

    Parameters
    ----------
    splice_type : 'SE' | 'MXE'
        String specifying the splice type. Currently only SE (skipped exon) and
        MXE (mutually exclusive exon) are supported
    exons : list
        List of exons or CDS's (ids, strings, you name it) in the exact order
        of the splice type, e.g. (exon1_id, exon2_id, exon3_id) for SE

    Returns
    -------
    isoform1_exons : tuple
        Tuple of exons corresponding to isoform 1
    isoform2_exons : tuple
        Tuple of exons corresponding to isoform 2
    """
    isoform1, isoform2 = None, None
    if splice_type == 'SE':
        isoform1 = exons[0], exons[2]
        isoform2 = exons[0], exons[1], exons[2]
    elif splice_type == 'MXE':
        isoform1 = exons[0], exons[1], exons[3]
        isoform2 = exons[0], exons[2], exons[3]
    return isoform1, isoform2


def exon_seqs_to_isoform_seqs(exon_fastas, event_ids, splice_type):
    parsed = [SeqIO.parse(exon_fasta, 'fasta') for exon_fasta in exon_fastas]

    isoform_seqs = defaultdict(list)
    for event_id, seqs in zip(event_ids, parsed):
        isoforms = splice_type_exons(splice_type, seqs)
        for i, isoform in enumerate(isoforms):
            i = i+1
            seq = ''.join(x.seq for x in isoform)
            seqrecord = SeqRecord(seq, id=event_id, description='isoform{0}'.format(i))
            isoform_seqs[i].append(seqrecord)
    return isoform_seqs


with open(exon1_filename) as infile1, open(exon2_filename) as infile2, open(exon3_filename) as infile3:
    parsed1 = SeqIO.parse(infile1, 'fasta')
    parsed2 = SeqIO.parse(infile2, 'fasta')
    parsed3 = SeqIO.parse(infile3, 'fasta')
    for i, (seq1, seq2, seq3) in enumerate(zip(parsed1, parsed2, parsed3)):
        event_name = attribute_df.index[i]
#         print seq1.name
        exon_id1 = seq_name_to_exon_id(seq1.name)
        exon_id2 = seq_name_to_exon_id(seq2.name)
        exon_id3 = seq_name_to_exon_id(seq3.name)

        cds_id1 = ':'.join(['CDS', exon_id1.split('exon:')[1]])
        cds_id2 = ':'.join(['CDS', exon_id2.split('exon:')[1]])
        cds_id3 = ':'.join(['CDS', exon_id3.split('exon:')[1]])

        try:
#             print seq1.name, exon_id1
            exon1 = v19db[exon_id1]
            exon2 = v19db[exon_id2]
            exon3 = v19db[exon_id3]
        except gffutils.FeatureNotFoundError:
            continue
#         print exon1
#         print

        transcripts1 = v19db.parents(exon1, featuretype='transcript')
        transcripts2 = v19db.parents(exon2, featuretype='transcript')
        transcripts3 = v19db.parents(exon3, featuretype='transcript')

        isoform1s = set(transcripts1) & set(transcripts3)
        isoform2s = set(transcripts2) & isoform1s

        isoform1s = isoform1s.difference(isoform2s)
#         print 'isoform1', isoform1s
#         print 'isoform2', isoform2s

        for isoform in isoform1s:
            event_isoform = '{}_isoform1'.format(event_name)
            name = '{}_{}'.format(event_isoform, isoform.id)
            reverse = exon_id1[-1] == '-'
            cds = v19db.children(isoform, featuretype='CDS', reverse=reverse, order_by='start')

            cds_in_splice_form = [c for c in cds if c.id.startswith(cds_id1) or c.id.startswith(cds_id3)]
#             print 'cds_in_splice_form', cds_in_splice_form

            if len(cds_in_splice_form) == 2:
                frame1 = cds_in_splice_form[0].frame
                frame3 = cds_in_splice_form[1].frame

                seq = seq1.seq + seq2.seq

                seq_translated = seq[int(frame1):].translate()
                if seq_translated in isoform1_translations[event_isoform]:
                    continue

#                 print 'name', name
                result_seq = SeqRecord(seq_translated, id=name, description='')
                isoform1_seqs.append(result_seq)
                isoform1_translations[event_isoform].append(seq_translated)
            else:
                isoform1_translations[event_isoform].append('no translation')

        for isoform in isoform2s:
            event_isoform = '{}_isoform2'.format(event_name)
            name = '{}_{}'.format(event_isoform, isoform.id)
            reverse = exon_id1[-1] == '-'
            cds = v19db.children(isoform, featuretype='CDS', reverse=reverse, order_by='start')

            cds_in_splice_form = [c for c in cds if c.id.startswith(cds_id1) or c.id.startswith(cds_id2) or c.id.startswith(cds_id3)]
#             print 'cds_in_splice_form', cds_in_splice_form

            if len(cds_in_splice_form) == 3:
                frame1 = cds_in_splice_form[0].frame
                frame2 = cds_in_splice_form[1].frame
                frame3 = cds_in_splice_form[2].frame

                seq = seq1.seq + seq2.seq + seq3.seq
                seq_translated = seq[int(frame1):].translate()
                if seq_translated in isoform2_translations[event_isoform]:
                    continue

#                 print 'name', name
                result_seq = SeqRecord(seq_translated, id=name, description='')
                isoform2_seqs.append(result_seq)
                isoform2_translations[event_isoform].append(seq_translated)
#                 isoform1_names.add(name)
            else:
                isoform2_translations[event_isoform].append('no translation')


#         for f in frame.split(','):
#             f = int(f)
# #                 print f,
#             seq_translated = seq.seq[f:].translate()
#             name = '{}_frame{}'.format(ind, f)
#             result_seq = SeqRecord(seq_translated, id=name, description='')
#             translated.ix[ind, f] = str(result_seq.seq)
#             result_seqs.append(result_seq)

filename = '/projects/ps-yeolab/obotvinnik/miso_helpers/hg19/se_isoform1_translated.fa'
with open(filename, 'w') as outfile:
    SeqIO.write(isoform1_seqs, outfile, 'fasta')

filename = '/projects/ps-yeolab/obotvinnik/miso_helpers/hg19/se_isoform2_translated.fa'
with open(filename, 'w') as outfile:
    SeqIO.write(isoform2_seqs, outfile, 'fasta')





events_to_domain_disruptions = pd.Series()

# class Isoform(object):

no_translation = 'no translation'
translation_no_domain_match = 'translation but no domain match'

no_domains = (no_translation, translation_no_domain_match)

def event_to_domain_disruption(event_name, row):
    isoform1, isoform2 = None, None
#     print row

    if pd.isnull(row.isoform1_translation):
        isoform1 = 'no translation'
    if pd.isnull(row.isoform2_translation):
        isoform2 = 'no translation'


    isoform1_domains = isoform1_pfam_domain_name.ix[event_name].dropna()
    isoform2_domains = isoform2_pfam_domain_name.ix[event_name].dropna()
    intersection = len(isoform1_domains.index.intersection(isoform2_domains.index))
    union = len(isoform1_domains.index.union(isoform2_domains.index))

#     if isoform1_domains.count() == 0 and isoform2_domains.count() == 0:
#         return 'translation but no domain match'
    if isoform1_domains.count() == 0 and isoform1 != no_translation:
        isoform1 = translation_no_domain_match
    if isoform2_domains.count() == 0 and isoform2 != no_translation:
        isoform2 = translation_no_domain_match

    if isoform1 in no_domains:
        if isoform2_domains.count() > 0:
            isoform2 = 'domain'
            return '{} --> {}'.format(isoform1, isoform2)
        else:
            return '{} --> {}'.format(isoform1, isoform2)
    elif isoform2 in no_domains:
        if isoform1_domains.count() > 0:
            isoform1 = 'domain'
            return '{} --> {}'.format(isoform1, isoform2)
        else:
            return '{} --> {}'.format(isoform1, isoform2)

#     if event_name == 'chr10:103345619:103345913:-@chr10:103344359:103344676:-@chr10:103343265:103343438:-':
#         import pdb; pdb.set_trace()

    if isoform1 != 'no translation' and isoform2 != 'no translation':
        if isoform1_domains.equals(isoform2_domains):
            isoform1, isoform2 = 'same domains', 'same domains'
            return 'same domains'
        elif union == intersection:
            diff = isoform1_domains-isoform2_domains
            if (diff < 0).any():
                return 'same domains, increased count'
            elif (diff > 0).any():
                isoform1 = 'same domains, more count'
                isoform2 = 'domain'
                return 'same domains, decreased count'
        elif union > intersection:
            if intersection == 0:
                return 'non-overlapping domains'
            else:
                if len(isoform1_domains.index.difference(isoform2_domains.index)) > 0:
                    if len(isoform2_domains.index.difference(isoform1_domains.index)) == 0:
                        isoform1 = 'more domains'
                        isoform2 = 'fewer domains'
                        return 'overlapping domains, loss of unique domains'
                    else:
                        isoform1 = 'overlapping domains, unique domains'
                        isoform2 = 'overlapping domains, unique domains'
                        return 'overlapping domains, each isoform has unique'
                if len(isoform2_domains.index.difference(isoform1_domains.index)) > 0:
                    return 'overlapping domains, gain of unique domains'

# translations = study.splicing.feature_data.ix[:, ['isoform1_translation', 'isoform2_translation']].dropna(how='all')
# translations.head()
#
# for event_name, row in translations.iterrows():
#     events_to_domain_disruptions[event_name] = event_to_domain_disruption(event_name, row)
#
# events_to_domain_disruptions