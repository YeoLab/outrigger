from collections import Counter

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
    columns = ['target_name', 'target_accession', 'target_length', 'query_name', 'query_accession', 'query_length', 'sequence_e_value',
             'sequence_score', 'sequence_bias', 'domain_number', 'domain_total', 'domain_conditional_e_value', 'domain_independent_e_value',
             'domain_score', 'domain_bias', 'target_start', 'target_stop', 'query_start', 'query_stop', 'query_domain_envelope_start',
             'query_domain_envelope_stop', 'mean_posterior_probability', 'target_description']
    df = pd.DataFrame.from_records(entries, columns=columns)
    df = df.convert_objects(convert_numeric=True)
    return df



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

translations = study.splicing.feature_data.ix[:, ['isoform1_translation', 'isoform2_translation']].dropna(how='all')
translations.head()

for event_name, row in translations.iterrows():
    events_to_domain_disruptions[event_name] = event_to_domain_disruption(event_name, row)

events_to_domain_disruptions