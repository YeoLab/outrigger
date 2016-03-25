
import glob
import os

import itertools
import pandas as pd

from .region import Region
from ..io.common import STRAND

SPLICE_TYPE_ISOFORM_EXONS = {'SE': {'isoform1': ['exon1', 'exon3'],
                                    'isoform2': ['exon1', 'exon2', 'exon3']},
                             'MXE': {'isoform1': ['exon1', 'exon3', 'exon4'],
                                     'isoform2': ['exon1', 'exon2', 'exon4']}
                             }


class SplicingAnnotator(object):
    """Annotates basic features of splicing events: gene ids and names"""

    def __init__(self, db, events, splice_type):
        self.db = db
        self.events = events
        self.splice_type = splice_type
        self.isoform_exons = SPLICE_TYPE_ISOFORM_EXONS[self.splice_type]
        self.exons = list(set(itertools.chain(*self.isoform_exons.values())))
        self.exons.sort()

    def attributes(self):
        """Retrieve all GTF attributes for each isoform's event"""

        lines = []

        for row_index, row in self.events.iterrows():
            for isoform, exons in self.isoform_exons.items():
                attributes = {}

                exon1 = self.db[row[exons[0]]]
                other_exons = row[exons[1:]]
                for (key, value) in exon1.attributes.items():
                    # v2 = set.intersection(*[set(self.db[e].attributes[key])
                    #                            for e in other_exons])
                    other_values = set([])
                    for i, e in enumerate(other_exons):
                        try:
                            if i == 0:
                                other_values = set(self.db[e].attributes[key])
                            else:
                                other_values.intersection_update(
                                    self.db[e].attributes[key])
                        except KeyError:
                            # i -= 1
                            continue

                    intersection = set(value) & other_values
                    if len(intersection) > 0:
                        attributes[key] = ','.join(sorted(list(intersection)))
                attributes = pd.Series(attributes, name=row['exons'])
                attributes.index = isoform + '_' + attributes.index
                lines.append(attributes)
        # import pdb; pdb.set_trace()
        return pd.concat(lines, axis=1).T

    def lengths(self):
        """Retrieve exon and intron lengths for an event"""
        df = pd.DataFrame(index=self.events.index)

        for exon_col in self.exons:
            df['{}_region'.format(exon_col)] = self.events[exon_col].map(Region)
            df['{}_length'.format(exon_col)] = df['{}_region'.format(
                exon_col)].map(len)

        first_exon = '{}_region'.format(self.exons[0])
        last_exon = '{}_region'.format(self.exons[-1])
        positive = self.events[STRAND] == '+'
        df.loc[positive, 'intron_length'] = \
            df.loc[positive].apply(
                lambda x: x[last_exon].start - x[first_exon].stop - 1, axis=1)
        df.loc[~positive, 'intron_length'] = \
            df.loc[~positive].apply(
                lambda x: x[first_exon].start - x[last_exon].stop - 1, axis=1)

        regions = [x for x in df if x.endswith('_region')]
        df = df.drop(regions, axis=1)
        df.index = self.events['exons']
        return df


