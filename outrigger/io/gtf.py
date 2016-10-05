"""
Functions for creating GTF databases using gffutils and using those databases
to annotate alternative events.
"""

import itertools
import os

import gffutils
import pandas as pd

from ..common import SPLICE_TYPE_ISOFORM_EXONS
from ..region import Region

# Annotations from:
# ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz

gene_transcript = set(('gene', 'transcript'))


def transform(f):
    if f.featuretype in gene_transcript:
        return f
    else:
        exon_location = '{}:{}:{}-{}:{}'.format(
            f.featuretype, f.seqid, f.start, f.stop, f.strand)
        exon_id = exon_location
        if f.featuretype == 'CDS':
            exon_id += ':' + f.frame
        f.attributes['location_id'] = [exon_id]
        return f


def create_db(gtf_filename, db_filename=None):
    db_filename = ':memory:' if db_filename is None else db_filename

    return gffutils.create_db(
        gtf_filename,
        db_filename,
        merge_strategy='merge',
        id_spec={'gene': 'gene_id', 'transcript': 'transcript_id',
                 'exon': 'location_id', 'CDS': 'location_id',
                 'start_codon': 'location_id',
                 'stop_codon': 'location_id', 'UTR': 'location_id'},
        transform=transform,
        force=True,
        verbose=True,
        disable_infer_genes=True,
        disable_infer_transcripts=True,
        force_merge_fields=['source'])


class SplicingAnnotator(object):
    """Annotates basic features of splicing events: gene ids and names"""

    def __init__(self, db, events, splice_type):
        self.db = db
        self.events = events
        self.splice_type = splice_type
        self.isoform_exons = SPLICE_TYPE_ISOFORM_EXONS[
            self.splice_type.lower()]
        self.exon_cols = list(set(itertools.chain(
            *self.isoform_exons.values())))
        self.exon_cols.sort()

        # Make a dataframe with outrigger.Region objects
        self.exon_regions = pd.DataFrame(index=self.events.index)
        self.region_cols = ['{}_region'.format(x) for x in self.exon_cols]
        for exon_col, region_col in zip(self.exon_cols, self.region_cols):
            self.exon_regions[region_col] = self.events[exon_col].map(Region)

    def attributes(self):
        """Retrieve all GTF attributes for each isoform's event"""

        lines = []

        for event_id, row in self.events.iterrows():
            for isoform, exons in self.isoform_exons.items():
                attributes = {}

                exon1 = self.db[row[exons[0]]]
                other_exons = row[exons[1:]]

                lengths = self.lengths(row[exons])
                attributes.update(lengths)
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
                attributes = pd.Series(attributes, name=event_id)
                attributes.index = isoform + '_' + attributes.index
                lines.append(attributes)
        event_attributes = pd.concat(lines, axis=1).T
        event_attributes = self.events.join(event_attributes)
        return event_attributes

    def exon_bedfiles(self, folder):
        for region_col in self.region_cols:
            column = self.exon_regions[region_col]
            lines = (region.to_bed_format(event_id)
                     for event_id, region in column.iteritems())

            exon_name = region_col.split('_')[0]
            basename = exon_name + '.bed'
            filename = os.path.join(folder, basename)

            with open(filename, 'w') as f:
                f.write('\n'.join(lines) + '\n')

    def lengths(self, exons):
        """Retrieve exon and intron lengths for an event

        Parameters
        ----------
        exons : list
            List of exon ids, e.g. ["exon:chr1:100-200:+",
            "exon:chr1:300-400:+"]

        Returns
        -------
        lengths : dict

        """
        regions = tuple(map(Region, exons))

        lengths = dict(('exon{}'.format(i+1), len(exon))
                       for i, exon in enumerate(regions))
        strand = regions[0].strand

        first_exon = regions[0]
        last_exon = regions[1]
        if strand == '-':
            intron_start = first_exon.start
            intron_stop = last_exon.stop
            intron_length = first_exon.start - last_exon.stop - 1

        else:
            # If strand is positive or undefined
            intron_length = last_exon.start - first_exon.stop - 1

        intron = 'intron:{chrom}:{start}:{stop}:{strand}'.format(
            chrom=first_exon.chrom, start=intron_start, stop=intron_stop,
            strand=strand)

        lengths.update(intron_length=int(intron_length), intron=intron)

        return lengths
