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


def maybe_analyze(db):
    try:
        # For gffutils >0.8.7.1
        db.analyze()
    except AttributeError:
        # For compatability with gffutils<=0.8.7.1
        db.execute('ANALYZE features')


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

    db = gffutils.create_db(
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
    maybe_analyze(db)
    return db


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
        self.regions = pd.DataFrame(index=self.events.index)
        self.region_cols = ['{}_region'.format(x) for x in self.exon_cols]

        for exon_col, region_col in zip(self.exon_cols, self.region_cols):
            self.regions[region_col] = self.events[exon_col].map(Region)

        # Make introns and copy-pastable genome locations for the whole event
        intron_regions = self.regions[self.region_cols].apply(
            self.event_introns_regions, axis=1)

        self.regions = pd.concat([self.regions, intron_regions], axis=1)
        self.region_cols.extend(['intron_region', 'event_region'])

        # Add the lengths of exons, introns, event region, and the genome
        # location ("name") of each intron
        self.lengths = self.regions.applymap(len)
        intron_names = intron_regions.applymap(lambda x: x.name)
        self.events = pd.concat([self.events, self.lengths, intron_names],
                                axis=1)

    def attributes(self):
        """Retrieve all GTF attributes for each isoform's event"""

        lines = []

        for event_id, row in self.events.iterrows():
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
                attributes = pd.Series(attributes, name=event_id)
                attributes.index = isoform + '_' + attributes.index
                lines.append(attributes)
        event_attributes = pd.concat(lines, axis=1).T
        event_attributes = self.events.join(event_attributes)
        return event_attributes

    def exon_bedfiles(self, folder):
        for region_col in self.region_cols:
            column = self.regions[region_col]
            lines = (region.to_bed_format(event_id)
                     for event_id, region in column.iteritems())

            name = region_col.split('_')[0]
            basename = name + '.bed'
            filename = os.path.join(folder, basename)

            with open(filename, 'w') as f:
                f.write('\n'.join(lines) + '\n')

    def event_introns_regions(self, exons):
        """Make intron and event regions for an event

        Parameters
        ----------
        exons : outrigger.Regions
            List of exon ids, e.g. ["exon:chr1:100-200:+",
            "exon:chr1:300-400:+"]

        Returns
        -------
        regions : dict

        """
        first_exon = exons[0]
        last_exon = exons[-1]

        chrom = first_exon.chrom
        strand = first_exon.strand

        if strand == '-':
            intron_stop = first_exon.start
            intron_start = last_exon.stop

            event_start = last_exon.start
            event_stop = first_exon.stop
        else:
            # If strand is positive or undefined
            intron_start = first_exon.stop
            intron_stop = last_exon.start

            event_start = first_exon.start
            event_stop = last_exon.stop

        intron = Region('intron:{chrom}:{start}-{stop}:{strand}'.format(
            chrom=chrom, start=intron_start, stop=intron_stop,
            strand=strand))
        event = Region('event:{chrom}:{start}-{stop}:{strand}'.format(
            chrom=chrom, start=event_start, stop=event_stop, strand=strand))

        regions = pd.Series(dict(intron_region=intron, event_region=event))

        return regions
