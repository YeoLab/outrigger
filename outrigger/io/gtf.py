"""
Functions for creating GTF databases using gffutils and using those databases
to annotate alternative events.
"""
from collections import Counter
import itertools
import os

import gffutils
from gffutils.helpers import merge_attributes
import pandas as pd

from ..common import SPLICE_TYPE_ISOFORM_EXONS, OUTRIGGER_DE_NOVO, NOVEL_EXON
from ..region import Region, STRANDS

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
        """Annotate splicing events with their respective genes

        Parameters
        ----------
        db : gffutils.FeatureDB
            Database including all the exons found in the events
        events : pandas.DataFrame
            Table of events, with the event ids as the index
        splice_type : 'se' | 'mxe'
            The type of alternative splicing, which informs the exon
            configurations for different isoforms
        """
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
        self.lengths.columns = [x.replace('_region', '_length')
                                for x in self.lengths]
        self.lengths = self.lengths.astype(int)

        intron_names = intron_regions.applymap(lambda x: x.name)
        intron_names.columns = [x.replace('_region', '_location')
                                for x in intron_names]
        self.events = pd.concat([self.events, self.lengths, intron_names],
                                axis=1)

    def maybe_get_feature(self, feature_id):
        try:
            return self.db[feature_id]
        except gffutils.FeatureNotFoundError:
            r = Region(feature_id)
            feature = location_to_feature(self.db, r.chrom, r.start, r.stop,
                                          r.strand, source=OUTRIGGER_DE_NOVO,
                                          featuretype=NOVEL_EXON)
            self.db.update([feature], make_backup=False,
                           id_spec={NOVEL_EXON: 'location_id'},
                           transform=transform)
            return feature

    def attributes(self):
        """Retrieve all GTF attributes for each isoform's event"""

        ignore_keys = 'location_id', 'exon_id', 'exon_number'

        lines = []

        for event_id, row in self.events.iterrows():
            attributes = pd.Series(name=event_id)
            for isoform, exons in self.isoform_exons.items():

                for e in exons:
                    attributes[e] = row[e]

                n_exons = len(exons)

                exon_ids = row[exons]
                exon_features = [self.maybe_get_feature(exon_id) for
                                 exon_id in exon_ids]

                keys = set(itertools.chain(
                    *[exon.attributes.keys() for exon in exon_features]))

                for key in keys:
                    # Skip the location IDs which is specific to the
                    # outrigger-built database, and the exon ids which will
                    # never match up across all exons
                    if key in ignore_keys:
                        continue
                    values = Counter()

                    for exon_id in exon_ids:
                        try:
                            values.update(
                                self.db[exon_id].attributes[key])
                        except KeyError:
                            continue
                    if len(values) > 0:
                        # Only use attributes that came up in for all exons
                        # of the isoform
                        values = [value for value, count in values.items()
                                  if count == n_exons]
                        new_key = isoform + '_' + key
                        attributes[new_key] = ','.join(sorted(values))
            lines.append(attributes)

        event_attributes = pd.concat(lines, axis=1).T
        df = pd.concat([self.events, event_attributes], axis=1)
        df = df.loc[:, ~df.columns.duplicated()]
        return df

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


def location_to_feature(db, chrom, start, stop, strand, source, featuretype):
    if strand not in STRANDS:
        strand = '.'
    overlapping_genes = db.region(seqid=chrom, start=start, end=stop,
                                  strand=strand, featuretype='gene')

    exon_id = 'exon:{chrom}:{start}-{stop}:{strand}'.format(
        chrom=chrom, start=start, stop=stop, strand=strand)

    attributes = {}
    for g in overlapping_genes:
        attributes = merge_attributes(attributes, g.attributes)

    exon = gffutils.Feature(chrom, source=source,
                            featuretype=featuretype, start=start,
                            end=stop, strand=strand, id=exon_id,
                            attributes=attributes)
    return exon
