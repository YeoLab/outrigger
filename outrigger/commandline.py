#!/usr/bin/env python

import argparse
import glob
import logging
import os
import sys
import warnings

import gffutils
import numpy as np

from outrigger.index import events, junctions, gtf
from outrigger.psi import compute
from outrigger.io import star
from outrigger import util

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import pandas as pd

SPLICE_JUNCTIONS_CSV = 'junctions/reads.csv'


class CommandLine(object):
    def __init__(self, input_options=None):
        self.parser = argparse.ArgumentParser(
            description='Calculate "percent-spliced in" (Psi) scores of '
                        'alternative splicing on a *de novo*, custom-built '
                        'splicing index')

        self.subparser = self.parser.add_subparsers(help='Sub-commands')

        # --- Subcommand to build the index of splicing events --- #
        index_parser = self.subparser.add_parser(
            'index', help='Build an index of splicing events using a graph '
                          'database on your junction reads and an annotation')

        index_parser.add_argument(
            '-o', '--output', required=False, type=str, action='store',
            default='./outrigger_output',
            help='Name of the folder where you saved the output from '
                 '"outrigger index" (default is ./outrigger_output, which is '
                 'relative to the directory where you called the program)". '
                 'You will need this file for the next step, "outrigger psi"'
                 ' (default="./outrigger_output")')

        junctions = index_parser.add_mutually_exclusive_group(required=True)
        junctions.add_argument(
            '-j', '--sj-out-tab', type=str, action='store',
            nargs='*', help='SJ.out.tab files from STAR aligner output')
        junctions.add_argument(
            '-c', '--junction-read-csv',
            help="Name of the splice junction files to calculate psi scores "
                 "on. If not provided, the compiled '{sj_csv}' file with all"
                 " the samples from the SJ.out.tab files that were used "
                 "during 'outrigger index' will be used. Not required if you "
                 "specify SJ.out.tab file with '--sj-out-tab'".format(
                        sj_csv=SPLICE_JUNCTIONS_CSV))
        index_parser.add_argument('-m', '--min-reads', type=int,
                                  action='store',
                                  required=False, default=10,
                                  help='Minimum number of reads per junction '
                                       'for that junction to count in creating'
                                       ' the index of splicing events '
                                       '(default=10)')

        gtf = index_parser.add_mutually_exclusive_group(required=True)
        gtf.add_argument('-g', '--gtf-filename', type=str, action='store',
                         help="Name of the gtf file you want to use. If "
                              "a gffutils feature database doesn't "
                              "already exist at this location plus '.db'"
                              " (e.g. if your gtf is gencode.v19."
                              "annotation.gtf, then the database is "
                              "inferred to be gencode.v19.annotation.gtf"
                              ".db), then a database will be auto-"
                              "created. Not required if you provide a "
                              "pre-built database with "
                              "'--gffutils-db'")
        gtf.add_argument('-d', '--gffutils-db', type=str,
                         action='store',
                         help="Name of the gffutils database file you "
                              "want to use. The exon IDs defined here "
                              "will be used in the function when creating"
                              " splicing event names. Not required if you"
                              " provide a gtf file with '--gtf-filename'")
        index_parser.add_argument('--debug', required=False, action='store_true',
                                help='If given, print debugging logging '
                                     'information to standard out')
        index_parser.set_defaults(func=self.index)

        # --- Subcommand to calculate psi on the built index --- #
        psi_parser = self.subparser.add_parser(
            'psi', help='Calculate "percent spliced-in" (Psi) values using '
                        'the splicing event index built with "outrigger '
                        'index"')
        psi_parser.add_argument('-i', '--index', required=False,
                                default='./outrigger_index',
                                help='Name of the folder where you saved the '
                                     'output from "outrigger index" (default '
                                     'is ./outrigger_index, which is relative '
                                     'to the directory where you called this '
                                     'program, assuming you have called '
                                     '"outrigger psi" in the same folder as '
                                     'you called "outrigger index")')
        splice_junctions = psi_parser.add_mutually_exclusive_group(
            required=False)
        splice_junctions.add_argument(
            '-c', '--junction-read-csv', required=False,
            help="Name of the splice junction files to calculate psi scores "
                 "on. If not provided, the compiled '{sj_csv}' file with all "
                 "the samples from the SJ.out.tab files that were used during "
                 "'outrigger index' will be used. Not required if you specify "
                 "SJ.out.tab file with '--sj-out-tab'".format(
                        sj_csv=SPLICE_JUNCTIONS_CSV))
        splice_junctions.add_argument(
            '-j', '--sj-out-tab', required=False,
            type=str, action='store', nargs='*',
            help='SJ.out.tab files from STAR aligner output. Not required if '
                 'you specify a file with "--splice-junction-csv"')
        psi_parser.add_argument('-m', '--min-reads', type=int, action='store',
                                required=False, default=10,
                                help='Minimum number of reads per junction for'
                                     ' calculating Psi (default=10)')
        psi_parser.add_argument('--reads-col', default='reads',
                                help="Name of column in --splice-junction-csv "
                                     "containing reads to use. "
                                     "(default='reads')")
        psi_parser.add_argument('--sample-id-col', default='sample_id',
                                help="Name of column in --splice-junction-csv"
                                     "containing sample ids to use. "
                                     "(default='sample_id')")
        psi_parser.add_argument('--junction-id-col',
                                default='junction_id',
                                help="Name of column in --splice-junction-csv"
                                     "containing the ID of the junction to use"
                                     ". Must match exactly with the junctions "
                                     "in the index."
                                     "(default='junction_id')")
        psi_parser.add_argument('--debug', required=False, action='store_true',
                                help='If given, print debugging logging '
                                     'information to standard out')
        psi_parser.set_defaults(func=self.psi)

        print(input_options)
        if input_options is None or len(input_options) == 0:
            self.parser.print_usage()
            self.args = None
        else:
            self.args = self.parser.parse_args(input_options)
        print(self.args)

        self.args.func()

    def index(self):
        index = Index(vars(self.args))
        index.execute()

    def psi(self):
        psi = Psi(vars(self.args))
        psi.execute()


class Subcommand(object):

    output = 'outrigger_output'

    def csv(self, sj_out_tab):
        """Create a csv file of compiled splice junctions"""

        util.progress('Reading SJ.out.files and creating a big splice junction'
                      ' matrix of reads spanning exon-exon junctions...')
        splice_junctions = star.read_multiple_sj_out_tab(sj_out_tab)
        splice_junctions['reads'] = splice_junctions['unique_junction_reads']

        filename = os.path.join(self.output, SPLICE_JUNCTIONS_CSV)
        util.progress('Writing {} ...\n'.format(filename))
        splice_junctions.to_csv(filename, index=False)
        util.done()
        return splice_junctions

    @staticmethod
    def maybe_make_db(gtf_filename, gffutils_db):
        """Get GFFutils database from file or create from a gtf"""
        if gffutils_db is not None:
            util.progress('Reading gffutils database from {} ...\n'.format(
                gffutils_db))
            db = gffutils.FeatureDB(gffutils_db)
            util.done()
        else:
            db_filename = '{}.db'.format(gtf_filename)
            try:
                db = gffutils.FeatureDB(db_filename)
            except ValueError:
                util.progress(
                    'Creating a "gffutils" '
                    'database {} ...\n'.format(db_filename))
                db = gtf.create_db(gtf_filename, db_filename)
                util.done()
        return db


class Index(Subcommand):

    def __init__(self, output=None, splice_junctions=None, min_reads=None,
                 gtf_filename=None, gffutils_db=None, debug=False):
        self.output = output
        self.splice_junctions = splice_junctions
        self.min_reads = min_reads
        self.gtf_filename = gtf_filename
        self.gffutils_db = gffutils_db
        self.debug = debug


    @staticmethod
    def junction_metadata(spliced_reads):
        """Get just the juunction info from the concatenated read files"""
        util.progress('Creating splice junction metadata of merely where '
                      'junctions start and stop')
        metadata = junctions.make_metadata(spliced_reads)
        util.done()
        return metadata


    @staticmethod
    def make_exon_junction_adjacencies(metadata, db):
        """Get annotated exons next to junctions in data"""
        util.progress('Getting junction-direction-exon triples for graph '
                      'database ...')
        exon_junction_adjacencies = junctions.ExonJunctionAdjacencies(
            metadata, db)
        junction_exon_triples = exon_junction_adjacencies.find_adjacencies()
        util.done()
        return junction_exon_triples

    def make_graph(self, junction_exon_triples, db):
        """Create graph database of exon-junction adjacencies"""
        util.progress('Populating graph database of the '
                      'junction-direction-exon triples ...')

        event_maker = events.EventMaker(junction_exon_triples, db)
        util.done()
        return event_maker

    def make_events_by_traversing_graph(self, event_maker):
        for name, abbrev in events.EVENT_TYPES:
            name_with_spaces = name.replace('_', ' ')
            # Find event junctions
            util.progress('Finding all {name} ({abbrev}) event ...'.format(
                name=name_with_spaces, abbrev=abbrev.upper()))
            events_of_type = getattr(event_maker, name)()
            util.done()

            # Write to a file
            csv = os.path.join(self.output, ['index', abbrev, 'junctions'])
            util.progress('Writing {abbrev} events to {csv} '
                          '...'.format(abbrev=abbrev.upper(), csv=csv))
            events_of_type.to_csv(csv, index=False)
            util.done()

    def execute(self):
        # Must output the junction exon triples
        logger = logging.getLogger('outrigger.index')

        if not os.path.exists(self.output):
            os.mkdir(self.output)

        if self.debug:
            logger.setLevel(10)

        spliced_reads = self.csv()
        metadata = self.junction_metadata(spliced_reads)

        db = self.maybe_make_db(self.gtf_filename, self.gffutils_db)

        junction_exon_triples = self.make_exon_junction_adjacencies(metadata, db)

        util.progress('Building exon-junction graph ...')
        event_maker = events.EventMaker(junction_exon_triples, db=db)
        util.done()
        self.make_events_by_traversing_graph(event_maker)



class Psi(Subcommand):

    def __init__(self, index=None, junction_read_csv=None, sj_out_tab=None,
                 min_reads=None, sample_id_col=None, reads_col=None,
                 junction_id_col=None, debug=None):
        self.index = index
        self.junction_read_csv = junction_read_csv
        self.sj_out_tab = sj_out_tab
        self.min_reads = min_reads
        self.sample_id_col = sample_id_col
        self.reads_col = reads_col
        self.junction_id_col = junction_id_col
        self.debug = debug

    def execute(self):
        """Calculate percent spliced in (psi) of splicing events"""

        logger = logging.getLogger('outrigger.psi')
        if self.debug:
            logger.setLevel(10)

        try:
            util.progress(
                'Reading splice junction reads from {} ...'.format(
                    self.junction_read_csv))
            dtype = {self.reads_col: np.float32}
            splice_junction_reads = pd.read_csv(
                self.junction_read_csv, dtype=dtype)

            try:
                assert self.reads_col in splice_junction_reads
            except AssertionError:
                raise (AssertionError(
                    'The column specifying reads, '
                    '"{}", is not contained in {}'.format(
                        self.reads_col), self.junction_read_csv))
            try:
                assert self.sample_id_col in splice_junction_reads
            except AssertionError:
                raise (AssertionError('The column specifying the sample '
                                      'ids, "{}", is not contained in '
                                      '{}'.format(
                    self.sample_id_col,
                    self.junction_read_csv)))
            try:
                assert self.junction_id_col in splice_junction_reads
            except AssertionError:
                raise (AssertionError(
                    'The column specifying the junction  '
                    'location, "{}", is not contained in '
                    '{}'.format(
                        self.junction_id_col,
                        self.junction_read_csv)))

            util.done()
        except KeyError:
            splice_junction_reads = self.csv(self.sj_out_tab, )

        splice_junction_reads = splice_junction_reads.set_index(
            [self.junction_id_col, self.sample_id_col])
        splice_junction_reads.sort_index(inplace=True)
        logger.debug('\n--- Splice Junction reads ---')
        logger.debug(repr(splice_junction_reads.head()))

        events_folder = os.path.join(self.index, 'events')
        psis = []
        for filename in glob.iglob('{}/*.csv'.format(events_folder)):
            event_type = os.path.basename(filename).split('.csv')[0]
            util.progress('Reading {} events from {} ...'.format(
                event_type.upper(), filename))

            isoform_junctions = compute.ISOFORM_JUNCTIONS[event_type]
            event_annotation = pd.read_csv(filename, index_col=0)
            util.done()
            logger.debug('\n--- Splicing event annotation ---')
            logger.debug(repr(event_annotation.head()))

            util.progress('Calculating percent spliced-in (Psi) '
                          'scores on {} events ...'.format(
                event_type.upper()))
            event_psi = compute.calculate_psi(
                event_annotation, splice_junction_reads,
                min_reads=self.min_reads, debug=self.debug,
                reads_col=self.reads_col, **isoform_junctions)
            event_psi.to_csv(os.path.join(self.index,
                                          '{}_psi.csv'.format(event_type)))
            psis.append(event_psi)
            util.done()

        util.progress('Concatenating all calculated psi scores '
                      'into one big matrix...')
        splicing = pd.concat(psis, axis=1)
        util.done()
        splicing = splicing.T
        csv = os.path.join(self.index, 'psi.csv')
        util.progress('Writing a samples x features matrix of Psi '
                      'scores to {} ...'.format(csv))
        splicing.to_csv(csv)
        util.done()



def main():
    CommandLine(sys.argv[1:])


if __name__ == '__main__':
    # try:
    main()
