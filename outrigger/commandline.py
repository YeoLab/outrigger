#!/usr/bin/env python

import argparse
import glob
import logging
import os
import sys
import warnings

import gffutils
import numpy as np

from outrigger import events, gtf, junctions, psi, star, util


with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import pandas as pd

SPLICE_JUNCTIONS_CSV = 'splice_junctions.csv'


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
            '-i', '--index', required=False, type=str, action='store',
            default='./outrigger_index',
            help='Name of the folder where you saved the output from '
                 '"outrigger index" (default is ./outrigger_index, which is '
                 'relative to the directory where you called the program)". '
                 'You will need this file for the next step, "outrigger psi"')

        junctions = index_parser.add_mutually_exclusive_group(required=True)
        junctions.add_argument(
            '-j', '--sj-out-tab', type=str, action='store',
            nargs='*', help='SJ.out.tab files from STAR aligner output')
        junctions.add_argument(
            '-c', '--splice-junction-csv',
            help="Name of the splice junction files to calculate psi scores "
                 "on. If not provided, the compiled '{sj_csv}' file with all"
                 " the samples from the SJ.out.tab files thatm were used "
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
        gtf.add_argument('-g', '--gtf', type=str, action='store',
                         help="Name of the gtf file you want to use. If "
                              "a gffutils feature database doesn't "
                              "already exist at this location plus '.db'"
                              " (e.g. if your gtf is gencode.v19."
                              "annotation.gtf, then the database is "
                              "inferred to be gencode.v19.annotation.gtf"
                              ".db), then a database will be auto-"
                              "created. Not required if you provide a "
                              "pre-built database with "
                              "'--gffutils-database'")
        gtf.add_argument('-d', '--gffutils-database', type=str,
                         action='store',
                         help="Name of the gffutils database file you "
                              "want to use. The exon IDs defined here "
                              "will be used in the function when creating"
                              " splicing event names. Not required if you"
                              " provide a gtf file with '--gtf'")
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
            '-c', '--splice-junction-csv', required=False,
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
        psi_parser.add_argument('--junction-location-col',
                                default='junction_location',
                                help="Name of column in --splice-junction-csv"
                                     "containing the ID of the junction to use"
                                     ". Must match exactly with the junctions "
                                     "in the index."
                                     "(default='junction_location')")
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

    def csv(self):
        """Create a csv file of compiled splice junctions"""

        util.progress('Reading SJ.out.files and creating a big splice junction'
                      ' matrix of reads spanning exon-exon junctions...')
        splice_junctions = star.read_multiple_sj_out_tab(self.args.sj_out_tab)
        splice_junctions['reads'] = splice_junctions['unique_junction_reads']

        filename = os.path.join(self.args.index, SPLICE_JUNCTIONS_CSV)
        util.progress('{Writing {} ...\n'.format(filename))
        splice_junctions.to_csv(filename, index=False)
        util.done()
        return splice_junctions

    def junction_metadata(self):
        """Get just the juunction info from the concatenated read files"""
        util.progress('Creating splice junction metadata of merely where '
                      'junctions start and stop')

        util.done()

    def make_db(self):
        if self.args.gffutils_database is not None:
            util.progress('Reading gffutils database from {} ...\n'.format(
                self.args.gffutils_database))
            db = gffutils.FeatureDB(self.args.gffutils_database)
            util.done()
        else:
            db_filename = '{}.db'.format(self.args.gtf)
            try:
                db = gffutils.FeatureDB(db_filename)
            except ValueError:
                util.progress(
                    'Creating a "gffutils" '
                    'database {} ...\n'.format(db_filename))
                db = gtf.create_db(self.args.gtf, db_filename)
                util.done()
        return db

    def make_exon_junction_adjacencies(self):
        """Get annotated exons next to junctions in data"""
        eja = junctions.ExonJunctionAdjacencies(self.splice_junctions, self.db)
        pass

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
            csv = os.path.join(self.args.index, ['index', abbrev, 'junctions'])
            util.progress('Writing {abbrev} events to {csv} '
                          '...'.format(abbrev=abbrev.upper(), csv=csv))
            events_of_type.to_csv(csv, index=False)
            util.done()

    def index(self):
        # Must output the junction exon triples
        logger = logging.getLogger('outrigger.index')

        if not os.path.exists(self.args.index):
            os.mkdir(self.args.index)

        if self.args.debug:
            logger.setLevel(10)

        spliced_reads = self.csv()


        db = self.make_db()

        util.progress('Getting junction-direction-exon triples for graph '
                      'database ...')
        exon_junction_adjacencies = junctions.ExonJunctionAdjacencies(
            spliced_reads, db)
        junction_exon_triples = exon_junction_adjacencies.get_adjacent_exons()

        util.progress('Building exon-junction graph ...')
        event_maker = events.EventMaker(junction_exon_triples, db=db)
        util.done()
        self.make_events_by_traversing_graph(event_maker)


    def psi(self):
        """Calculate percent spliced in (psi) of splicing events"""

        logger = logging.getLogger('outrigger.psi')
        if self.args.debug:
            logger.setLevel(10)

        try:
            util.progress(
                'Reading splice junction reads from {} ...'.format(
                    self.args.splice_junction_csv))
            dtype = {self.args.reads_col: np.float32}
            splice_junction_reads = pd.read_csv(
                self.args.splice_junction_csv, dtype=dtype)

            try:
                assert self.args.reads_col in splice_junction_reads
            except AssertionError:
                raise (AssertionError(
                    'The column specifying reads, '
                    '"{}", is not contained in {}'.format(
                        self.args.reads_col), self.args.splice_junction_csv))
            try:
                assert self.args.sample_id_col in splice_junction_reads
            except AssertionError:
                raise (AssertionError('The column specifying the sample '
                                      'ids, "{}", is not contained in '
                                      '{}'.format(
                                        self.args.sample_id_col,
                                        self.args.splice_junction_csv)))
            try:
                assert self.args.junction_location_col in splice_junction_reads
            except AssertionError:
                raise (AssertionError(
                    'The column specifying the junction  '
                    'location, "{}", is not contained in '
                    '{}'.format(
                        self.args.junction_location_col,
                        self.args.splice_junction_csv)))

            util.done()
        except KeyError:
            util.progress('Creating consolidated splice junction '
                             'file "{sj_csv}" from SJ.out.tab files '
                          '...'.format(sj_csv=SPLICE_JUNCTIONS_CSV))
            splice_junction_reads = self.csv()
            util.done()

        splice_junction_reads = splice_junction_reads.set_index(
            [self.args.junction_location_col, self.args.sample_id_col])
        splice_junction_reads.sort_index(inplace=True)
        logger.debug('\n--- Splice Junction reads ---')
        logger.debug(repr(splice_junction_reads.head()))

        events_folder = os.path.join(self.args.index, 'events')
        psis = []
        for filename in glob.iglob('{}/*.csv'.format(events_folder)):
            event_type = os.path.basename(filename).split('.csv')[0]
            util.progress('Reading {} events from {} ...'.format(
                event_type.upper(), filename))

            isoform_junctions = psi.ISOFORM_JUNCTIONS[event_type]
            event_annotation = pd.read_csv(filename, index_col=0)
            util.done()
            logger.debug('\n--- Splicing event annotation ---')
            logger.debug(repr(event_annotation.head()))

            util.progress('Calculating percent spliced-in (Psi) '
                          'scores on {} events ...'.format(
                event_type.upper()))
            event_psi = psi.calculate_psi(
                event_annotation, splice_junction_reads,
                min_reads=self.args.min_reads, debug=self.args.debug,
                reads_col=self.args.reads_col, **isoform_junctions)
            event_psi.to_csv(os.path.join(self.args.index,
                                          '{}_psi.csv'.format(event_type)))
            psis.append(event_psi)
            util.done()

        util.progress('Concatenating all calculated psi scores '
                         'into one big matrix...')
        splicing = pd.concat(psis, axis=1)
        util.done()
        splicing = splicing.T
        csv = os.path.join(self.args.index, 'psi.csv')
        util.progress('Writing a samples x features matrix of Psi '
                      'scores to {} ...'.format(csv))
        splicing.to_csv(csv)
        util.done()


def main():
    CommandLine(sys.argv[1:])


if __name__ == '__main__':
    # try:
    main()
