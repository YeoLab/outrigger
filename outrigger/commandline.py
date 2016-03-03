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

    @staticmethod
    def _done():
        """Write timestamp plus 'Done.' to stdout"""
        sys.stdout.write('{}\t\tDone.\n'.format(util.timestamp()))

    def csv(self):
        """Create a csv file of compiled splice junctions"""
        splice_junctions = star.read_multiple_sj_out_tab(self.args.sj_out_tab)
        splice_junctions['reads'] = splice_junctions['uniquely_mapped_reads']

        filename = os.path.join(self.args.index, SPLICE_JUNCTIONS_CSV)
        sys.stdout.write('{}\tWriting {} ...\n'.format(util.timestamp(),
                                                       filename))
        splice_junctions.to_csv(filename, index=False)
        return splice_junctions

    def index(self):
        # Must output the junction exon triples
        logger = logging.getLogger('outrigger.index')

        if self.args.deug:
            logger.setLevel(10)
        sys.stdout.write('{}\tReading SJ.out.files and creating a big splice '
                         'junction matrix ...\n'.format(util.timestamp()))
        splice_junctions = self.csv()
        self._done()

        if self.args.gffutils_database is not None:
            sys.stdout.write(
                '{}\tReading gffutils database from {} ...\n'.format(
                                util.timestamp(), self.args.gffutils_database))
            db = gffutils.FeatureDB(self.args.gffutils_database)
            self._done()
        else:
            db_filename = '{}.db'.format(self.args.gtf)
            try:
                db = gffutils.FeatureDB(db_filename)
            except ValueError:
                sys.stdout.write(
                    '{}\tCreating a "gffutils" '
                    'database {} ...\n'.format(util.timestamp(), db_filename))
                db = gtf.create_db(self.args.gtf, db_filename)
                self._done()

        sys.stdout.write('{}\tGetting junction-direction-exon triples for '
                         'graph database ...\n'.format(util.timestamp()))
        junction_annotator = junctions.JunctionAnnotator(splice_junctions, db)
        junction_exon_triples = junction_annotator.get_adjacent_exons()
        self._done()

        sys.stdout.write('{}\tPopulating graph database of the '
                         'junction-direction-exon triples '
                         '...'.format(util.timestamp()))
        event_maker = events.EventMaker(junction_exon_triples, db)
        self._done()

        for name, abbrev in events.EVENT_TYPES:
            name_with_spaces = name.replace('_', ' ')
            # Find event junctions
            sys.stdout.write(
                '{timestamp}\tFinding all {name} ({abbrev}) event ...'
                '\n'.format(timestamp=util.timestamp(),
                            name=name_with_spaces, abbrev=abbrev.upper()))
            events_of_type = getattr(event_maker, name)()
            self._done()

            # Write to a file
            csv = os.path.join(self.args.index, ['index', abbrev, 'junctions'])
            sys.stdout.write('{timestamp}\tWriting {abbrev} events to {csv} '
                             '...\n'.format(timestamp=util.timestamp(),
                                            abbrev=abbrev.upper(), csv=csv))
            events_of_type.to_csv(csv, index=False)
            self._done()

    def psi(self):
        """Calculate percent spliced in (psi) of splicing events"""

        logger = logging.getLogger('outrigger.psi')
        if self.args.debug:
            logger.setLevel(10)

        try:
            sys.stdout.write(
                '{}\tReading splice junction reads from {} ...'
                '\n'.format(util.timestamp(), self.args.splice_junction_csv))
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

            sys.stdout.write('{}\t\tDone.\n'.format(util.timestamp()))
        except KeyError:
            sys.stdout.write('{}\tCreating consolidated splice junction '
                             'file "{sj_csv}" from SJ.out.tab files ...'
                             '\n'.format(util.timestamp(),
                                         sj_csv=SPLICE_JUNCTIONS_CSV))
            splice_junction_reads = self.csv()
            self._done()

        splice_junction_reads = splice_junction_reads.set_index(
            [self.args.junction_location_col, self.args.sample_id_col])
        splice_junction_reads.sort_index(inplace=True)
        logger.debug('\n--- Splice Junction reads ---')
        logger.debug(repr(splice_junction_reads.head()))

        events_folder = os.path.join(self.args.index, 'events')
        psis = []
        for filename in glob.iglob('{}/*.csv'.format(events_folder)):
            event_type = os.path.basename(filename).split('.csv')[0]
            sys.stdout.write('{}\tReading {} events from {} ...\n'.format(
                util.timestamp(), event_type.upper(), filename))

            isoform_junctions = psi.ISOFORM_JUNCTIONS[event_type]
            event_annotation = pd.read_csv(filename, index_col=0)
            self._done()
            logger.debug('\n--- Splicing event annotation ---')
            logger.debug(repr(event_annotation.head()))

            sys.stdout.write('{}\tCalculating percent spliced-in (Psi) '
                             'scores on {} events ...\n'.format(
                                util.timestamp(), event_type.upper()))
            event_psi = psi.calculate_psi(
                event_annotation, splice_junction_reads,
                min_reads=self.args.min_reads, debug=self.args.debug,
                reads_col=self.args.reads_col, **isoform_junctions)
            event_psi.to_csv(os.path.join(self.args.index,
                                          '{}_psi.csv'.format(event_type)))
            psis.append(event_psi)
            sys.stdout.write('{}\t\tDone.\n'.format(util.timestamp()))

        sys.stdout.write('{}\tConcatenating all calculated psi scores '
                         'into one big matrix...\n'.format(util.timestamp()))
        splicing = pd.concat(psis, axis=1)
        sys.stdout.write('{}\t\tDone.\n')
        splicing = splicing.T
        csv = os.path.join(self.args.index, 'psi.csv')
        sys.stdout.write('{}\tWriting a samples x features matrix of Psi '
                         'scores to {} ...\n'.format(util.timestamp(), csv))
        splicing.to_csv(csv)
        sys.stdout.write('{}\t\tDone.\n'.format(util.timestamp()))


def main():
    CommandLine(sys.argv[1:])


if __name__ == '__main__':
    # try:
    main()
