#!/usr/bin/env python

import argparse
import glob
import os
import sys
import warnings

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import pandas as pd

from outrigger import psi
from outrigger import star
from outrigger import util

class CommandLine(object):
    def __init__(self, input_options=None):
        self.parser = argparse.ArgumentParser(
            description='Calculate "percent-spliced in" (Psi) scores of '
                        'alternative splicing on a *de novo*, custom-built '
                        'splicing index')

        self.subparser = self.parser.add_subparsers(help='Sub-commands')

        # Subcommand to build the index of splicing events
        index_parser = self.subparser.add_parser(
            'index', help='Build an index of splicing events using a graph '
                          'database on your junction reads and an annotation')

        index_parser.add_argument('-i', '--index', required=False, type=str,
                                  action='store', default='./outrigger_index',
                                  help='Name of the folder where you saved the '
                                     'output from "outrigger index" (default '
                                     'is ./outrigger_index, which is relative '
                                     'to the directory where you called the '
                                     'program)". You will need this file for '
                                       'the next step, "outrigger psi"')

        junctions = index_parser.add_mutually_exclusive_group(required=True)
        junctions.add_argument(
            '-j', '--sj-out-tab', type=str, action='store',
            nargs='*', help='SJ.out.tab files from STAR aligner output')
        junctions.add_argument('-c', '--splice-junction-csv',
                               help="Name of the splice junction files to calculate psi scores "
                                    "on. If not provided, the compiled 'sj.csv' file with all the"
                                    " samples from the SJ.out.tab files that were used during "
                                    "'outrigger index' will be used. Not required if you specify "
                                    "SJ.out.tab file with '--sj-out-tab'")

        index_parser.add_argument('-m', '--min-reads', type=int, action='store',
                                  required=False, default=10,
                                  help='Minimum number of reads per junction for '
                                    'that junction to count in creating the '
                                    'index of splicing events (default=10)')

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

        # Subcommand to calculate psi on the built index
        psi_parser = self.subparser.add_parser(
            'psi', help='Calculate "percent spliced-in" (Psi) values using the '
                        'splicing event index built with "outrigger index"')
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
                 "on. If not provided, the compiled 'sj.csv' file with all the"
                 " samples from the SJ.out.tab files that were used during "
                 "'outrigger index' will be used. Not required if you specify "
                 "SJ.out.tab file with '--sj-out-tab'")
        splice_junctions.add_argument(
            '-j', '--sj-out-tab', required=False,
            type=str, action='store', nargs='*',
            help='SJ.out.tab files from STAR aligner output. Not required if '
                 'you specify a file with "--splice-junction-csv"')
        psi_parser.add_argument('-m', '--min-reads', type=int, action='store',
                                required=False,
                                help='Minimum number of reads per junction for '
                                     'calculating Psi (default=10)')
        psi_parser.add_argument('--debug', required=False, action='store_true',
                                help='If given, print debugging logging '
                                     'information to standard out')

        print(input_options)
        if input_options is None or len(input_options) == 0:
            self.parser.print_usage()
            self.args = None
        else:
            self.args = vars(self.parser.parse_args(input_options))
        print(self.args)

    def csv(self):
        """Create a csv file of compiled splice junctions

        Parameters
        ----------
        var1 : array_like
            Array_like means all those objects -- lists, nested lists, etc. --
            that can be converted to an array.  We can also refer to
            variables like `var1`.
        var2 : int
            The type above can either refer to an actual Python type
            (e.g. ``int``), or describe the type of the variable in more
            detail, e.g. ``(N,) ndarray`` or ``array_like``.
        Long_variable_name : {'hi', 'ho'}, optional
            Choices in brackets, default first when optional.

        Returns
        -------
        type
            Explanation of anonymous return value of type ``type``.
        describe : type
            Explanation of return value named `describe`.
        out : type
            Explanation of `out`.
        """
        splice_junctions = star.read_multiple_sj_out_tab(self.args['sj_out_tab'])
        splice_junctions.to_csv(os.path.join(self.args['index'], 'sj.csv'), index=False)
        return splice_junctions

    def index(self):
        pass

    def psi(self):

        try:
            csv = self.args['splice_junction_csv']
            sys.stdout.write('{}\tReading splice junction reads from {} ...\n'.format(
                util.timestamp(), csv))
            splice_junction_reads = pd.read_csv(csv)
            sys.stdout.write('{}\t\tDone.\n'.format(util.timestamp()))
        except KeyError:
            sys.stdout.write('{}\tCreating consolidated splice junction file "sj.csv" '
                             'from SJ.out.tab files ...\n'.format(util.timestamp()))
            splice_junction_reads = self.csv()
            sys.stdout.write('{}\t\tDone.\n'.format(util.timestamp()))

        splice_junction_reads = splice_junction_reads.set_index(
            ['junction_location', 'sample_id'])

        events_folder = os.path.join(self.args['index'], 'events')
        psis = []
        for filename in glob.iglob('{}/*.csv'.format(events_folder)):
            event_type = os.path.basename(filename).split('.csv')[0]
            sys.stdout.write('{}\tReading {} events from {} ...\n'.format(
                util.timestamp(), event_type.upper(), filename))

            isoform_junctions = psi.ISOFORM_JUNCTIONS[event_type]
            event_annotation = pd.read_csv(filename)
            sys.stdout.write('{}\t\tDone.\n'.format(util.timestamp()))

            sys.stdout.write('{}\tCalculating percent spliced-in (Psi) '
                             'scores on {} events ...\n'.format(
                util.timestamp(), event_type.upper()))
            event_psi = psi.calculate_psi(event_annotation, splice_junction_reads,
                                      min_reads=self.args['min_reads'],
                                      debug=self.args['debug'],
                                      **isoform_junctions)
            sys.stdout.write('{}\t\tDone.\n'.format(util.timestamp()))
            psis.append(event_psi)

        sys.stdout.write('{}\tConcatenating all calculated psi scores '
                         'into one big matrix...\n'.format(util.timestamp()))
        splicing = pd.concat(psis)
        splicing = splicing.T
        csv = os.path.join(self.args['index'], 'psi.csv')
        sys.stdout.write('{}\tWriting a samples x features matrix of Psi '
                         'scores to {} ...\n'.format(util.timestamp(), csv))
        splicing.to_csv(csv)
        sys.stdout.write('{}\t\tDone.\n'.format(util.timestamp()))

def main():
    cl = CommandLine(sys.argv[1:])

if __name__ == '__main__':
    # try:
    main()
