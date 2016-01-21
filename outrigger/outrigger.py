#!/usr/bin/env python

import argparse
import datetime
import os
import sys

import gffutils


class CommandLine(object):
    def __init__(self, inOpts=None):
        self.parser = argparse.ArgumentParser(
            description='Calculate "percent-spliced in" (Psi) scores of '
                        'alternative splicing on a *de novo*, custom-built '
                        'splicing index')
        self.subparser = self.parser.add_subparsers(help='Sub-commands')


    # def build(self):
        """Subcommand to build the index of splicing events"""
        build_parser = self.subparser.add_parser(
            'build', help='Build an index of splicing events using a graph '
                          'database on your junction reads and an annotation')

        build_parser.add_argument('-j', '--sj-out-tab', required=True,
                                 type=str, action='store', nargs='*',
                                 help='SJ.out.tab files from STAR aligner '
                                      'output')
        gtf = build_parser.add_mutually_exclusive_group(required=True)
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
        build_parser.add_argument('-o', '--output', required=False, type=str,
                                 action='store', default='./outrigger_index',
                                 help='Where you want to save the index of '
                                      'splicing events')

    # def psi(self):
        """Subcommand to calculate psi on the built index"""
        psi_parser = self.subparser.add_parser(
            'psi', help='Calculate "percent spliced-in" (Psi) values using the '
                        'splicing event index built with "outrigger index"')
        psi_parser.add_argument('-i', '--index', required=True,
                                help='Name of the folder where you saved the '
                                     'output from "outrigger index" (default '
                                     'is ./outrigger_index, which is relative '
                                     'to the directory where you called the '
                                     'program)"')
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
                 'you specify')

        if inOpts is None or len(inOpts) == 0:
            self.args = self.parser.print_usage()
        else:
            self.args = vars(self.parser.parse_args(inOpts))
        print(self.args)

        # # parse_args defaults to [1:] for args, but you need to
        # # exclude the rest of the args too, or validation will fail
        # args = self.parser.parse_args(sys.argv[1:2])
        # if not hasattr(self, self.args[0]):
        #     sys.stderr.write('Unrecognized command\n')
        #     self.parser.print_help()
        #     exit(1)
        # # use dispatch pattern to invoke method with same name
        # getattr(self, args.command)()

    def build(self):
        pass

    def psi(self):
        pass



if __name__ == '__main__':
    # try:
    cl = CommandLine(sys.argv[1:])
