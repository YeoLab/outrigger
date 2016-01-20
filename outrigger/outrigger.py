#!/usr/bin/env python

import argparse
import datetime
import os
import sys

import gffutils


class CommandLine(object):
    def __init__(self, inOpts=None):
        self.parser = argparse.ArgumentParser(
            description='Calculate percent-spliced in on a *de novo* '
                        'custom-built splicing index')
        self.subparser = self.parser.add_subparsers(help='Sub-commands')



        # # parse_args defaults to [1:] for args, but you need to
        # # exclude the rest of the args too, or validation will fail
        # args = self.parser.parse_args(sys.argv[1:2])
        # if not hasattr(self, args.command):
        #     sys.stderr.write('Unrecognized command\n')
        #     self.parser.print_help()
        #     exit(1)
        # # use dispatch pattern to invoke method with same name
        # getattr(self, args.command)()

    # def build(self):
        """Subcommand to build the index of splicing events"""
        build_parser = self.subparser.add_parser(
            'build', help='Build an index of splicing events using a graph '
                          'database on your junction reads and an annotation')

        build_parser.add_argument('-j', '--sj-out-tab', required=True,
                                 type=str, action='store',
                                 help='SJ.out.tab files from STAR aligner output')
        gtf = build_parser.add_mutually_exclusive_group(required=True)
        gtf.add_argument('-g', '--gtf', type=str, action='store',
                              help="Name of the gtf file you want to use. If a "
                                   "gffutils feature database doesn't already exist"
                                   " at this location plus '.db' (e.g. if your gtf"
                                   " is gencode.v19.annotation.gtf, then the "
                                   "database is inferred to be gencode.v19."
                                   "annotation.gtf.db), then a database will be "
                                   "auto-created")
        gtf.add_argument('-db', '--gffutils-database', type=str,
                              action='store',
                              help="Name of the gffutils database file you want to "
                                   "use. The exon IDs defined here will be used in "
                                   "the function")
        build_parser.add_argument('-o', '--output', required=False, type=str,
                                 action='store', default='./outrigger_output',
                                 help='Where you want to save the generated files')

    # def psi(self):
        """Subcommand to calculate psi on the built index"""
        psi_parser = self.subparser.add_parser(
            'psi', help='Calculate "percent spliced-in" (Psi) values using the '
                        'splicing event index built with "outrigger index"')
        psi_parser.add_argument('-j')

        if inOpts is None:
            self.args = self.parser.print_usage()
        else:
            self.args = vars(self.parser.parse_args(inOpts))

    # def do_usage_and_die(self, str):
    #     '''
    #     If a critical error is encountered, where it is suspected that the
    #     program is not being called with consistent parameters or data, this
    #     method will write out an error string (str), then terminate execution
    #     of the program.
    #     '''
    #     sys.stderr.write(str)
    #     self.parser.print_usage()
    #     return 2


# Class: Usage
class Usage(Exception):
    '''
    Used to signal a Usage error, evoking a usage statement and eventual
    exit when raised
    '''

    def __init__(self, msg):
        self.msg = msg


if __name__ == '__main__':
    # try:
    cl = CommandLine()

        # gtf = cl.args['gtf']
        # db_filename = gtf + '.db'
        # if gtf is not None:
        #     if os.path.exists(db_filename):
        #         db = gffutils.FeatureDB(db_filename)
        #     else:
        #         sys.stdout.write('Starting writing the gffutils FeatureDB of '
        #                          '{} ...\n\tStart:{}\n'.format(
        #             gtf, datetime.datetime.now()))
        #         from .gtf import create_db
        #         db = create_db(gtf, db_filename)
        #         sys.stdout.write('\tDone:{}\n'.format(
        #             datetime.datetime.now()))
        # else:
        #     db = gffutils.FeatureDB(cl.args['db'])
        #
        # splice_junctions = pd.read_csv(cl.args['junctions'],
        #                                sep=cl.args['sep'])
        #
        # ja = JunctionAnnotator(splice_junctions, db, **cl.args)

    # except Usage as err:
    #     cl.do_usage_and_die()
