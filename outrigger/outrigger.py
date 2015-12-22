import argparse
import datetime
import os

import gffutils


class CommandLine(object):
    def __init__(self, inOpts=None):
        self.parser = parser = argparse.ArgumentParser(
            description='Annotate splice junctions with adjacent exons by '
                        'creating a CSV of exon,direction,junction triples '
                        '(ready for a graph database)')
        parser.add_argument('-j', '--junctions', required=True,
                            type=str, action='store',
                            help='Table of splice junctions, which must have '
                                 'the columns specified by "--exon-start", '
                                 '"--exon-stop", "--chrom", "--junction-id", '
                                 'and "--strand"')
        parser.add_argument('-s', '--sep', type=str, action='store',
                            default=',',
                            help="Separator between items in '--junctions'. "
                                 "Default is a comma, indicating a "
                                 "comma-separated variable (CSV) file. For "
                                 "tab-separated, specify '\\t'")
        gtf = parser.add_mutually_exclusive_group(required=True)
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
        parser.add_argument('-o', '--output', required=False, type=str,
                            action='store', default='./outrigger_output',
                            help='Where you want to save the generated files')
        parser.add_argument('--chrom', type=str, action='store',
                            default='chrom',
                            help="Name of the column in the file specified by "
                                 "'--junctions' which contains the "
                                 "chromosome name of the location of the "
                                 "junctions")
        parser.add_argument('--exon-start', type=str, action='store',
                            default='exon_start',
                            help="Name of the column in the file specified by "
                                 "'--junctions' which contains the "
                                 "genomic start site of where exons would be, "
                                 "relative to the junction")
        parser.add_argument('--exon-stop', type=str, action='store',
                            default='exon_stop',
                            help="Name of the column in the file specified by "
                                 "'--splice-junctions' which contains the "
                                 "genomic stop site of where exons would be, "
                                 "relative to the junction")
        parser.add_argument('--strand', type=str, action='store',
                            default='strand',
                            help="Name of the column in the file specified by "
                                 "'--splice-junctions' which contains the "
                                 "strand of the junction")
        parser.add_argument('--junction-id', type=str, action='store',
                            default='junction_id',
                            help="Name of the column in the file specified by "
                                 "'--splice-junctions' which contains the "
                                 "unique identifier of the junction")
        if inOpts is None:
            self.args = vars(self.parser.parse_args())
        else:
            self.args = vars(self.parser.parse_args(inOpts))

    def do_usage_and_die(self, str):
        '''
        If a critical error is encountered, where it is suspected that the
        program is not being called with consistent parameters or data, this
        method will write out an error string (str), then terminate execution
        of the program.
        '''
        import sys

        print >> sys.stderr, str
        self.parser.print_usage()
        return 2


# Class: Usage
class Usage(Exception):
    '''
    Used to signal a Usage error, evoking a usage statement and eventual
    exit when raised
    '''

    def __init__(self, msg):
        self.msg = msg


if __name__ == '__main__':
    try:
        cl = CommandLine()

        gtf = cl.args['gtf']
        db_filename = gtf + '.db'
        if gtf is not None:
            if os.path.exists(db_filename):
                db = gffutils.FeatureDB(db_filename)
            else:
                sys.stdout.write('Starting writing the gffutils FeatureDB of '
                                 '{} ...\n\tStart:{}\n'.format(
                    gtf, datetime.datetime.now()))
                from .gtf import create_db
                db = create_db(gtf, db_filename)
                sys.stdout.write('\tDone:{}\n'.format(
                    datetime.datetime.now()))
        else:
            db = gffutils.FeatureDB(cl.args['db'])

        splice_junctions = pd.read_csv(cl.args['junctions'],
                                       sep=cl.args['sep'])

        ja = JunctionAnnotator(splice_junctions, db, **cl.args)

    except Usage, err:
        cl.do_usage_and_die()
