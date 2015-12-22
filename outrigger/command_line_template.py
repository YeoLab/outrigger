#!/usr/bin/env python

__author__ = 'olga'

import argparse
import os


class CommandLine(object):
    def __init__(self, inOpts=None):
        self.parser = parser = argparse.ArgumentParser(
            description='Generate a STAR genome index from fasta files')
        parser.add_argument('-g', '--genomeFastaFiles', required=True,
                            type=str, action='store',
                            help='Fasta files of the genome you want to index using '
                                 'STAR.')
        parser.add_argument('-o', '--genomeDir', required=True, type=str,
                            action='store',
                            help='Where you want to save the generated genome')
        parser.add_argument('-n', '--name', default='genomeGenerate',
                            action='store', type=str,
                            help='The name of the submitted job in the queue')
        sjdb = parser.add_mutually_exclusive_group(required=False)
        sjdb.add_argument('--sjdbFileChrStartEnd', default='',
                          type=str, action='store',
                          help='A bed-file-like splice junction file, for example the '
                               'SJ.out.tab file produced by STAR')
        sjdb.add_argument('--sjdbGTFfile', default='',
                          type=str, action='store',
                          help='A GTF file to create a splice junction database from')
        parser.add_argument('--sjdbOverhang', default=100, type=str,
                            action='store',
                            help='Number of bases to overhang for the splice '
                                 'junctions. Ideally should be the (length of '
                                 'one read)-1')
        parser.add_argument('--out-sh', action='store', type=str,
                            required=False,
                            help='The sh file written and submitted to the '
                                 'cluster')
        parser.add_argument('--do-not-submit', required=False,
                            action='store_true', default=False,
                            help='Flag to not actually submit the job but '
                                 'just write the sh file (for testing)')
        parser.add_argument('--queue-type', required=False, type=str,
                            action='store', default='PBS',
                            help='Type of the queue to submit to. For testing '
                                 'purposes on non-server devices, e.g. laptops')
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


class GenomeGenerate(object):
    def __init__(self, genomeDir, genomeFastaFiles, sjdb,
                 sjdbOverhang, job_name, out_sh=None, submit=True):
        """Any CamelCase here is directly copied from the STAR inputs for
        complete compatibility
        """
        # Make the directory if it's not there already
        try:
            os.mkdir(genomeDir)
        except OSError:
            pass

        commands = []
        commands.append('STAR --runMode genomeGenerate --genomeDir {0} '
                        '--genomeFastaFiles {1} --runThreadN 16 {2} '
                        '--sjdbOverhang {3}'.format(
            genomeDir, genomeFastaFiles, sjdb, sjdbOverhang))

        sub = Submitter(queue_type='PBS', sh_filename=out_sh,
                        commands=commands,
                        job_name=job_name, nodes=1, ppn=16, queue='home',
                        walltime='4:00:00')
        sub.job(submit=submit)


if __name__ == '__main__':
    try:
        cl = CommandLine()

        job_name = cl.args['name']
        out_sh = name = job_name + '.sh' if cl.args['out_sh'] is None \
            else cl.args['out_sh']
        submit = not cl.args['do_not_submit']

        sjdb_arguments = ['sjdbGTFfile', 'sjdbFileChrStartEnd']

        sjdb = ''.join('--{} {}'.format(k, cl.args[k]) for k in sjdb_arguments
                       if cl.args[k])
        GenomeGenerate(cl.args['genomeDir'], cl.args['genomeFastaFiles'], sjdb,
                       cl.args['sjdbOverhang'], job_name, out_sh,
                       submit=submit)
    except Usage, err:
        cl.do_usage_and_die()
