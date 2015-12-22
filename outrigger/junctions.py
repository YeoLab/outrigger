import datetime
import sys

import gffutils
import pandas as pd

UPSTREAM = 'upstream'
DOWNSTREAM = 'downstream'
DIRECTIONS = UPSTREAM, DOWNSTREAM


class JunctionAnnotator(object):
    """Annotate junctions with adjacent exons"""

    def __init__(self, splice_junctions, db, junction_id='junction_id',
                 exon_start='exon_start', exon_stop='exon_stop',
                 chrom='chrom', strand='strand'):
        """Initialize class to get upstream/downstream exons of junctions

        Parameters
        ----------
        splice_junctions : pandas.DataFrame
            A table of splice junctions with the columns indicated by the
            variables `junction_id`, `exon_start`, `exon_stop`, `chrom`,
            `strand`
        db : gffutils.FeatureDB
            Gffutils Database of gene, transcript, and exon features.
        junction_id, exon_start, exon_stop, chrom, strand : str
            Columns in `splice_junctions`
        """
        columns = self.junction_id, self.exon_start, self.exon_stop, \
                  self.chrom, self.strand

        for column in columns:
            if column not in splice_junctions:
                raise ValueError('The required column {} is not in the splice '
                                 'junction dataframe'.format(column))

        self.splice_junctions = splice_junctions.set_index(junction_id)
        self.splice_junctions = self.splice_junctions.sort_index()

        self.junction_id = junction_id
        self.exon_start = exon_start
        self.exon_stop = exon_stop
        self.chrom = chrom
        self.strand = strand

        self.db = db

    @staticmethod
    def _single_junction_exon_triple(direction_ind, direction, exon_id):
        """Create exon, direction, junction triples for an exon + its junctions

        Parameters
        ----------
        direction_ind : pandas.Series
            A boolean series of the indices of the junctions matching with the
            provided exon. The index of the series must be the junction ID
        direction : str
            The direction of the exon relative to the junctions, either
            "upstream" or "downstream"
        exon_id : str
            Unique identifier of the exon

        Returns
        -------


        """
        length = direction_ind.sum()

        exons = [exon_id] * length
        directions = [direction] * length
        junctions = direction_ind[direction_ind].index
        return pd.DataFrame(zip(exons, directions, junctions),
                            columns=['exon', 'direction', 'junction'])

    @staticmethod
    def genome_to_transcript_adjacency(adjacent_in_genome, strand):
        """If negative strand, swap the upstream/downstream adjacency"""
        if strand == '+':
            return {UPSTREAM: adjacent_in_genome[UPSTREAM],
                    DOWNSTREAM: adjacent_in_genome[DOWNSTREAM]}
        elif strand == '-':
            return {UPSTREAM: adjacent_in_genome[DOWNSTREAM],
                    DOWNSTREAM: adjacent_in_genome[UPSTREAM]}


    def genome_adjacent(self, exon, exon_start='exon_start',
                        exon_stop='exon_stop', chrom='chrom', strand='strand'):
        """Get indices of junctions next to an exon, in genome coordinates"""
        chrom_ind = self.splice_junctions[chrom] == exon.chrom

        strand_ind = self.splice_junctions[strand] == exon.strand

        upstream_in_genome = chrom_ind & strand_ind \
                             & (self.splice_junctions[exon_stop] == exon.stop)
        downstream_in_genome = chrom_ind & strand_ind \
                               & (self.splice_junctions[exon_start] == exon.start)
        return {UPSTREAM: upstream_in_genome, DOWNSTREAM: downstream_in_genome}


    def adjacent_junctions(self, exon, exon_start='exon_start',
                           exon_stop='exon_stop', chrom='chrom',
                           strand='strand'):
        """Get junctions adjacent to this exon"""
        dfs = []
        adjacent_in_genome = self.genome_adjacent(exon, self.splice_junctions,
                                                  exon_start, exon_stop, chrom,
                                                  strand)
        adjacent_in_transcriptome = self.genome_to_transcript_adjacency(
            adjacent_in_genome, exon.strand)

        exon_id = exon.id
        for direction, ind in adjacent_in_transcriptome.items():
            if ind.any():
                df = self._single_junction_exon_triple(ind, direction, exon_id)
                dfs.append(df)

        if len(dfs) > 0:
            return pd.concat(dfs, ignore_index=True)
        else:
            return pd.DataFrame()

    def get_adjacent_exons(self, exon_start='exon_start',
                           exon_stop='exon_stop', chrom='chrom'):
        """Get upstream and downstream exons of each junction

        Use junctions defined in ``sj_metadata`` and exons in ``db`` to create
        triples of (exon, direction, junction), which are read like
        (subject, object, verb) e.g. ('exon1', 'upstream', 'junction12'), for
        creation of a graph database.

        Parameters
        ----------
        sj_metadata : pandas.DataFrame
            A splice junction metadata dataframe with the junction id as the
            index, with  columns defined by variables ``exon_start`` and
            ``exon_stop``.
        db : gffutils.FeatureDB
            A database of gene annotations created by gffutils. Must have features
            of type "exon"
        exon_start : str, optional
            Name of the column in sj_metadata corresponding to the start of the
            exon
        exon_stop : str, optional
            Name of the column in sj_metadata corresponding to the end of the exon

        Returns
        -------
        junction_exon_triples : pandas.DataFrame
            A three-column dataframe describing the relationship of where an exon
            is relative to junctions
        """
        n_exons = sum(1 for _ in self.db.features_of_type('exon'))

        dfs = []

        sys.stdout.write('Starting annotation of all junctions with known '
                         'exons...\n')
        for i, exon in enumerate(self.db.features_of_type('exon')):
            if (i + 1) % 10000 == 0:
                sys.stdout.write('\t{}/{} exons completed\n'.format(i + 1,
                                                                    n_exons))
            df = self.adjacent_junctions(exon)
        junction_exon_triples = pd.concat(dfs, ignore_index=True)
        sys.stdout.write('Done.\n')
        return junction_exon_triples


import argparse
import os


class CommandLine(object):
    def __init__(self, inOpts=None):
        self.parser = parser = argparse.ArgumentParser(
            description='Annotate splice junctions with adjacent exons by '
                        'creating a CSV of exon,direction,junction triples '
                        '(ready for a graph database)')
        parser.add_argument('-j', '--junctions', required=True,
                            type=str, action='store',
                            help='Table of splice junctions, which must have the'
                                 ' columns specified by "--exon-start", '
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
                            action='store', default='./',
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
