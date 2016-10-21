#!/usr/bin/env python

import argparse
import logging
import os
import pdb
import shutil
import sys
import traceback

import gffutils
import numpy as np
import pandas as pd

from outrigger import __version__
import outrigger.common
from outrigger import util, common
from outrigger.index import events, adjacencies
from outrigger.io import star, gtf
from outrigger.psi import compute
from outrigger.validate import check_splice_sites


OUTPUT = os.path.join('.', 'outrigger_output')
JUNCTION_PATH = os.path.join(OUTPUT, 'junctions')
JUNCTION_READS_PATH = os.path.join(JUNCTION_PATH, 'reads.csv')
JUNCTION_METADATA_PATH = os.path.join(JUNCTION_PATH, 'metadata.csv')
INDEX = os.path.join(OUTPUT, 'index')
EVENTS_CSV = 'events.csv'
METADATA_CSV = 'metadata.csv'


class CommandLine(object):
    def __init__(self, input_options=None):
        self.parser = argparse.ArgumentParser(
            description='outrigger {version}\nCalculate "percent-spliced in" '
                        '(Psi) scores of alternative splicing on a *de novo*, '
                        'custom-built splicing index -- '
                        'just for you!'.format(version=__version__))
        self.parser.add_argument(
            '--version', action='version',
            version='outrigger {version}'.format(version=__version__))
        self.subparser = self.parser.add_subparsers(help='Sub-commands')

        # --- Subcommand to build the index of splicing events --- #
        index_parser = self.subparser.add_parser(
            'index', help='Build an index of splicing events using a graph '
                          'database on your junction reads and an annotation')

        index_parser.add_argument(
            '-o', '--output', required=False, type=str, action='store',
            default=OUTPUT,
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
            '-c', '--junction-reads-csv', required=False,
            help="Name of the splice junction files to calculate psi scores "
                 "on. If not provided, the compiled '{sj_csv}' file with all "
                 "the samples from the SJ.out.tab files that were used during "
                 "'outrigger index' will be used. Not required if you specify "
                 "SJ.out.tab file with '--sj-out-tab'".format(
                        sj_csv=JUNCTION_READS_PATH))
        index_parser.add_argument('-m', '--min-reads', type=int,
                                  action='store',
                                  required=False, default=10,
                                  help='Minimum number of reads per junction '
                                       'for that junction to count in creating'
                                       ' the index of splicing events '
                                       '(default=10)')
        index_parser.add_argument('--ignore-multimapping', action='store_true',
                                  help='Applies to STAR SJ.out.tab files only.'
                                       ' If this flag is used, then do not '
                                       'include reads that mapped to multiple '
                                       'locations in the genome, not uniquely '
                                       'to a locus, in the read count for a '
                                       'junction. By default, this is off, and'
                                       ' all reads are used.')
        index_parser.add_argument(
            '-l', '--max-de-novo-exon-length',
            default=adjacencies.MAX_DE_NOVO_EXON_LENGTH, action='store',
            help='Maximum length of an exon detected '
                 '*de novo* from the dataset. This is to'
                 ' prevent multiple kilobase long exons '
                 'from being accidentally created. '
                 '(default={})'.format(adjacencies.MAX_DE_NOVO_EXON_LENGTH))

        gtf_parser = index_parser.add_mutually_exclusive_group(required=True)
        gtf_parser.add_argument('-g', '--gtf-filename', type=str,
                                action='store',
                                help="Name of the gtf file you want to use. If"
                                     " a gffutils feature database doesn't "
                                     "already exist at this location plus "
                                     "'.db' (e.g. if your gtf is gencode.v19."
                                     "annotation.gtf, then the database is "
                                     "inferred to be gencode.v19.annotation."
                                     "gtf.db), then a database will be auto-"
                                     "created. Not required if you provide a "
                                     "pre-built database with "
                                     "'--gffutils-db'")
        gtf_parser.add_argument('-d', '--gffutils-db', type=str,
                                action='store',
                                help="Name of the gffutils database file you "
                                     "want to use. The exon IDs defined here "
                                     "will be used in the function when "
                                     "creating splicing event names. Not "
                                     "required if you provide a gtf file with"
                                     " '--gtf-filename'")
        index_parser.add_argument('--debug', required=False,
                                  action='store_true',
                                  help='If given, print debugging logging '
                                       'information to standard out (Warning:'
                                       ' LOTS of output. Not recommended '
                                       'unless you think something is going '
                                       'wrong)')
        overwrite_parser = index_parser.add_mutually_exclusive_group(
            required=False)
        overwrite_parser.add_argument('--force', action='store_true',
                                      help="If the 'outrigger index' command "
                                           "was interrupted, there will be "
                                           "intermediate files remaining. If "
                                           "you wish to restart outrigger and "
                                           "overwrite them all, use this flag."
                                           " If you want to continue from "
                                           "where you left off, use the "
                                           "'--resume' flag. If neither is "
                                           "specified, the program exits and "
                                           "complains to the user.")
        overwrite_parser.add_argument('--resume', action='store_true',
                                      help="If the 'outrigger index' command "
                                           "was interrupted, there will be "
                                           "intermediate files remaining. If "
                                           "you want to continue from "
                                           "where you left off, use this flag."
                                           " The default action is to do "
                                           "nothing and ask the user for "
                                           "input.")
        index_parser.set_defaults(func=self.index)

        # -- Subcommand to validate exons of built index --- #
        validate_parser = self.subparser.add_parser(
            'validate', help='Ensure that the splicing events found all have '
                             'the correct splice sites')
        validate_parser.add_argument('-f', '--fasta', required=True,
                                     help='Location of the genome fasta file '
                                          'for which to get the splice site '
                                          'sequences from')
        validate_parser.add_argument('-g', '--genome', required=True,
                                     help='Either the genome name (e.g. '
                                          '"mm10" or "hg19") or location of '
                                          'the genome chromosome '
                                          'sizes file for "bedtools flank" to '
                                          'make sure we do not accidentally '
                                          'ask for genome positions that are '
                                          'outside of the defined range')
        validate_parser.add_argument('-i', '--index', required=False,
                                     default=INDEX,
                                     help='Name of the folder where you saved '
                                          'the output from "outrigger index" '
                                          '(default is {}, which is relative '
                                          'to the directory where you called '
                                          'this program, assuming you have '
                                          'called "outrigger psi" in the same '
                                          'folder as you called "outrigger '
                                          'index")'.format(INDEX))
        validate_parser.add_argument(
            '-s', '--valid-splice-sites', required=False,
            default=check_splice_sites.MAMMALIAN_SPLICE_SITES,
            help="The intron-definition based splice sites that are allowed in"
                 " the data, which is in the format 5'/3' of the intron, and "
                 "separated by commas for different types. Default is {}, "
                 "which are the major and minor spliceosome splice sites in "
                 "mammalian "
                 "systems.".format(check_splice_sites.MAMMALIAN_SPLICE_SITES))
        validate_parser.add_argument('--debug', required=False,
                                     action='store_true',
                                     help='If given, print debugging logging '
                                          'information to standard out')
        validate_parser.set_defaults(func=self.validate)

        # --- Subcommand to calculate psi on the built index --- #
        psi_parser = self.subparser.add_parser(
            'psi', help='Calculate "percent spliced-in" (Psi) values using '
                        'the splicing event index built with "outrigger '
                        'index"')
        psi_parser.add_argument('-i', '--index', required=False,
                                default=INDEX,
                                help='Name of the folder where you saved the '
                                     'output from "outrigger index" (default '
                                     'is {}, which is relative '
                                     'to the directory where you called this '
                                     'program, assuming you have called '
                                     '"outrigger psi" in the same folder as '
                                     'you called "outrigger '
                                     'index")'.format(INDEX))
        splice_junctions = psi_parser.add_mutually_exclusive_group(
            required=False)
        splice_junctions.add_argument(
            '-c', '--junction-reads-csv', required=False,
            help="Name of the splice junction files to calculate psi scores "
                 "on. If not provided, the compiled '{sj_csv}' file with all "
                 "the samples from the SJ.out.tab files that were used during "
                 "'outrigger index' will be used. Not required if you specify "
                 "SJ.out.tab file with '--sj-out-tab'".format(
                        sj_csv=JUNCTION_READS_PATH))
        splice_junctions.add_argument(
            '-j', '--sj-out-tab', required=False,
            type=str, action='store', nargs='*',
            help='SJ.out.tab files from STAR aligner output. Not required if '
                 'you specify a file with "--junction-reads-csv"')
        psi_parser.add_argument('-m', '--min-reads', type=int, action='store',
                                required=False, default=10,
                                help='Minimum number of reads per junction for'
                                     ' calculating Psi (default=10)')
        psi_parser.add_argument('--ignore-multimapping', action='store_true',
                                help='Applies to STAR SJ.out.tab files only.'
                                     ' If this flag is used, then do not '
                                     'include reads that mapped to multiple '
                                     'locations in the genome, not uniquely '
                                     'to a locus, in the read count for a '
                                     'junction. By default, this is off, and'
                                     ' all reads are used.')
        psi_parser.add_argument('--reads-col', default='reads',
                                help="Name of column in --splice-junction-csv "
                                     "containing reads to use. "
                                     "(default='reads')")
        psi_parser.add_argument('--sample-id-col', default='sample_id',
                                help="Name of column in --splice-junction-csv "
                                     "containing sample ids to use. "
                                     "(default='sample_id')")
        psi_parser.add_argument('--junction-id-col',
                                default='junction_id',
                                help="Name of column in --splice-junction-csv "
                                     "containing the ID of the junction to use"
                                     ". Must match exactly with the junctions "
                                     "in the index."
                                     "(default='junction_id')")
        psi_parser.add_argument('--debug', required=False, action='store_true',
                                help='If given, print debugging logging '
                                     'information to standard out')
        psi_parser.set_defaults(func=self.psi)

        if input_options is None or len(input_options) == 0:
            self.parser.print_usage()
            self.args = None
        else:
            self.args = self.parser.parse_args(input_options)

        if self.args is not None:
            if self.args.debug:
                print(self.args)
                print(input_options)

            self.args.func()

    def index(self):
        index = Index(**vars(self.args))
        index.execute()

    def validate(self):
        validate = Validate(**vars(self.args))
        validate.execute()

    def psi(self):
        psi = Psi(**vars(self.args))
        psi.execute()

    def do_usage_and_die(self, str):
        '''Cleanly exit if incorrect parameters are given

        If a critical error is encountered, where it is suspected that the
        program is not being called with consistent parameters or data, this
        method will write out an error string (str), then terminate execution
        of the program.
        '''
        sys.stderr.write(str)
        self.parser.print_usage()
        type, value, tb = sys.exc_info()
        traceback.print_exc()
        debug = os.getenv('PYTHONDEBUG')
        if debug or self.args.debug:
            print(os.environ)
            # If the PYTHONDEBUG shell environment variable exists, then
            # launch the python debugger
            pdb.post_mortem(tb)

        return 2


# Class: Usage
class Usage(Exception):
    '''
    Used to signal a Usage error, evoking a usage statement and eventual
    exit when raised
    '''

    def __init__(self, msg):
        self.msg = msg


class Subcommand(object):

    output_folder = OUTPUT
    sj_out_tab = None
    junction_reads_csv = JUNCTION_READS_PATH
    ignore_multimapping = False
    min_reads = common.MIN_READS
    reads_col = common.READS
    gtf_filename = None
    gffutils_db = None
    debug = False
    force = False
    resume = False

    def __init__(self, **kwargs):

        # Read all arguments and set as attributes of this class
        for key, value in kwargs.items():
            setattr(self, key, value)

        sys.stdout.write('self.junction_reads_csv: {}\n'.format(
            self.junction_reads_csv))

        # Since you can specify either junction reads csv or sj out tab,
        # the other one might get overwritten as None
        if self.junction_reads_csv is None:
            self.junction_reads_csv = JUNCTION_READS_PATH

        for folder in self.folders:
            self.maybe_make_folder(folder)

    def maybe_make_folder(self, folder):
        util.progress("Creating folder {} ...".format(folder))
        if not os.path.exists(folder):
            os.mkdir(folder)
        util.done()

    @property
    def folders(self):
        return self.output_folder, self.index_folder, self.gtf_folder, \
               self.junctions_folder

    @property
    def index_folder(self):
        return os.path.join(self.output_folder, 'index')

    @property
    def gtf_folder(self):
        return os.path.join(self.index_folder, 'gtf')

    @property
    def junctions_folder(self):
        return os.path.join(self.output_folder, 'junctions')

    def csv(self):
        """Create a csv file of compiled splice junctions"""
        if not os.path.exists(self.junction_reads_csv):
            util.progress(
                'Reading SJ.out.files and creating a big splice junction'
                ' table of reads spanning exon-exon junctions...')
            splice_junctions = star.read_multiple_sj_out_tab(
                self.sj_out_tab, ignore_multimapping=self.ignore_multimapping)

            dirname = os.path.dirname(self.junction_reads_csv)
            if not os.path.exists(dirname):
                os.makedirs(dirname)
            util.progress('Writing {} ...\n'.format(self.junction_reads_csv))
            splice_junctions.to_csv(self.junction_reads_csv, index=False)
            util.done()
        else:
            util.progress('Found compiled junction reads file in {} and '
                          'reading it in ...'.format(self.junction_reads_csv))
            splice_junctions = pd.read_csv(self.junction_reads_csv)
            util.done()
        return splice_junctions

    def maybe_make_db(self):
        """Get GFFutils database from file or create from a gtf"""
        if self.gffutils_db is not None:
            copied_db = os.path.join(self.gtf_folder,
                                     os.path.basename(self.gffutils_db))
            util.progress('Copying gffutils database from {} to {} '
                          '...'.format(self.gffutils_db, copied_db))
            shutil.copyfile(self.gffutils_db, copied_db)
            util.done()

            util.progress('Reading gffutils database from {} ...'.format(
                copied_db))
            db = gffutils.FeatureDB(copied_db)
            util.done()
        else:
            basename = os.path.basename(self.gtf_filename)
            db_filename = os.path.join(self.gtf_folder,
                                       '{}.db'.format(basename))
            util.progress("Found GTF file in {}".format(self.gtf_filename))
            try:
                db = gffutils.FeatureDB(db_filename)
                util.progress(
                    "Found existing built outrigger-built gffutils database "
                    "file in {}".format(db_filename))
            except (ValueError, TypeError):
                util.progress(
                    'Creating a "gffutils" '
                    'database {} ...'.format(db_filename))
                db = gtf.create_db(self.gtf_filename, db_filename)
                util.done()
        return db

    def maybe_overwrite(self, filename):
        """Ensures that filename is not overwritten unless user-specified

        - If the file doesn't exist, return "True", as in "yes, please
          overwrite" even though there wasn't really a file there, but the
          action is still to create the file in the next step
        - If the file exists and the user specified nothing, exit the program
          and complain that either --force or --resume must be specified
        - If the file exists and the user specified --force, then return True
        - If the file exists and the user specified --resume, then return False
          so that the file is not overwritten

        Parameters
        ----------
        filename : str
            Path to a file that you may want to overwrite

        Returns
        -------
        if_overwrite : bool
            If True, then the next step in the program has the go-ahead to
            "overwrite" or create the file. If False, the file exists and the
            user doesn't want to overwrite it
        """
        if not os.path.exists(filename):
            return True
        if os.path.exists(filename):
            if self.force:
                util.progress("Found existing {filename}, overwriting with "
                              "--force flag".format(filename=filename))
                return True
            if self.resume:
                util.progress(
                    "With the flag '--resume', Found an existing file "
                    "containing novel exons,"
                    "{filename}, not re-calculating. To force overwriting, "
                    "use the flag ''--force'.".format(filename=filename))
                return False
        else:
            raise ValueError("Found existing {filename} "
                             "but don't "
                             "know whether you want me to continue where I "
                             "stopped ('--resume') or force overwrite "
                             "and restart from "
                             "scratch ('--force')! Exiting."
                             ".".format(filename=filename))


class Index(Subcommand):

    max_de_novo_exon_length = adjacencies.MAX_DE_NOVO_EXON_LENGTH

    def filter_junctions_on_reads(self, spliced_reads):
        # Filter junction metadata to only get junctions with minimum reads
        util.progress('Filtering for only junctions with minimum {} reads '
                      '...'.format(self.min_reads))
        original = len(spliced_reads.groupby(common.JUNCTION_ID))
        enough_reads_rows = spliced_reads[self.reads_col] >= self.min_reads
        spliced_reads = spliced_reads.loc[enough_reads_rows]
        enough_reads = len(spliced_reads.groupby(common.JUNCTION_ID))
        filtered = original - enough_reads
        util.progress('\t{enough}/{original} junctions remain after '
                      'filtering out {filtered} junctions with < '
                      '{min_reads} '
                      'reads.'.format(filtered=filtered, enough=enough_reads,
                                      original=original,
                                      min_reads=self.min_reads))
        util.done(2)
        return spliced_reads

    @staticmethod
    def junction_metadata(spliced_reads):
        """Get just the junction info from the concatenated read files"""
        util.progress('Creating splice junction metadata of merely where '
                      'junctions start and stop')

        metadata = star.make_metadata(spliced_reads)
        util.done()
        return metadata

    def make_exon_junction_adjacencies(self, metadata, db):
        """Get annotated exon_cols next to junctions in data"""
        exon_junction_adjacencies = adjacencies.ExonJunctionAdjacencies(
            metadata, db, max_de_novo_exon_length=self.max_de_novo_exon_length)

        novel_exons_gtf = os.path.join(self.gtf_folder, 'novel_exons.gtf')
        if self.maybe_overwrite(novel_exons_gtf):
            util.progress('Detecting de novo exons based on gaps between '
                          'junctions ...')
            exon_junction_adjacencies.detect_exons_from_junctions()
            util.done()

            util.progress('Writing novel exons to {} ...'.format(
                novel_exons_gtf))
            exon_junction_adjacencies.write_de_novo_exons(novel_exons_gtf)
            util.done()

        csv = os.path.join(self.index_folder,
                           'junction_exon_direction_triples.csv')
        if not os.path.exists(csv) or self.force:
            util.progress('Getting junction-direction-exon triples for graph '
                          'database ...')
            junction_exon_triples = \
                exon_junction_adjacencies.upstream_downstream_exons()
            util.done()

            util.progress('Writing junction-exon-direction triples'
                          ' to {}...'.format(csv))
            junction_exon_triples.to_csv(csv, index=False)
            util.done()
        elif self.resume:
            junction_exon_triples = pd.read_csv(csv)
        else:
            raise ValueError("Found existing junction-exon-triples file "
                             "({csv}) but don't "
                             "know whether you want me to continue where I "
                             "stopped ('--resume') or force restart from "
                             "scratch ('--force')! Exiting."
                             ".".format(csv=csv))

        return junction_exon_triples

    @staticmethod
    def make_graph(junction_exon_triples, db):
        """Create graph database of exon-junction adjacencies"""
        util.progress('Populating graph database of the '
                      'junction-direction-exon triples ...')

        event_maker = events.EventMaker(junction_exon_triples, db)
        util.done()
        return event_maker

    def make_events_by_traversing_graph(self, event_maker, db):
        for splice_name, splice_abbrev in common.SPLICE_TYPES:
            csv = os.path.join(self.index_folder, splice_abbrev.lower(),
                               EVENTS_CSV)
            if not os.path.exists(csv) or self.force:
                name_with_spaces = splice_name.replace('_', ' ')
                # Find event junctions
                util.progress(
                    'Finding all {name} ({abbrev}) events ...'.format(
                        name=name_with_spaces, abbrev=splice_abbrev.upper()))
                events_of_type = getattr(event_maker, splice_name)()
                util.done()

                # Write to a file

                dirname = os.path.dirname(csv)
                if not os.path.exists(dirname):
                    os.makedirs(dirname)

                n_events = len(events_of_type.groupby(level=0, axis=0))
                if n_events > 0:
                    util.progress(
                        'Found {n} {abbrev} events.'.format(
                            n=n_events, abbrev=splice_abbrev.upper(), csv=csv))
                    self.get_event_attributes(db, events_of_type,
                                              splice_abbrev)
                else:
                    util.progress(
                        'No {abbrev} events found in the junction and exon '
                        'data.'.format(abbrev=splice_abbrev.upper()))
            else:
                util.progress('Found existing {name} ({abbrev}) splicing '
                              'events file ({csv}), so not searching. To force'
                              ' re-finding these splicing events, use the flag'
                              ' "--force".'.format(name=splice_name,
                                                   abbrev=splice_abbrev,
                                                   csv=csv))

    def get_event_attributes(self, db, event_df, splice_type):
        util.progress(
            'Making metadata file of {splice_type} events, '
            'annotating them with GTF attributes ...'.format(
                splice_type=splice_type.upper()))

        sa = gtf.SplicingAnnotator(db, event_df, splice_type.upper())
        util.progress('Making ".bed" files for exons in each event ...')
        folder = os.path.join(self.index_folder, splice_type)
        sa.exon_bedfiles(folder=folder)
        util.done()

        attributes = sa.attributes()
        util.done()

        # Write to a file
        csv = os.path.join(self.index_folder, splice_type, EVENTS_CSV)
        util.progress('Writing {splice_type} events to {csv} '
                      '...'.format(splice_type=splice_type.upper(), csv=csv))
        attributes.to_csv(csv, index=True,
                          index_label=outrigger.common.EVENT_ID_COLUMN)
        util.done()

    def write_new_gtf(self, db):
        gtf = os.path.join(self.gtf_folder,
                           os.path.basename(self.gtf_filename))
        util.progress('Write new GTF to {} ...'.format(gtf))
        with open(gtf, 'w') as f:
            for feature in db.all_features():
                f.write(str(feature) + '\n')
        util.done()

    def execute(self):
        # Must output the junction exon triples
        logger = logging.getLogger('outrigger.index')

        if self.debug:
            logger.setLevel(10)

        spliced_reads = self.csv()
        spliced_reads = self.filter_junctions_on_reads(spliced_reads)

        metadata = self.junction_metadata(spliced_reads)
        metadata_csv = os.path.join(self.junctions_folder, METADATA_CSV)
        util.progress('Writing metadata of junctions to {csv}'
                      ' ...'.format(csv=metadata_csv))
        metadata.to_csv(metadata_csv, index=False)
        util.done()

        db = self.maybe_make_db()

        junction_exon_triples = self.make_exon_junction_adjacencies(
            metadata, db)

        event_maker = self.make_graph(junction_exon_triples, db=db)
        self.make_events_by_traversing_graph(event_maker, db)

        self.write_new_gtf(db)


class Validate(Subcommand):

    def exon_pair_splice_sites(self, exonA, exonB, splice_abbrev):
        name = '{exonA}-{exonB}_splice_site'.format(
            exonA=exonA, exonB=exonB)

        exonA_splice_sites = self.individual_exon_splice_sites(
            exonA, splice_abbrev, 'downstream')
        exonB_splice_sites = self.individual_exon_splice_sites(
            exonB, splice_abbrev, 'upstream')

        intron_splice_site = exonA_splice_sites + '/' \
            + exonB_splice_sites
        intron_splice_site.name = name
        return intron_splice_site

    def individual_exon_splice_sites(self, exon, splice_abbrev, direction):
        exon_bed = os.path.join(self.index_folder, splice_abbrev,
                                '{}.bed'.format(exon))
        splice_sites = check_splice_sites.read_splice_sites(
            exon_bed, self.genome, self.fasta, direction)
        return splice_sites

    def execute(self):
        valid_splice_sites = check_splice_sites.splice_site_str_to_tuple(
            self.valid_splice_sites)

        for splice_name, splice_abbrev in common.SPLICE_TYPES:
            splice_name_spaces = splice_name.replace('_', ' ').title()
            util.progress('Finding valid splice sites in {} ({}) '
                          'splice type ...'.format(splice_name_spaces,
                                                   splice_abbrev.upper()))
            isoform_exons = common.SPLICE_TYPE_ISOFORM_EXONS[splice_abbrev]

            validated_folder = os.path.join(self.index_folder, splice_abbrev,
                                            'validated')
            self.maybe_make_folder(validated_folder)

            splice_sites_seriess = []

            for isoform, exons in isoform_exons.items():
                valid_str = ' or '.join(valid_splice_sites)
                util.progress('\tFinding valid splice sites for {isoform} of'
                              ' {splice_name} events which match '
                              '{valid_splice_sites}'
                              '...'.format(isoform=isoform,
                                           splice_name=splice_name_spaces,
                                           valid_splice_sites=valid_str))
                exon_pairs = zip(exons, exons[1:])
                for exonA, exonB in exon_pairs:
                    util.progress('\t\tFinding splice sites for {exonA} and '
                                  '{exonB} ...'.format(exonA=exonA,
                                                       exonB=exonB))
                    intron_splice_site = self.exon_pair_splice_sites(
                        exonA, exonB, splice_abbrev)
                    splice_sites_seriess.append(intron_splice_site)
                    util.done(4)
                util.done(3)
            splice_sites = pd.concat(splice_sites_seriess, axis=1)

            csv = os.path.join(self.index_folder, splice_abbrev,
                               'splice_sites.csv')
            util.progress('\tWriting splice sites to {csv} ...'.format(
                csv=csv))
            splice_sites.to_csv(csv)
            util.done(3)

            n_total = len(splice_sites.groupby(level=0, axis=0))
            splice_sites_is_valid = splice_sites.isin(valid_splice_sites)
            valid_events_rows = splice_sites_is_valid.all(axis=1)
            splice_sites_validated = splice_sites.loc[valid_events_rows]
            n_valid = len(splice_sites_validated.groupby(level=0, axis=0))

            util.progress("\tValidated {valid}/{total} {splice_name} "
                          "({splice_abbrev}) events. "
                          "".format(valid=n_valid, total=n_total,
                                    splice_name=splice_name_spaces,
                                    splice_abbrev=splice_abbrev.upper()))

            original_events_csv = os.path.join(self.index_folder,
                                               splice_abbrev, EVENTS_CSV)
            validated_events_csv = os.path.join(validated_folder, EVENTS_CSV)
            util.progress('\tWriting validated events to {csv} ...'.format(
                csv=validated_events_csv))

            with open(validated_events_csv, 'w') as f_validated:
                with open(original_events_csv) as f_original:
                    for i, line in enumerate(f_original):
                        if i == 0:
                            f_validated.write(line)
                            continue
                        if line.split(',')[0] in splice_sites_validated.index:
                            f_validated.write(line)
            util.done(3)


class Psi(Subcommand):

    # Instantiate empty variables here so PyCharm doesn't get mad at me
    index = INDEX
    reads_col = None
    sample_id_col = None
    junction_id_col = None

    required_cols = {'--reads-col': reads_col,
                     '--sample-id-col': sample_id_col,
                     '--junction-id-col': junction_id_col}

    def __init__(self, **kwargs):
        # Read all arguments and set as attributes of this class
        for key, value in kwargs.items():
            setattr(self, key, value)

        # If no junction reads csv was specified, use the default one in
        # The outrigger output folder
        if self.junction_reads_csv is None:
            self.junction_reads_csv = JUNCTION_READS_PATH

        if not os.path.exists(self.index_folder):
            raise OSError("The index folder ({}) doesn't exist! Cowardly "
                          "exiting because I don't know what events to "
                          "calcaulate psi on :(".format(self.index_folder))

        for splice_name, splice_folder in self.splice_type_folders.items():
            if not os.path.exists(splice_folder):
                raise OSError(
                    "The splicing index of {} ({}) splice types "
                    "doesn't exist! Cowardly existing because I "
                    "don't know how to define events :(".format(
                        splice_name, splice_folder))

        if not os.path.exists(self.junction_reads_csv):
            raise OSError(
                "The junction reads csv file ({}) doesn't exist! "
                "Cowardly exiting because I don't have the junction "
                "counts calcaulate psi on :(".format(
                    self.junction_reads_csv))

        for folder in self.folders:
            self.maybe_make_folder(folder)

    @property
    def psi_folder(self):
        return os.path.join(self.output_folder, 'psi')

    @property
    def index_folder(self):
        if not hasattr(self, 'index'):
            return os.path.join(self.output_folder, 'index')
        else:
            return self.index

    @property
    def splice_type_folders(self):
        return dict((splice_name, os.path.join(self.index_folder,
                                               splice_abbrev))
                    for splice_name, splice_abbrev in
                    outrigger.common.SPLICE_TYPES)

    @property
    def folders(self):
        return self.output_folder, self.psi_folder

    def maybe_read_junction_reads(self):
        try:
            dtype = {self.reads_col: np.float32}
            if self.junction_reads_csv is None:
                self.junction_reads_csv = JUNCTION_READS_PATH
            util.progress(
                'Reading splice junction reads from {} ...'.format(
                    self.junction_reads_csv))
            junction_reads = pd.read_csv(
                self.junction_reads_csv, dtype=dtype)
            util.done()
        except OSError:
            raise IOError(
                "There is no junction reads file at the expected location"
                " ({csv}). Are you in the correct directory?".format(
                    csv=self.junction_reads_csv))
        return junction_reads

    def validate_junction_reads_data(self, junction_reads):
        for flag, col in self.required_cols.items():
            if col not in junction_reads:
                raise ValueError(
                    "The required column name {col} does not exist in {csv}. "
                    "You can change this with the command line flag, "
                    "{flag}".format(col=col, csv=self.junction_reads_csv,
                                    flag=flag))

    def maybe_get_validated_events(self, splice_abbrev):
        splice_folder = os.path.join(self.index_folder, splice_abbrev)
        events = os.path.join(splice_folder, EVENTS_CSV)
        validated_events = os.path.join(splice_folder, 'validated', EVENTS_CSV)
        if os.path.exists(validated_events):
            return validated_events
        else:
            return events

    def execute(self):
        """Calculate percent spliced in (psi) of splicing events"""

        logger = logging.getLogger('outrigger.psi')
        if self.debug:
            logger.setLevel(10)

        junction_reads = self.maybe_read_junction_reads()

        junction_reads = junction_reads.set_index(
            [self.junction_id_col, self.sample_id_col])
        junction_reads.sort_index(inplace=True)
        logger.debug('\n--- Splice Junction reads ---')
        logger.debug(repr(junction_reads.head()))

        psis = []
        for splice_name, splice_abbrev in outrigger.common.SPLICE_TYPES:
            filename = self.maybe_get_validated_events(splice_abbrev)
            # event_type = os.path.basename(filename).split('.csv')[0]
            util.progress('Reading {name} ({abbrev}) events from {filename}'
                          ' ...'.format(name=splice_name, abbrev=splice_abbrev,
                                        filename=filename))

            event_annotation = pd.read_csv(filename, index_col=0)
            util.done()

            isoform_junctions = outrigger.common.ISOFORM_JUNCTIONS[
                splice_abbrev]
            logger.debug('\n--- Splicing event annotation ---')
            logger.debug(repr(event_annotation.head()))

            util.progress(
                'Calculating percent spliced-in (Psi) scores on '
                '{name} ({abbrev}) events ...'.format(
                    name=splice_name, abbrev=splice_abbrev))
            event_psi = compute.calculate_psi(
                event_annotation, junction_reads,
                min_reads=self.min_reads, debug=self.debug,
                reads_col=self.reads_col, **isoform_junctions)
            csv = os.path.join(self.psi_folder, splice_abbrev,
                               'psi.csv'.format(splice_abbrev))
            self.maybe_make_folder(os.path.dirname(csv))
            event_psi.to_csv(csv)
            psis.append(event_psi)
            util.done()

        util.progress('Concatenating all calculated psi scores '
                      'into one big matrix...')
        splicing = pd.concat(psis, axis=1)
        util.done()
        splicing = splicing.T
        csv = os.path.join(self.psi_folder, 'outrigger_psi.csv')
        util.progress('Writing a samples x features matrix of Psi '
                      'scores to {} ...'.format(csv))
        splicing.to_csv(csv)
        util.done()


def main():
    try:
        cl = CommandLine(sys.argv[1:])
    except Usage:
        cl.do_usage_and_die()


if __name__ == '__main__':
    main()
