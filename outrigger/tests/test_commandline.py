
import filecmp
import os

import pandas as pd
import pandas.util.testing as pdt
import pytest


class TestSubcommand(object):

    def test___init__(self, tmpdir):
        from outrigger.commandline import Subcommand

        output = os.path.join(tmpdir.strpath, 'blue_ivy')

        kwargs = dict(asdf="beyonce", jkl=1234, output=output)

        subcommand = Subcommand(**kwargs)

        for key, value in kwargs.items():
            assert getattr(subcommand, key) == value
        for folder in subcommand.folders:
            assert os.path.exists(folder)


def assert_directories_equal(dir1, dir2, ignore=None,
                             sortables=('exon', 'junction')):
    """Compare contents of subdirectories to assert they are equal"""

    directory_comparison = filecmp.dircmp(dir1, dir2, ignore=ignore)

    directory_comparison.report_full_closure()
    assert len(directory_comparison.left_only) == 0
    assert len(directory_comparison.right_only) == 0

    for subdir in directory_comparison.subdirs.values():
        print(subdir.common_files)

        assert len(subdir.left_only) == 0
        assert len(subdir.right_only) == 0
        for filename in subdir.common_files:
            filename1 = os.path.join(subdir.left, filename)
            filename2 = os.path.join(subdir.right, filename)

            df1, df2 = None, None

            # If the files are csv or bed tables, check that they're equal
            if filename.endswith('.csv'):
                df1 = pd.read_csv(filename1)
                df2 = pd.read_csv(filename2)

                df1.sort_values(df1.columns.tolist(), inplace=True)
                df2.sort_values(df2.columns.tolist(), inplace=True)
            elif filename.endswith('.bed'):
                df1 = pd.read_table(filename1, header=None)
                df2 = pd.read_table(filename2, header=None)

                df1.sort_values([3, 0, 1, 2], kind='mergesort', inplace=True)
                df2.sort_values([3, 0, 1, 2], kind='mergesort', inplace=True)

            if df1 is not None:
                df1.index = range(len(df1.index))
                df2.index = range(len(df2.index))

                df1.sort_index(axis=1, inplace=True)
                df2.sort_index(axis=1, inplace=True)

                pdt.assert_frame_equal(df1, df2)
                continue

            # Otherwise, just use the file sizes
            size1 = os.stat(filename1).st_size
            size2 = os.stat(filename2).st_size
            try:
                assert size1 == size2
            except AssertionError:
                raise AssertionError('{f1} ({size1}) and {f2} ({size2}) have '
                                 'different sizes'.format(f1=filename1,
                                                          f2=filename2,
                                                          size1=size1,
                                                          size2=size2))


class TestCommandLine(object):

    def test_no_arguments(self, capsys):
        """
        User passes no args, should fail with SystemExit
        """

        from outrigger.commandline import CommandLine

        CommandLine()

        text = '[-h] [--version] {index,validate,psi} ...'
        out, err = capsys.readouterr()

        # Argparse for Python2 sends the version info to stderr, but Python3
        # argparse sends the info to stdout so we concatenate here
        outerr = out + err
        assert 'usage' in outerr
        assert text in outerr

    def test_help(self, capsys):
        """
        User asks for help, should SystemExit and give helpful output
        """
        from outrigger.commandline import CommandLine

        with pytest.raises(SystemExit):
            CommandLine(['--help'])

        out, err = capsys.readouterr()
        assert 'outrigger' in out
        assert 'psi' in out
        assert 'validate' in out
        assert 'usage' in out

    def test_main_version(self, capsys):
        from outrigger.commandline import CommandLine
        from outrigger import __version__

        with pytest.raises(SystemExit):
            CommandLine(['--version'])

        out, err = capsys.readouterr()

        # Argparse for Python2 sends the version info to stderr, but Python3
        # argparse sends the info to stdout so we concatenate here
        outerr = out + err
        assert 'outrigger' in outerr
        assert __version__ in outerr

    def test_main_index(self, tmpdir, tasic2016_unprocessed, sj_filenames,
                        tasic2016_outrigger_output_index):
        from outrigger.commandline import CommandLine

        output_folder = tmpdir.strpath

        gtf = os.path.join(tasic2016_unprocessed, 'gtf',
                           'gencode.vM10.annotation.subset.gtf')
        arguments = ['index', '--sj-out-tab']
        arguments.extend(sj_filenames)
        arguments.extend(['--gtf', gtf, '--output', output_folder,
                          '--n-jobs', '1'])
        # import pdb; pdb.set_trace()
        # assert False
        CommandLine(arguments)

        dir1 = os.path.join(output_folder, 'index')
        dir2 = tasic2016_outrigger_output_index
        ignore = ['psi', '.DS_Store', 'validated', 'splice_sites.csv',
                  # Databases get stored in a weird random way... we're still
                  # checking that the final gtfs are the same
                  'gencode.vM10.annotation.subset.gtf.db']
        assert_directories_equal(dir1, dir2, ignore)

    def test_main_index_reads_csv(self, tmpdir, tasic2016_unprocessed,
                                  tasic2016_outrigger_output,
                                  tasic2016_outrigger_output_index):
        from outrigger.commandline import CommandLine

        output_folder = tmpdir.strpath

        gtf = os.path.join(tasic2016_unprocessed, 'gtf',
                           'gencode.vM10.annotation.subset.gtf')
        arguments = ['index', '--junction-reads-csv']
        arguments.append(os.path.join(tasic2016_outrigger_output, 'junctions',
                                      'reads.csv'))
        arguments.extend(['--gtf', gtf, '--output', output_folder,
                          '--n-jobs', '1'])
        # import pdb; pdb.set_trace()
        # assert False
        CommandLine(arguments)

        dir1 = os.path.join(output_folder, 'index')
        dir2 = tasic2016_outrigger_output_index
        ignore = ['psi', '.DS_Store', 'validated', 'splice_sites.csv',
                  # Databases get stored in a weird random way... we're still
                  # checking that the final gtfs are the same
                  'gencode.vM10.annotation.subset.gtf.db']
        assert_directories_equal(dir1, dir2, ignore)

    def test_main_index_parallelized(self, tmpdir, tasic2016_unprocessed,
                                     sj_filenames,
                                     tasic2016_outrigger_output_index):
        from outrigger.commandline import CommandLine

        output_folder = tmpdir.strpath

        gtf = os.path.join(tasic2016_unprocessed, 'gtf',
                           'gencode.vM10.annotation.subset.gtf')
        arguments = ['index', '--sj-out-tab']
        arguments.extend(sj_filenames)
        arguments.extend(['--gtf', gtf, '--output', output_folder,
                          '--n-jobs', '-1'])
        # import pdb; pdb.set_trace()
        # assert False
        CommandLine(arguments)

        dir1 = os.path.join(output_folder, 'index')
        dir2 = tasic2016_outrigger_output_index
        ignore = ['psi', '.DS_Store', 'validated', 'splice_sites.csv',
                  # Databases get stored in a weird random way... we're still
                  # checking that the final gtfs are the same
                  'gencode.vM10.annotation.subset.gtf.db']
        assert_directories_equal(dir1, dir2, ignore)

    def test_main_index_bam(self, tmpdir, tasic2016_unprocessed,
                            bam_filenames, tasic2016_outrigger_output_bam):
        from outrigger.commandline import CommandLine

        output_folder = tmpdir.strpath

        gtf = os.path.join(tasic2016_unprocessed, 'gtf',
                           'gencode.vM10.annotation.subset.gtf')
        arguments = ['index', '--bam']
        arguments.extend(bam_filenames)
        arguments.extend(['--gtf', gtf, '--output', output_folder])
        # import pdb; pdb.set_trace()
        # assert False
        CommandLine(arguments)

        dir1 = os.path.join(output_folder, 'index')
        dir2 = os.path.join(tasic2016_outrigger_output_bam, 'index')
        ignore = ['psi', '.DS_Store', 'validated', 'splice_sites.csv',
                  # Databases get stored in a weird random way... we're still
                  # checking that the final gtfs are the same
                  'gencode.vM10.annotation.subset.gtf.db']
        assert_directories_equal(dir1, dir2, ignore)

    def test_main_validate(self, tmpdir, negative_control_folder,
                           negative_control_output):
        from outrigger.commandline import CommandLine

        args = ['validate', '--genome',
                '{folder}/chromsizes'.format(
                    folder=negative_control_folder),
                '--fasta',
                '{folder}/genome.fasta'.format(
                    folder=negative_control_folder),
                '--output', tmpdir.strpath,
                '--index', os.path.join(negative_control_output, 'index')]
        CommandLine(args)

        dir1 = tmpdir.strpath
        dir2 = negative_control_output
        assert_directories_equal(dir1, dir2,
                                 ignore=['.DS_Store', 'junctions', 'gtf'])

    def test_main_psi(self, tmpdir, tasic2016_unprocessed,
                      tasic2016_outrigger_output, sj_filenames):
        from outrigger.commandline import CommandLine

        output_folder = tmpdir.strpath

        gtf = os.path.join(tasic2016_unprocessed, 'gtf',
                           'gencode.vM10.annotation.subset.gtf')
        arguments = ['index', '--sj-out-tab']
        arguments.extend(sj_filenames)
        arguments.extend(['--gtf', gtf, '--output', output_folder])
        CommandLine(arguments)

        args = ['psi', '--output', output_folder, '--n-jobs', '1']
        CommandLine(args)

        dir1 = output_folder
        dir2 = tasic2016_outrigger_output
        assert_directories_equal(dir1, dir2, ignore=['.DS_Store'])

    def test_main_psi_parallelized(self, tmpdir, tasic2016_unprocessed,
                                   tasic2016_outrigger_output, sj_filenames):
        from outrigger.commandline import CommandLine

        output_folder = tmpdir.strpath

        gtf = os.path.join(tasic2016_unprocessed, 'gtf',
                           'gencode.vM10.annotation.subset.gtf')
        arguments = ['index', '--sj-out-tab']
        arguments.extend(sj_filenames)
        arguments.extend(['--gtf', gtf, '--output', output_folder])
        CommandLine(arguments)

        args = ['psi', '--output', output_folder, '--n-jobs', '-1']
        CommandLine(args)

        dir1 = output_folder
        dir2 = tasic2016_outrigger_output
        assert_directories_equal(dir1, dir2, ignore=['.DS_Store'])

    def test_main_psi_bam(self, tmpdir, tasic2016_outrigger_output_index,
                          tasic2016_outrigger_output_bam, bam_filenames):
        from outrigger.commandline import CommandLine

        output_folder = tmpdir.strpath

        args = ['psi', '--output', output_folder, '--n-jobs', '1',
                '--index', tasic2016_outrigger_output_index,
                '--bam']
        args.extend(bam_filenames)
        CommandLine(args)

        dir1 = output_folder
        dir2 = tasic2016_outrigger_output_bam
        assert_directories_equal(dir1, dir2, ignore=['.DS_Store', 'index'])
