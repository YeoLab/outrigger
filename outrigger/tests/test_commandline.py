
import filecmp
import glob
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


class TestCommandLine(object):

    def test_no_arguments(self):
        """
        User passes no args, should fail with SystemExit
        """

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
        assert 'outrigger' in out
        assert __version__ in out


    def test_main_index(self, tmpdir, capsys, tasic2016_unprocessed):
        from outrigger.commandline import CommandLine

        sj_out_tab_globber = os.path.join(tasic2016_unprocessed, 'sj_out_tab',
                                          '*SJ.out.tab')
        gtf = os.path.join(tasic2016_unprocessed, 'gtf',
                           'gencode.vM10.annotation.subset.gtf')
        arguments = ['index', '--sj-out-tab']
        arguments.extend(glob.iglob(sj_out_tab_globber))
        arguments.extend(['--gtf', gtf, '--output', tmpdir.strpath, '--debug'])
        # import pdb; pdb.set_trace()
        # assert False
        CommandLine(arguments)

        out, err = capsys.readouterr()

        dir1 = os.path.join(tmpdir.strpath, 'index')
        dir2 = os.path.join('outrigger', 'tests', 'data', 'tasic2016',
                            'outrigger_output', 'index')
        ignore = ['psi', '.DS_Store', 'validated', 'splice_sites.csv',
                  # Databases get stored in a weird random way... we're still
                  # checking that the final gtfs are the same
                  'gencode.vM10.annotation.subset.gtf.db']
        directory_comparison = filecmp.dircmp(dir1, dir2, ignore=ignore)

        directory_comparison.report_full_closure()
        assert len(directory_comparison.left_only) == 0
        assert len(directory_comparison.right_only) == 0
        # assert False

        for subdir in directory_comparison.subdirs.values():
            print(subdir.common_files)

            assert len(subdir.left_only) == 0
            assert len(subdir.right_only) == 0
            for filename in subdir.common_files:
                filename1 = os.path.join(subdir.left, filename)
                filename2 = os.path.join(subdir.right, filename)

                # If the files are csv or bed tables, check that they're equal
                df1, df2 = None, None

                if filename.endswith('.csv'):
                    df1 = pd.read_csv(filename1)
                    df2 = pd.read_csv(filename2)

                    df1.sort_index(inplace=True)
                    df2.sort_index(inplace=True)
                elif filename.endswith('.bed'):
                    df1 = pd.read_table(filename1, header=None)
                    df2 = pd.read_table(filename2, header=None)

                    df1.sort_values([3, 1, 2], kind='mergesort', inplace=True)
                    df2.sort_values([3, 1, 2], kind='mergesort', inplace=True)

                    df1.index = range(len(df1.index))
                    df2.index = range(len(df2.index))

                if df1 is not None:
                    pdt.assert_frame_equal(df1, df2)

                # Otherwise, just use the file sizes
                stat1 = os.stat(filename1)
                stat2 = os.stat(filename2)
                assert stat1.st_size == stat2.st_size

    def test_main_psi(self, tmpdir, tasic_unprocessed):
        from outrigger.commandline import CommandLine

        args = 'index --sj-out-tab {folder}/sj_out_tab/* ' \
               '--gtf {folder}/gtf/gencode.vM10.annotation.subset.gtf' \
               ''.format(tasic_unprocessed).split()
        CommandLine(args)

        args = ['psi']
        CommandLine(args)

        dir1 = tmpdir.strpath
        dir2 = os.path.join('outrigger', 'test_data', 'tasic2016',
                            'outrigger_output')
        directory_comparison = filecmp.dircmp(dir1, dir2,
                                              ignore=['.DS_Store'])
        assert len(directory_comparison.left_only) == 0
        assert len(directory_comparison.right_only) == 0


    def test_main_validate(self, tmpdir, negative_control_folder):
        from outrigger.commandline import CommandLine

        args =['validate', '--genome',
               '{folder}/chromsizes'.format(folder=negative_control_folder),
               '--fasta',
               '{folder}/genome.fasta'.format(folder=negative_control_folder)]
        CommandLine(args)

        dir1 = tmpdir.strpath
        dir2 = os.path.join('outrigger', 'test_data', 'tasic2016',
                            'outrigger_output')
        directory_comparison = filecmp.dircmp(dir1, dir2,
                                              ignore=['.DS_Store'])
        assert len(directory_comparison.left_only) == 0
        assert len(directory_comparison.right_only) == 0
