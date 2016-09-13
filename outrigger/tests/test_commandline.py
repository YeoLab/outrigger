
import subprocess
import filecmp
import glob
import os


class TestSubcommand(object):

    def test___init(self, tmpdir):
        from outrigger.commandline import Subcommand

        output = os.path.join(tmpdir.strpath, 'blue_ivy')

        kwargs = dict(asdf="beyonce", jkl=1234, output=output)

        subcommand = Subcommand(**kwargs)

        for key, value in kwargs.items():
            assert getattr(subcommand, key) == value
        for folder in subcommand.folders:
            assert os.path.exists(folder)


def test_main_help_from_commandline(tmpdir):
    os.chdir(tmpdir.strpath)

    command = 'outrigger -h'
    args = command.split()

    outrigger_output = str(subprocess.check_output(args))
    assert 'outrigger' in outrigger_output
    assert 'psi' in outrigger_output
    assert 'validate' in outrigger_output
    assert 'help' in outrigger_output
    assert 'usage' in outrigger_output


def test_main_index(tmpdir):
    from outrigger.commandline import CommandLine

    sj_out_tab = ' '.join(glob.iglob('outrigger/tests/data/tasic2016/unprocessed/sj_out_tab/*'))  # noqa
    gtf = ' outrigger/tests/data/tasic2016/unprocessed/gtf/gencode.vM10.annotation.subset.gtf'  # noqa
    args = 'index --sj-out-tab {splice_junctions} --gtf {gtf} --output {output}'.format(  # noqa
        splice_junctions=sj_out_tab, gtf=gtf, output=tmpdir.strpath).split()
    # import pdb; pdb.set_trace()
    CommandLine(args)

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

            stat1 = os.stat(filename1)
            stat2 = os.stat(filename2)
            file1 =  open(filename1)
            file2 =  open(filename2)

            assert stat1.st_size == stat2.st_size
            # for line1, line2 in zip(file1, file2):
            #     assert line1 == line2

#
# def test_main_psi(tmpdir):
#     from outrigger.commandline import CommandLine
#
#     args = 'index --sj-out-tab outrigger/test_data/tasic2016/unprocessed/sj_out_tab/* --gtf outrigger/test_data/tasic2016/unprocessed/gtf/gencode.vM10.annotation.snap25.myl6.gtf'.split()  # noqa
#     CommandLine(args)
#
#     args = ['psi']
#     CommandLine(args)
#
#     dir1 = tmpdir.strpath
#     dir2 = os.path.join('outrigger', 'test_data', 'tasic2016',
#                         'outrigger_output')
#     directory_comparison = filecmp.dircmp(dir1, dir2,
#                                           ignore=['.DS_Store'])
#     assert len(directory_comparison.left_only) == 0
#     assert len(directory_comparison.right_only) == 0

#
# def test_main_validate(tmpdir, negative_control_folder):
#     from outrigger.commandline import CommandLine
#
#     args =['validate', '--genome',
#            '{folder}/chromsizes'.format(folder=negative_control_folder),
#            '--fasta',
#            '{folder}/genome.fasta'.format(folder=negative_control_folder)]
#     CommandLine(args)
#
#     dir1 = tmpdir.strpath
#     dir2 = os.path.join('outrigger', 'test_data', 'tasic2016',
#                         'outrigger_output')
#     directory_comparison = filecmp.dircmp(dir1, dir2,
#                                           ignore=['.DS_Store'])
#     assert len(directory_comparison.left_only) == 0
#     assert len(directory_comparison.right_only) == 0
