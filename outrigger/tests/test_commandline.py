
# import filecmp
# import glob
import os
import subprocess

import pytest


slow = pytest.mark.skipif(
    not pytest.config.getoption("--runslow"),
    reason="need --run-slow option to run"
)


class TestSubcommand(object):

    def test___init(self):
        from outrigger.commandline import Subcommand

        kwargs = dict(asdf="beyonce", jkl=1234)

        s = Subcommand(**kwargs)

        for key, value in kwargs.items():
            assert getattr(s, key) == value


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

@slow
def test_make_arabdopsis(outrigger_folder):
    os.chdir(outrigger_folder)

    command = 'make arabdopsis'
    args = command.split()

    outrigger_output = str(subprocess.check_output(args))
    assert 'Found 14 SE events' in outrigger_output
    assert '11371/21995 junctions remain' in outrigger_output
    assert '25 novel exons on chromosome 4' in outrigger_output
    assert 'No MXE events found' in outrigger_output

# def test_main_index(tmpdir):
#     from outrigger.commandline import CommandLine
#
#     sj_out_tab = ' '.join(glob.iglob('outrigger/test_data/tasic2016/unprocessed/sj_out_tab/*'))  # noqa
#     gtf = ' outrigger/test_data/tasic2016/unprocessed/gtf/gencode.vM10.annotation.snap25.myl6.gtf'  # noqa
#     args = 'index --sj-out-tab {splice_junctions} --gtf {gtf} --output {output}'.format(  # noqa
#         splice_junctions=sj_out_tab, gtf=gtf, output=tmpdir.strpath).split()
#     # import pdb; pdb.set_trace()
#     CommandLine(args)
#
#     dir1 = os.path.join(tmpdir.strpath, 'index')
#     dir2 = os.path.join('outrigger', 'test_data', 'tasic2016',
#                         'outrigger_output', 'index')
#     directory_comparison = filecmp.dircmp(dir1, dir2,
#                                           ignore=['psi', '.DS_Store'])
#     assert len(directory_comparison.left_only) == 0
#     assert len(directory_comparison.right_only) == 0
#
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
