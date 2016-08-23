
import filecmp
import glob
import os


class TestSubcommand(object):

    def test___init(self):
        from outrigger.commandline import Subcommand

        kwargs = dict(asdf="beyonce", jkl=1234)

        s = Subcommand(**kwargs)

        for key, value in kwargs.items():
            assert getattr(s, key) == value


# def test_main_index(tmpdir):
#     from outrigger.commandline import CommandLine
#
#     sj_out_tab = ' '.join(glob.iglob('outrigger/test_data/tasic2016/unprocessed/sj_out_tab/*'))
#     gtf = ' outrigger/test_data/tasic2016/unprocessed/gtf/gencode.vM10.annotation.snap25.myl6.gtf'
#     args = 'index --sj-out-tab {splice_junctions} --gtf {gtf} --output {output}'.format(
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
#     args = 'index --sj-out-tab outrigger/test_data/tasic2016/unprocessed/sj_out_tab/* --gtf outrigger/test_data/tasic2016/unprocessed/gtf/gencode.vM10.annotation.snap25.myl6.gtf'.split()
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
