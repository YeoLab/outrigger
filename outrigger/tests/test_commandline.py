import os
import sys


class TestCommandLine(object):
    def test___init(self):
        pass


class TestSubcommand(object):

    def test___init(self):
        from outrigger.commandline import Subcommand

        kwargs = dict(asdf="beyonce", jkl=1234)

        s = Subcommand(**kwargs)

        for key, value in kwargs.items():
            assert getattr(s, key) == value

    def test_csv(self):
        pass
