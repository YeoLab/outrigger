import datetime
import sys


def timestamp():
    return str(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))


def done(n_tabs=2):
    """Write timestamp plus 'Done.' to stdout

    Parameters
    ----------
    n_tabs : int
        Number of tabs to include. Default is 2
    """
    sys.stdout.write('{}{}Done.\n'.format(timestamp(), '\t' * n_tabs))


def progress(message):
    """Write a timestamped progress message to standard output"""
    sys.stdout.write('{}\t{}\n'.format(timestamp(), message))
