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


def strip_strand(location):
    """Remove strand information from a string location

    If not a string, return the original thing

    >>> strip_strand('junction:chr1:1309826-1310084:-')
    'junction:chr1:1309826-1310084'
    """
    try:
        return location.replace(':+', '').replace(':-', '').replace(
            ':undefined', '')
    except:
        return location
