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


def progress(message, n_tabs=1):
    """Write a timestamped progress message to standard output

    Parameters
    ----------
    n_tabs : int
        Integer number of tabs to pad before the message. Default is 1
    """
    sys.stdout.write('{time}{tabs}{message}\n'.format(time=timestamp(),
                                                      message=message,
                                                      tabs='\t'*n_tabs))
