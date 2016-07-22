import datetime
import sys


def timestamp():
    return str(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))


def done():
    """Write timestamp plus 'Done.' to stdout"""
    sys.stdout.write('{}\t\tDone.\n'.format(timestamp()))


def progress(message):
    """Write a timestamped progress message to standard output"""
    sys.stdout.write('{}\t{}\n'.format(timestamp(), message))
