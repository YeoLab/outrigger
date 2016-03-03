import datetime
import sys


def timestamp():
    return str(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))


def done():
    """Write timestamp plus 'Done.' to stdout"""
    sys.stdout.write('{}\t\tDone.\n'.format(timestamp()))
