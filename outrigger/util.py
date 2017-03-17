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


def extract_alternative_constitutive(psi):
    """Separate psi matrix to events that are alternative vs constitutive

    Parameters
    ----------
    psi : pandas.DataFrame
        This is a (samples, features) shaped dataframe of the percent
        spliced-in values

    Returns
    -------
    alternative, constitutively0, constitutively1 : pandas.DataFrame
        Slices of the input dataframe that are alternative or constitutive
    """

    notnull = psi.notnull()

    constitutively0 = (psi == 0)[notnull].all()
    constitutively1 = (psi == 1)[notnull].all()
    alternative = psi.columns[(~constitutively0) & (~constitutively1)]

    constitutively0 = constitutively0[constitutively0].index
    constitutively1 = constitutively1[constitutively1].index

    return psi[alternative], psi[constitutively0], psi[constitutively1]
