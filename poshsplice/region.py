"""Define locations in the genome"""

import six

class Region(object):

    __slots__ = ('name', 'region', 'chrom', 'start', 'stop', '_start',
                 '_stop', 'strand')

    def __init__(self, name):
        """A location in the genome

        Parameters
        ----------
        name : str
            A string in the form, region:chrom:start-stop:strand, e.g.
            "exon:chr1:100-200:+" Start must always be smaller than stop.
        """
        self.name = name

        try:
            region, chrom, startstop, strand = name.split(':')
        except ValueError:
            # There is no "region"
            chrom, startstop, strand = name.split(':')
            region = None
        start, stop = map(int, startstop.split('-'))
        if start > stop:
            raise ValueError('Start ({0}) cannot be smaller than stop'
                             ' ({1})'.format(start, stop))

        self.region = region
        self.chrom = six.u(chrom)

        # Use "private" variables to store the true genome location
        self._start = start
        self._stop = stop
        self.strand = six.u(strand)

    # Use non-private start, stop to save negative strand coordinates with
    # negative numbers to make the calculation of "before" and "after" in
    # the transcript easier
    @property
    def start(self):
        if self.strand == '-':
            return -self._start
        else:
            return self._start

    @property
    def stop(self):
        if self.strand == '-':
            return -self._stop
        else:
            return self._stop

    @property
    def name(self):
        base = '{0}:{1}-{2}:{3}'.format(self.chrom, self._start,
                                        self._stop, self.strand)
        if self.region is not None:
            base = self.region + ':' + base
        return six.u(base)

    def __repr__(self):
        return 'poshsplice.Region <{0}>'.format(self.name)

    def __str__(self):
        return six.u(self.name)

    def __eq__(self, other):
        return all(getattr(self, attr) == getattr(other, attr)
                   for attr in self.__slots__)

    def overlaps(self, other):
        """Returns true if any part of other region is contained in this one"""
        if other._start > self._stop or other._stop < self._start:
            return False
        else:
            return True