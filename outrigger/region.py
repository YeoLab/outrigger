"""Define locations in the genome"""


STRANDS = '+', '-', '.'


class Region(object):

    __slots__ = ('region', 'chrom', 'start', 'stop', 'strand')

    def __init__(self, name=None, region=None, chrom=None, start=None,
                 stop=None, strand=None):
        """A location in the genome

        Parameters
        ----------
        name : str
            A string of either of the two forms:
                - chrom:start-stop:strand, e.g. "chr1:100-200:-"
                - region:chrom:start-stop:strand, e.g. "exon:chr1:100-200:+"
            Start must always be smaller than stop.
        """
        if name is not None:
            region = None
            try:
                region, chrom, startstop, strand = name.split(':')
            except ValueError:
                # There is no "region"
                chrom, startstop, strand = name.split(':')
            start, stop = map(int, startstop.split('-'))
            if start > stop:
                raise ValueError('Start ({0}) cannot be larger than stop'
                                 ' ({1})'.format(start, stop))
        self.region = region
        self.chrom = chrom

        self.start = start
        self.stop = stop
        self.strand = strand

    # Use non-private start, stop to save negative strand coordinates with
    # negative numbers to make the calculation of "before" and "after" in
    # the transcript easier
    @property
    def _start(self):
        """Relative start of the feature

        If strand is negative, this is negative so relativity is easy"""
        if self.strand == '-':
            return -self.start
        else:
            return self.start

    @property
    def _stop(self):
        """Relative start of the feature

        If strand is negative, this is negative so relativity is easy"""
        if self.strand == '-':
            return -self.stop
        else:
            return self.stop

    @property
    def name(self):
        base = '{0}:{1}-{2}:{3}'.format(self.chrom, self.start,
                                        self.stop, self.strand)
        if self.region is not None:
            base = self.region + ':' + base
        return base

    def __len__(self):
        """Length of region. Add 1 to include last base of stop"""
        return self.stop - self.start + 1

    def __repr__(self):
        return 'outrigger.Region ({0})'.format(self.name)

    def __str__(self):
        return self.__repr__()

    def __eq__(self, other):
        if isinstance(other, Region):
            return all(getattr(self, attr) == getattr(other, attr)
                       for attr in self.__slots__)
        else:
            return False

    def __neq__(self, other):
        return not self.__eq__(other)

    def overlaps(self, other):
        """Returns true if any part of other region is contained in this one"""
        if other.chrom == self.chrom:
            if other.start > self.stop or other.stop < self.start:
                return False
            else:
                return True
        else:
            return False

    def to_zero_based(self):
        """Convert genome coordinates to 0-based

        Assumes that this event is one-based.
        """
        region = self.region
        chrom = self.chrom
        start = self.start - 1
        stop = self.stop
        strand = self.strand

        base = '{chrom}:{start}-{stop}:{strand}'.format(
            chrom=chrom, start=start, stop=stop, strand=strand)

        if region is not None:
            base = '{region}:{base}'.format(region=region, base=base)
        return Region(base)

    def to_bed_format(self, name=None):
        name = self.__repr__() if name is None else name

        s = '{chrom}\t{start}\t{stop}\t{name}\t{score}\t{strand}'.format(
            chrom=self.chrom, start=self.start-1, stop=self.stop, name=name,
            score='.', strand=self.strand)
        return s
