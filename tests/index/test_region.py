import pytest


@pytest.fixture
def location(strand):
    return 'chr1:100-200:{}'.format(strand)


@pytest.fixture
def location_with_region_name():
    return 'junction:chr1:100-200:+'


class TestRegion(object):

    def test___init(self, location):
        from outrigger.index.region import Region

        r = Region(location)

        assert r.region is None
        assert r.chrom == 'chr1'
        assert r.start == 100
        assert r.stop == 200
        assert r.strand == location[-1]
        assert r.name == location

    def test___init_region_name(self, location_with_region_name):
        from outrigger.index.region import Region

        r = Region(location_with_region_name)
        assert r.region == 'junction'
        assert r.chrom == 'chr1'
        assert r.start == 100
        assert r.stop == 200
        assert r.strand == '+'
        assert r.name == location_with_region_name

    @pytest.mark.xfail
    def test___init_start_larger_than_stop(self):
        from outrigger.index.region import Region

        Region('chr1:200-100:+')

    def test__start(self, location):
        from outrigger.index.region import Region

        r = Region(location)

        true__start = -r.start if r.strand == '-' else r.start
        assert r._start == true__start

    def test__stop(self, location):
        from outrigger.index.region import Region

        r = Region(location)

        true__stop = -r.stop if r.strand == '-' else r.stop
        assert r._stop == true__stop

    def test___len(self, location):
        from outrigger.index.region import Region

        r = Region(location)

        assert len(r) == 101

    def test___str(self, location):
        from outrigger.index.region import Region

        r = Region(location)

        true = 'outrigger.Region ({0})'.format(location)
        assert str(r) == true

    def test___str_with_region_name(self, location_with_region_name):
        from outrigger.index.region import Region

        r = Region(location_with_region_name)

        true = 'outrigger.Region ({0})'.format(location_with_region_name)
        assert str(r) == true

    def test___eq(self, location_with_region_name):
        from outrigger.index.region import Region

        r1 = Region(location_with_region_name)
        r2 = Region(location_with_region_name)

        assert r1 == r2

    def test___eq_not_region(self, location_with_region_name):
        from outrigger.index.region import Region

        r1 = Region(location_with_region_name)

        assert not r1 == location_with_region_name

    def test___neq(self, location_with_region_name, location):
        from outrigger.index.region import Region

        r1 = Region(location_with_region_name)
        r2 = Region(location)

        assert r1.__neq__(r2)

    def test_overlaps_true(self, location_with_region_name, location):
        from outrigger.index.region import Region

        r1 = Region(location_with_region_name)
        r2 = Region(location)

        assert r1.overlaps(r2)
        assert r2.overlaps(r1)

    def test_overlaps_false_same_chrom(self, location):
        from outrigger.index.region import Region

        r1 = Region(location)
        r2 = Region('chr1:400-500:-')

        assert not r1.overlaps(r2)
        assert not r2.overlaps(r1)

    def test_overlaps_false_different_chrom(self, location):
        from outrigger.index.region import Region

        r1 = Region(location)
        r2 = Region('chr2:400-500:-')

        assert not r1.overlaps(r2)
        assert not r2.overlaps(r1)
