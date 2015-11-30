from graphlite import connect, V
import numpy as np
import pandas as pd
import pandas.util.testing as pdt
import pytest


@pytest.fixture(params=['positive', 'negative'])
def strand(request):
    if request.param == 'positive':
        return '+'
    else:
        return '-'


@pytest.fixture()
def chrom():
    return 'chr1'


@pytest.fixture
def exon_start_stop():
    return {'exon1alt': (100, 125), 'exon1': (150, 175),
            'exon2a3ss': (200, 250), 'exon2': (225, 250),
            'exon2a5ss': (225, 275),
            'exon3': (300, 350),
            'exon4': (400, 425), 'exon4alt': (475, 500)}


@pytest.fixture
def transcripts():
    return (
        ('Transcript 1', ('exon1', 'exon2', 'exon3', 'exon4')),

        # Alt 1st exon, relative to transcript1
        ('Transcript 2', ('exon1alt', 'exon2', 'exon3', 'exon4')),

        # skipped exon, relative to transcript1
        ('Transcript 3', ('exon1', 'exon3', 'exon4')),

        # Alt 3' splice site, relative to transcript1
        ('Transcript 4', ('exon1', 'exon2a3ss', 'exon3', 'exon4')),

        # Alt 5' splice site, relative to transcript1
        ('Transcript 5', ('exon1', 'exon2a5ss', 'exon3', 'exon4')),

        # MXE, relative to transcript1
        ('Transcript 6', ('exon1', 'exon2', 'exon4')),

        # Twin Cassette, relative to transcript1
        ('Transcript 7', ('exon1', 'exon4')),

        # Alt last exon, relative to transcript1
        ('Transcript 8', ('exon1', 'exon2', 'exon3', 'exon4alt'))
    )


@pytest.fixture(params=[None, 'exon'])
def region(request):
    return request.param


@pytest.fixture
def junction_to_exons(chrom, exon_start_stop, transcripts, strand):
    from collections import defaultdict
    from outrigger.junctions_to_events import stringify_location

    data = defaultdict(lambda: {'upstream': set([]), 'downstream': set([])})

    for transcript, exons in transcripts:
        for exon1, exon2 in zip(exons, exons[1:]):

            start1, stop1 = exon_start_stop[exon1]
            start2, stop2 = exon_start_stop[exon2]
            exon1_location = stringify_location(chrom, start1, stop1,
                                                strand, 'exon')
            exon2_location = stringify_location(chrom, start2, stop2,
                                                strand, 'exon')

            # if strand == '-':
            #     start = stop2 + 1
            #     stop = start1 - 1
            # else:
            start = stop1 + 1
            stop = start2 - 1

            junction_location = stringify_location(chrom, start,
                                                   stop, strand, 'junction')

            if strand == '-':
                data[junction_location]['downstream'].add(exon1_location)
                data[junction_location]['upstream'].add(exon2_location)
            else:
                data[junction_location]['upstream'].add(exon1_location)
                data[junction_location]['downstream'].add(exon2_location)
    data = pd.DataFrame(data).T
    data = data.applymap(lambda x: ','.join(x))
    data = data.reset_index()
    data = data.rename(columns={'index': 'junction'})
    return data


@pytest.fixture
def junction_exon_triples(chrom, exon_start_stop, transcripts, strand):
    from outrigger.junctions_to_events import stringify_location
    data = []

    for transcript, exons in transcripts:
        for exon1, exon2 in zip(exons, exons[1:]):

            start1, stop1 = exon_start_stop[exon1]
            start2, stop2 = exon_start_stop[exon2]
            exon1_location = stringify_location(chrom, start1, stop1,
                                                strand, 'exon')
            exon2_location = stringify_location(chrom, start2, stop2,
                                                strand, 'exon')

            # if strand == '-':
            #     start = stop2 + 1
            #     stop = start1 - 1
            # else:
            start = stop1 + 1
            stop = start2 - 1

            junction_location = stringify_location(chrom, start, stop,
                                                   strand, 'junction')

            if strand == '-':
                data.append(
                    [exon1_location, 'downstream', junction_location])
                data.append(
                    [exon2_location, 'upstream', junction_location])
            else:
                data.append(
                    [exon1_location, 'upstream', junction_location])
                data.append(
                    [exon2_location, 'downstream', junction_location])
    data = pd.DataFrame(data, columns=['exon', 'direction', 'junction'])
    data = data.drop_duplicates()
    return data


def test_stringify_location(chrom, strand, region):
    from outrigger.junctions_to_events import stringify_location

    test = stringify_location(chrom, 100, 200, strand, region)

    if region is None:
        true = '{0}:{1}-{2}:{3}'.format(chrom, 100, 200, strand)
    else:
        true = '{0}:{1}:{2}-{3}:{4}'.format(region, chrom, 100, 200, strand)
    assert test == true


def assert_graph_items_equal(graph1, int_to_item1, item_to_int1, graph2,
                             int_to_item2, item_to_int2):
    from outrigger.junctions_to_events import DIRECTIONS

    for item1, number1 in item_to_int1.iteritems():
        for direction in DIRECTIONS:
            test = int_to_item1.loc[list(
                graph1.find(getattr(V(number1), direction)))].values

            number2 = item_to_int2.loc[item1]
            true = int_to_item2.loc[
                list(graph2.find(getattr(V(number2), direction)))].values

            test.sort()
            true.sort()

            pdt.assert_numpy_array_equal(test, true)

    for item2, number2 in item_to_int2.iteritems():
        for direction in DIRECTIONS:
            test = int_to_item2.loc[list(
                graph2.find(getattr(V(number2), direction)))].values

            number1 = item_to_int1.loc[item2]
            true = int_to_item1.loc[
                list(graph1.find(getattr(V(number1), direction)))].values

            test.sort()
            true.sort()

            pdt.assert_numpy_array_equal(test, true)


class TestAggregateJunctions(object):

    @pytest.fixture
    def junction_aggregator(self, junction_exon_triples):
        from outrigger.junctions_to_events import JunctionAggregator
        return JunctionAggregator(junction_exon_triples)

    def test_init(self, junction_exon_triples, graph):
        from outrigger.junctions_to_events import JunctionAggregator

        true_graph, true_int_to_item, true_item_to_int = graph

        test = JunctionAggregator(junction_exon_triples)
        pdt.assert_frame_equal(test.junction_exon_triples,
                               junction_exon_triples)
        assert test.db is None
        exons = junction_exon_triples.exon.unique()
        junctions = junction_exon_triples.junction.unique()
        items = np.concatenate([exons, junctions])
        int_to_item = pd.Series(items)
        item_to_int = pd.Series(dict((v, k) for k, v in
                                     int_to_item.iteritems()))

        pdt.assert_numpy_array_equal(test.exons, exons)
        pdt.assert_numpy_array_equal(test.junctions, junctions)
        pdt.assert_numpy_array_equal(test.items, items)
        pdt.assert_dict_equal(test.int_to_item, int_to_item)
        pdt.assert_dict_equal(test.item_to_int, item_to_int)

        assert_graph_items_equal(test.graph, test.int_to_item,
                                 test.item_to_int, true_graph,
                                 true_int_to_item,
                                 true_item_to_int)

    def test_from_junction_to_exons(self, junction_to_exons,
                                    junction_aggregator):
        from outrigger.junctions_to_events import JunctionAggregator

        test = JunctionAggregator.from_junction_to_exons(junction_to_exons)

        assert_graph_items_equal(test.graph, test.int_to_item,
                                 test.item_to_int, junction_aggregator.graph,
                                 junction_aggregator.int_to_item,
                                 junction_aggregator.item_to_int)

    def test_skipped_exon(self, junction_aggregator, strand):
        test = junction_aggregator.skipped_exon()

        if strand == '+':
            true = {
                ('exon:chr1:150-175:+',   # Exon 1
                 'exon:chr1:225-250:+',   # Exon 2
                 'exon:chr1:300-350:+'):  # Exon 3
                ('junction:chr1:176-224:+',
                 'junction:chr1:251-299:+',
                 'junction:chr1:176-299:+'),
                ('exon:chr1:150-175:+',   # Exon 1
                 'exon:chr1:225-275:+',   # Exon 2, alt 5' splice site
                 'exon:chr1:300-350:+'):  # Exon 3
                ('junction:chr1:176-224:+',
                 'junction:chr1:276-299:+',
                 'junction:chr1:176-299:+'),
                ('exon:chr1:150-175:+',   # Exon 1
                 'exon:chr1:300-350:+',   # Exon 3
                 'exon:chr1:400-425:+'):  # Exon 4
                ('junction:chr1:176-299:+',
                 'junction:chr1:351-399:+',
                 'junction:chr1:176-399:+'),
                ('exon:chr1:150-175:+',   # Exon 1
                 'exon:chr1:225-250:+',   # Exon 2
                 'exon:chr1:400-425:+'):  # Exon 4
                ('junction:chr1:176-224:+',
                 'junction:chr1:251-399:+',
                 'junction:chr1:176-399:+'),
                ('exon:chr1:225-250:+',   # Exon 2
                 'exon:chr1:300-350:+',   # Exon 3
                 'exon:chr1:400-425:+'):  # Exon 4
                ('junction:chr1:251-299:+',
                 'junction:chr1:351-399:+',
                 'junction:chr1:251-399:+'),
                ('exon:chr1:150-175:+',   # Exon 1
                 'exon:chr1:200-250:+',   # Exon 2, alt 3' splice site
                 'exon:chr1:300-350:+'):  # Exon 4
                ('junction:chr1:176-199:+',
                 'junction:chr1:251-299:+',
                 'junction:chr1:176-299:+')}
        else:
            true = {
                ('exon:chr1:300-350:-',
                 'exon:chr1:225-250:-',
                 'exon:chr1:150-175:-'):
                ('junction:chr1:251-299:-',
                 'junction:chr1:176-224:-',
                 'junction:chr1:176-299:-'),
                ('exon:chr1:300-350:-',
                 'exon:chr1:225-275:-',
                 'exon:chr1:150-175:-'):
                ('junction:chr1:276-299:-',
                 'junction:chr1:176-224:-',
                 'junction:chr1:176-299:-'),
                ('exon:chr1:400-425:-',
                 'exon:chr1:300-350:-',
                 'exon:chr1:150-175:-'):
                ('junction:chr1:351-399:-',
                 'junction:chr1:176-299:-',
                 'junction:chr1:176-399:-'),
                ('exon:chr1:400-425:-',
                 'exon:chr1:300-350:-',
                 'exon:chr1:225-250:-'):
                ('junction:chr1:351-399:-',
                 'junction:chr1:251-299:-',
                 'junction:chr1:251-399:-'),
                ('exon:chr1:300-350:-',
                 'exon:chr1:200-250:-',
                 'exon:chr1:150-175:-'):
                ('junction:chr1:251-299:-',
                 'junction:chr1:176-199:-',
                 'junction:chr1:176-299:-'),
                ('exon:chr1:400-425:-',
                 'exon:chr1:225-250:-',
                 'exon:chr1:150-175:-'):
                ('junction:chr1:251-399:-',
                 'junction:chr1:176-224:-',
                 'junction:chr1:176-399:-')}
        pdt.assert_dict_equal(test, true)

    def test_mutually_exclusive_exon(self, junction_aggregator, strand):
        test = junction_aggregator.mutually_exclusive_exon()

        if strand == '+':
            true = {
                ('exon:chr1:150-175:+',   # Exon 1
                 'exon:chr1:225-250:+',   # Exon 2
                 'exon:chr1:300-350:+',   # Exon 3
                 'exon:chr1:400-425:+'):  # Exon 4
                ('junction:chr1:176-299:+',   # Exon1-Exon2 junction
                 'junction:chr1:351-399:+',   # Exon1-Exon3 junction
                 'junction:chr1:176-224:+',   # Exon2-Exon3 junction
                 'junction:chr1:251-399:+')}  # Exon3-Exon4 junction
        else:
            true = {
                ('exon:chr1:400-425:-',
                 'exon:chr1:300-350:-',
                 'exon:chr1:225-250:-',
                 'exon:chr1:150-175:-'):
                ('junction:chr1:251-399:-',
                 'junction:chr1:176-224:-',
                 'junction:chr1:351-399:-',
                 'junction:chr1:176-299:-')}
        pdt.assert_dict_equal(test, true)

    def test_twin_cassette(self, junction_aggregator, strand):
        # test = junction_aggregator.twin_cassette()

        if strand == '+':
            true = {('exon:chr1:150-175:+',  # Exon 1
                     'exon:chr1:225-250:+',  # Exon 2
                     'exon:chr1:300-350:+',  # Exon 3
                     'exon:chr1:400-425:+'):  # Exon 4
                    ('junction:chr1:176-299:+',  # Exon1-Exon2 junction
                     'junction:chr1:251-299:+',  # Exon2-Exon3 junction
                     'junction:chr1:351-399:+',  # Exon3-Exon4 junction
                     'junction:chr1:176-399:+')  # Exon1-Exon4 junction
                    }
        else:
            true = {('exon:chr1:400-425:-',  # Exon 1
                     'exon:chr1:300-350:-',  # Exon 2
                     'exon:chr1:225-250:-',  # Exon 3
                     'exon:chr1:150-175:-'):  # Exon 4
                    ('junction:chr1:351-399:-',  # Exon1-Exon2 junction
                     'junction:chr1:251-299:-',  # Exon2-Exon3 junction
                     'junction:chr1:176-225:-',  # Exon3-Exon4 junction
                     'junction:chr1:176-399:-')  # Exon1-Exon4 junction
                    }
        return true
        # pdt.assert_dict_equal(test, true)

    def test_a5ss(self, junction_aggregator, strand):
        true = {('exon:chr1:225-250:+',   # Exon 2
                 'exon:chr1:225-275:+',   # Exon 2, Alt 5' splice site
                 'exon:chr1:300-350:+'):  # Exon 3
                ('junction:chr1:251-299:+',    # Exon2-Exon3 junction
                 'junction:chr1:276-299:+')}   # Exon2a5ss-Exon3 junction
        return true

    def test_a3ss(self, junction_aggregator, strand):
        true = {('exon:chr1:150-175:+',   # Exon 1
                 'exon:chr1:200-250:+',   # Exon 2, Alt 3' splice site
                 'exon:chr1:225-250:+'):  # Exon 2
                ('junction:chr1:176-199:+',   # Exon1-Exon2a3ss junction
                 'junction:chr1:176-224:+')}  # Exon1-Exon2 junction
        return true

    def test_afe(self, junction_aggregator, strand):
        true = {('exon:chr1:100-125:+',   # Exon 1 alt
                 'exon:chr1:150-175:+',   # Exon 1
                 'exon:chr1:225-250:+'):  # Exon 2
                ('junction:chr1:126-224:+',   # Exon1alt-Exon2 junction
                 'junction:chr1:176-224:+')}  # Exon1-Exon2 junction
        return true

    def test_ale(self, junction_aggregator, strand):
        true = {('exon:chr1:300-350:+',   # Exon 3
                 'exon:chr1:400-425:+',   # Exon 4
                 'exon:chr1:475-500:+'):  # Exon 4 alt
                ('junction:chr1:351-399:+',   # Exon3-Exon4 junction
                 'junction:chr1:351-474:+')}  # Exon3-Exon4alt junction
        return true


@pytest.fixture
def graph(exon_start_stop, transcripts, chrom, strand):
    from outrigger.junctions_to_events import stringify_location, opposite

    graph = connect(":memory:", graphs=['upstream', 'downstream'])

    items = []
    triples = set()

    for transcript, exons in transcripts:
        for exon1, exon2 in zip(exons, exons[1:]):

            start1, stop1 = exon_start_stop[exon1]
            start2, stop2 = exon_start_stop[exon2]
            exon1_location = stringify_location(chrom, start1, stop1, strand,
                                                'exon')
            exon2_location = stringify_location(chrom, start2, stop2, strand,
                                                'exon')

            # if strand == '-':
            #     start = stop2 + 1
            #     stop = start1 - 1
            # else:
            start = stop1 + 1
            stop = start2 - 1

            junction_location = stringify_location(chrom, start, stop, strand,
                                                   'junction')

            if exon1_location not in items:
                items.append(exon1_location)
            if exon2_location not in items:
                items.append(exon2_location)
            if junction_location not in items:
                items.append(junction_location)

            # Get unique integer for junction
            junction_i = items.index(junction_location)

            if strand == '-':
                exon1_triple = exon1_location, 'downstream', junction_location
                exon2_triple = exon2_location, 'upstream', junction_location
            else:
                exon1_triple = exon1_location, 'upstream', junction_location
                exon2_triple = exon2_location, 'downstream', junction_location

            exon_triples = exon1_triple, exon2_triple

            with graph.transaction() as tr:
                for exon_triple in exon_triples:
                    if exon_triple not in triples:
                        triples.add(exon_triple)

                        exon, direction, junction = exon_triple

                        # Get unique integer for exon
                        exon_i = items.index(exon)
                        tr.store(getattr(V(exon_i), direction)(junction_i))
                        tr.store(getattr(V(junction_i), opposite(direction))(
                            exon_i))
                    else:
                        continue
    int_to_item = pd.Series(items)
    item_to_int = pd.Series(dict((v, k) for k, v in int_to_item.iteritems()))
    return graph, int_to_item, item_to_int
