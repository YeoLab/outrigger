import logging
import os
import re

import pandas as pd
import pandas.util.testing as pdt
import pytest
from graphlite import connect, V

logging.basicConfig(level=logging.DEBUG)


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
    from outrigger.index.events import stringify_location

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
    from outrigger.index.events import stringify_location
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
    from outrigger.index.events import stringify_location

    test = stringify_location(chrom, 100, 200, strand, region)

    if region is None:
        true = '{0}:{1}-{2}:{3}'.format(chrom, 100, 200, strand)
    else:
        true = '{0}:{1}:{2}-{3}:{4}'.format(region, chrom, 100, 200, strand)
    assert test == true


def assert_graph_items_equal(graph1, items1, graph2, items2):
    """Checks all relationships in graph1 exist in graph2, and vice versa"""
    from outrigger.common import DIRECTIONS

    for number1, item1 in enumerate(items1):
        for direction in DIRECTIONS:
            test = [items1[i] for i in
                    graph1.find(getattr(V(number1), direction))]

            number2 = items2.index(item1)
            true = [items2[i] for i in
                    graph2.find(getattr(V(number2), direction))]

            test.sort()
            true.sort()

            pdt.assert_equal(test, true)

    for number2, item2 in enumerate(items2):
        for direction in DIRECTIONS:
            test = [items2[i] for i in
                    graph2.find(getattr(V(number2), direction))]

            number1 = items1.index(item2)
            true = [items1[i] for i in
                    graph1.find(getattr(V(number1), direction))]

            test.sort()
            true.sort()

            pdt.assert_equal(test, true)


class TestEventMaker(object):
    @pytest.fixture
    def event_maker(self, junction_exon_triples):
        from outrigger.index.events import EventMaker
        return EventMaker(junction_exon_triples)

    def test_init(self, junction_exon_triples):
        from outrigger.index.events import EventMaker, CHROM

        junction_exon_triples_chrom = junction_exon_triples.copy()
        junction_exon_triples_chrom[CHROM] = \
            junction_exon_triples_chrom['junction'].str.split(':').str[1]

        test = EventMaker(junction_exon_triples)
        pdt.assert_frame_equal(test.junction_exon_triples,
                               junction_exon_triples_chrom)
        assert test.db is None
        assert test.junction_col == 'junction'
        assert test.exon_col == 'exon'

    @pytest.fixture
    def strand_name(self, strand):
        if strand == '+':
            return "positive"
        else:
            return "negative"

    @pytest.fixture
    def events_csv(self, simulated_outrigger_index, strand_name):
        return os.path.join(simulated_outrigger_index, strand_name,
                            'events.csv')

    def test_finding_events(self, event_maker, capsys, strand_name,
                            splice_type, simulated_outrigger_index):
        """Test finding SE and MXE events in one function"""
        from outrigger.common import SPLICE_ABBREVS

        test_events = event_maker.find_events()

        out, err = capsys.readouterr()
        assert 'Combining all events into large dataframes' in out
        assert 'Done.' in out

        for splice_abbrev in SPLICE_ABBREVS:
            test = test_events[splice_abbrev]

            events_csv = os.path.join(
                simulated_outrigger_index, splice_abbrev,
                'events_{}_strand.csv'.format(strand_name))
            true = pd.read_csv(events_csv, index_col=0)

            sort_by = [x for x in true.columns if re.match('exon\d', x)]
            test.sort_values(by=sort_by, inplace=True)
            true.sort_values(by=sort_by, inplace=True)

            pdt.assert_frame_equal(test, true)

    def test_a5ss(self, event_maker, strand):
        true = {('exon:chr1:225-250:+',  # Exon 2
                 'exon:chr1:225-275:+',  # Exon 2, Alt 5' splice site
                 'exon:chr1:300-350:+'):  # Exon 3
                ('junction:chr1:251-299:+',  # Exon2-Exon3 junction
                 'junction:chr1:276-299:+')}  # Exon2a5ss-Exon3 junction
        return true

    def test_a3ss(self, event_maker, strand):
        true = {('exon:chr1:150-175:+',  # Exon 1
                 'exon:chr1:200-250:+',  # Exon 2, Alt 3' splice site
                 'exon:chr1:225-250:+'):  # Exon 2
                ('junction:chr1:176-199:+',  # Exon1-Exon2a3ss junction
                 'junction:chr1:176-224:+')}  # Exon1-Exon2 junction
        return true

    def test_afe(self, event_maker, strand):
        true = {('exon:chr1:100-125:+',  # Exon 1 alt
                 'exon:chr1:150-175:+',  # Exon 1
                 'exon:chr1:225-250:+'):  # Exon 2
                ('junction:chr1:126-224:+',  # Exon1alt-Exon2 junction
                 'junction:chr1:176-224:+')}  # Exon1-Exon2 junction
        return true

    def test_ale(self, event_maker, strand):
        true = {('exon:chr1:300-350:+',  # Exon 3
                 'exon:chr1:400-425:+',  # Exon 4
                 'exon:chr1:475-500:+'):  # Exon 4 alt
                ('junction:chr1:351-399:+',  # Exon3-Exon4 junction
                 'junction:chr1:351-474:+')}  # Exon3-Exon4alt junction
        return true


@pytest.fixture
def graph_items(exon_start_stop, transcripts, chrom, strand):
    from outrigger.index.events import stringify_location, opposite

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
    items = tuple(items)
    return graph, items


class TestSpliceGraph(object):

    @pytest.fixture
    def splice_graph(self, junction_exon_triples):
        from outrigger.index.events import SpliceGraph

        return SpliceGraph(junction_exon_triples)

    def test___init__(self, junction_exon_triples, graph_items):
        from outrigger.index.events import SpliceGraph

        test = SpliceGraph(junction_exon_triples)

        graph, items = graph_items

        exons = tuple(junction_exon_triples.exon.unique())
        junctions = tuple(junction_exon_triples.junction.unique())

        pdt.assert_equal(test.exons, exons)
        pdt.assert_equal(test.junctions, junctions)
        pdt.assert_equal(sorted(test.items), sorted(items))

        assert_graph_items_equal(test.graph, test.items, graph, items)

    @pytest.fixture
    def exon1_i(self, strand):
        if strand == '+':
            return 0
        if strand == '-':
            return 3

    @pytest.fixture
    def exon1_name(self, splice_graph, exon1_i):
        return splice_graph.exons[exon1_i]

    def test_exons_one_junction_upstream(self, splice_graph, exon1_i, strand):
        test = tuple(splice_graph.exons_one_junction_downstream(exon1_i))
        if strand == '+':
            true = (1, 2, 5, 6, 3)
        if strand == '-':
            true = (2, 1, 0)
        assert test == true

    def test_exons_two_junctions_downstream(self, splice_graph, exon1_i,
                                            strand):
        test = tuple(splice_graph.exons_two_junctions_downstream(exon1_i))
        if strand == '+':
            true = (2, 3, 2, 3, 7)
        if strand == '-':
            true = (0, 1, 4, 0, 5, 6)
        assert test == true

    def test_junctions_between_exons(self, splice_graph, strand, exon1_i):
        if strand == '+':
            test = tuple(splice_graph.junctions_between_exons(exon1_i, 1))
            true = (8,)
        if strand == '-':
            test = tuple(splice_graph.junctions_between_exons(exon1_i, 0))
            true = (16,)
        assert test == true

    @pytest.fixture
    def skipped_exon_events(self, strand):
        if strand == '+':
            return {('exon:chr1:150-175:+', 'exon:chr1:200-250:+', 'exon:chr1:300-350:+'): ['junction:chr1:176-299:+',  # noqa
                                                                         'junction:chr1:176-199:+',  # noqa
                                                                         'junction:chr1:251-299:+'],  # noqa
 ('exon:chr1:150-175:+', 'exon:chr1:225-250:+', 'exon:chr1:300-350:+'): ['junction:chr1:176-299:+',  # noqa
                                                                         'junction:chr1:176-224:+',  # noqa
                                                                         'junction:chr1:251-299:+'],  # noqa
 ('exon:chr1:150-175:+', 'exon:chr1:225-250:+', 'exon:chr1:400-425:+'): ['junction:chr1:176-399:+',  # noqa
                                                                         'junction:chr1:176-224:+',  # noqa
                                                                         'junction:chr1:251-399:+'],  # noqa
 ('exon:chr1:150-175:+', 'exon:chr1:225-275:+', 'exon:chr1:300-350:+'): ['junction:chr1:176-299:+',  # noqa
                                                                         'junction:chr1:176-224:+',  # noqa
                                                                         'junction:chr1:276-299:+'],  # noqa
 ('exon:chr1:150-175:+', 'exon:chr1:300-350:+', 'exon:chr1:400-425:+'): ['junction:chr1:176-399:+',  # noqa
                                                                         'junction:chr1:176-299:+',  # noqa
                                                                         'junction:chr1:351-399:+']}  # noqa
        if strand == '-':
            return {('exon:chr1:400-425:-', 'exon:chr1:225-250:-', 'exon:chr1:150-175:-'): ['junction:chr1:176-399:-',  # noqa
                                                                         'junction:chr1:251-399:-',  # noqa
                                                                         'junction:chr1:176-224:-'],  # noqa
 ('exon:chr1:400-425:-', 'exon:chr1:300-350:-', 'exon:chr1:150-175:-'): ['junction:chr1:176-399:-',  # noqa
                                                                         'junction:chr1:351-399:-',  # noqa
                                                                         'junction:chr1:176-299:-'],  # noqa
 ('exon:chr1:400-425:-', 'exon:chr1:300-350:-', 'exon:chr1:225-250:-'): ['junction:chr1:251-399:-',  # noqa
                                                                         'junction:chr1:351-399:-',  # noqa
                                                                         'junction:chr1:251-299:-']}  # noqa

    @pytest.fixture
    def mutually_exclusive_events(self, strand):
        if strand == '+':
            return {('exon:chr1:150-175:+', 'exon:chr1:225-250:+', 'exon:chr1:300-350:+', 'exon:chr1:400-425:+'): ['junction:chr1:176-299:+',  # noqa
                                                                                                'junction:chr1:351-399:+',  # noqa
                                                                                                'junction:chr1:176-224:+',  # noqa
                                                                                                'junction:chr1:251-399:+']}  # noqa
        if strand == '-':
            return {('exon:chr1:400-425:-', 'exon:chr1:300-350:-', 'exon:chr1:225-250:-', 'exon:chr1:150-175:-'): ['junction:chr1:251-399:-',  # noqa
                                                                                                'junction:chr1:176-224:-',  # noqa
                                                                                                'junction:chr1:351-399:-',  # noqa
                                                                                                'junction:chr1:176-299:-']}  # noqa

    def test__skipped_exon(self, splice_graph, exon1_i, exon1_name,
                           skipped_exon_events):
        test = splice_graph._skipped_exon(exon1_i, exon1_name)
        true = skipped_exon_events
        pdt.assert_dict_equal(test, true)

    def test__mutually_exclusive_exon(self, splice_graph, exon1_i, exon1_name,
                                      mutually_exclusive_events):
        test = splice_graph._mutually_exclusive_exon(exon1_i, exon1_name)
        true = mutually_exclusive_events
        pdt.assert_dict_equal(test, true)

    def test_single_exon_alternative_events(self, splice_graph, exon1_i,
                                            exon1_name,
                                            mutually_exclusive_events,
                                            skipped_exon_events):
        test = splice_graph.single_exon_alternative_events(
            exon1_i, exon1_name)
        true = {'se': skipped_exon_events, 'mxe': mutually_exclusive_events}
        pdt.assert_dict_equal(test, true)
