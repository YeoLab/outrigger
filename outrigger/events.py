import itertools
import logging
import sys

import graphlite
from graphlite import V
import numpy as np
import pandas as pd

from .region import Region
from .junctions import UPSTREAM, DOWNSTREAM, DIRECTIONS


def stringify_location(chrom, start, stop, strand, region=None):
    """"""
    if region is not None:
        return '{0}:{1}:{2}-{3}:{4}'.format(region, chrom, start, stop,
                                            strand)
    else:
        return '{0}:{1}-{2}:{3}'.format(chrom, start, stop, strand)


def opposite(direction):
    return UPSTREAM if direction == DOWNSTREAM else DOWNSTREAM


class EventMaker(object):

    def __init__(self, junction_exon_triples, db=None, junction_col='junction',
                 exon_col='exon'):
        """Combine splice junctions into splicing events

        Parameters
        ----------
        junction_exon_triples : pandas.DataFrame
            of "exon, direction, junction", e.g.:
            exon1, upstream, junction12

        db : gffutils.FeatureDB
            Gffutils Database of gene, transcript, and exon features. The exons
            must be accessible by the id provided on the `exon_col`
            columns. If not provided, certain splice types which require
            information about the transcript (AFE, ALE) cannot be annotated.
        """
        self.log = logging.getLogger('EventMaker')
        self.junction_exon_triples = junction_exon_triples
        self.db = db

        self.graph = graphlite.connect(":memory:", graphs=DIRECTIONS)
        self.exons = tuple(junction_exon_triples[exon_col].unique())
        self.n_exons = len(self.exons)
        self.junctions = tuple(junction_exon_triples[junction_col].unique())

        self.items = tuple(np.concatenate([self.exons, self.junctions]))
        self.item_to_region = pd.Series(map(Region, self.items),
                                        index=self.items)

        with self.graph.transaction() as tr:
            for i, row in self.junction_exon_triples.iterrows():
                junction = row[junction_col]
                exon = row[exon_col]

                junction_i = self.items.index(junction)
                exon_i = self.items.index(exon)

                self.log.debug('\n{} is {} of {}\n'.format(
                    exon, row.direction, junction))
                self.log.debug('{} is {} of {}\n'.format(
                    junction, opposite(row.direction), exon))

                tr.store(getattr(V(exon_i), row.direction)(junction_i))
                tr.store(getattr(V(junction_i),
                                 opposite(row.direction))(exon_i))

    @staticmethod
    def make_junction_exon_triples(junction_to_exons,
                                   junction_col='junction',
                                   upstream_col=UPSTREAM,
                                   downstream_col=DOWNSTREAM):
        """Create tidy table of exons upstream and downstream of a junction

        Parameters
        ----------
        sj_metadata : pandas.DataFrame
            A table with a column indicating "junction_location"
        junction_col : str
            Column name of the raw junction location (without |5p or |3p
            annotated)
        upstream_col : str
            Column name where exons upstream of the junction are stored
        downstream_col : str
            Column name where exons downstream of the junction are stored

        Returns
        -------
        junction_exon_triples
            A three-column table of junction_location, exon, and direction

        Examples
        --------
        >>> import pandas as pd
        >>> sj_metadata = pd.DataFrame(
        {'junction':['chr1:201-299:+', 'chr1:401:499:+'],
         'upstream': ['exon:chr1:100-200:+,exon:chr1:50-200:+',
                     'exon:chr1:300-400:+'],
         'downstream':['exon:chr1:300-400:+',
         'exon:chr1:500-600:+,exon:chr1:500-650:+']})
        >>> EventMaker.get_adjacent_exons(sj_metadata)

        """
        grouped = junction_to_exons.groupby(junction_col)
        direction_to_exon = {UPSTREAM: upstream_col,
                             DOWNSTREAM: downstream_col}
        dfs = []
        for direction, exon in direction_to_exon.items():
            df = grouped.apply(
                lambda x: x[exon].dropna().str.split(',').apply(pd.Series, 1))
            df = df.stack()
            df.index = df.index.droplevel((-2, -1))
            df = df.reset_index()
            df = df.rename(columns={0: 'exon'})
            df['direction'] = direction
            dfs.append(df)
        junction_exons = pd.concat(dfs, ignore_index=True)
        return junction_exons

    def event_dict_to_df(self, events, exon_names, junction_names):
        columns = list(exon_names) + list(junction_names) + ['event_id']
        data = pd.DataFrame(index=np.arange(len(events)), columns=columns)
        for i, (exons, junctions) in enumerate(events.items()):
            event_id = '@'.join(exons)
            data.loc[i, exon_names] = list(exons)
            data.loc[i, junction_names] = list(junctions)
            data.loc[i, 'event_id'] = event_id
        return data

    def skipped_exon(self):
        events = {}


        sys.stdout.write('Trying out {0} exons'
                         '...\n'.format(self.n_exons))
        for i, exon1_name in enumerate(self.exons):
            if (i + 1) % 10000 == 0:
                sys.stdout.write('\t{0}/{1} '
                                 'exons tested'.format(i + 1, self.n_exons))

            exon1_i = self.exons.index(exon1_name)
            exon23s = list(
                self.graph.find(
                    V().downstream(exon1_i)).traverse(V().upstream))
            exon23s = self.item_to_region[[self.items[i] for i in exon23s]]

            for exon_a, exon_b in itertools.combinations(exon23s, 2):
                if not exon_a.overlaps(exon_b):
                    exon2 = min((exon_a, exon_b), key=lambda x: x._start)
                    exon3 = max((exon_a, exon_b), key=lambda x: x._start)

                    exon2_i = self.exons.index(exon2.name)
                    exon3_i = self.exons.index(exon3.name)

                    exon23_junction_i = self.graph.find(
                        V(exon2_i).upstream).intersection(
                        V().upstream(exon3_i))
                    exon23_junction = [self.items[i] for i in
                                       set(exon23_junction_i)]
                    if len(exon23_junction) > 0:
                        # Isoform 1 - corresponds to Psi=0. Exclusion of exon2
                        exon13_junction = self.graph.find(
                            V(exon1_i).upstream) \
                            .intersection(V(exon3_i).downstream)

                        # Isoform 2 - corresponds to Psi=1. Inclusion of exon2
                        exon12_junction = self.graph.find(
                            V(exon1_i).upstream) \
                            .intersection(V(exon2_i).downstream)
                        exon23_junction = self.graph.find(
                            V(exon2_i).upstream) \
                            .intersection(V(exon3_i).downstream)

                        junctions_i = list(itertools.chain(
                            *[exon12_junction, exon23_junction,
                              exon13_junction]))
                        junctions = [self.items[i] for i in junctions_i]
                        exons = exon1_name, exon2.name, exon3.name

                        events[exons] = junctions
        events = self.event_dict_to_df(
            events, exon_names=['exon1', 'exon2', 'exon3'],
            junction_names=['junction12', 'junction23', 'junction13'])
        return events

    def mutually_exclusive_exon(self):
        events = {}

        for exon1_name in self.exons:
            exon1_i = self.items.index(exon1_name)

            exon23s_from1 = list(
                self.graph.find(V().downstream(
                    exon1_i)).traverse(V().upstream))
            exon4s = self.graph.find(V().downstream(exon1_i)).traverse(
                V().upstream).traverse(V().upstream).traverse(V().upstream)
            exon23s_from4 = exon4s.traverse(V().downstream).traverse(
                V().downstream)

            exon23s = set(exon23s_from4) & set(exon23s_from1)
            exon23s = [self.items[i] for i in exon23s]

            exon23s = self.item_to_region[exon23s]

            for exon_a, exon_b in itertools.combinations(exon23s, 2):
                if not exon_a.overlaps(exon_b):
                    exon2 = min((exon_a, exon_b), key=lambda x: x._start)
                    exon3 = max((exon_a, exon_b), key=lambda x: x._start)

                    exon2_i = self.items.index(exon2.name)
                    exon3_i = self.items.index(exon3.name)

                    exon4_from2 = set(
                        self.graph.find(V(exon2_i).upstream).traverse(
                            V().upstream))
                    exon4_from3 = set(
                        self.graph.find(V(exon3_i).upstream).traverse(
                            V().upstream))
                    try:
                        exon4_i = (exon4_from2 & exon4_from3).pop()
                        exon4_name = self.items[exon4_i]
                        # Isoform 1 - corresponds to Psi=0. Inclusion of exon3
                        exon13_junction = self.graph.find(
                            V(exon1_i).upstream).intersection(
                            V(exon3_i).downstream)
                        exon34_junction = self.graph.find(
                            V(exon3_i).upstream) \
                            .intersection(V(exon4_i).downstream)

                        # Isoform 2 - corresponds to Psi=1. Inclusion of exon2
                        exon12_junction = self.graph.find(
                            V(exon1_i).upstream).intersection(
                            V(exon2_i).downstream)
                        exon24_junction = self.graph.find(
                            V(exon2_i).upstream) \
                            .intersection(V(exon4_i).downstream)

                        exon_tuple = exon1_name, exon2.name, exon3.name, \
                            exon4_name
                        #             print exon12_junction.next()
                        junctions = list(
                            itertools.chain(*[exon13_junction,
                                              exon34_junction,
                                              exon12_junction,
                                              exon24_junction]))
                        junctions = [self.items[i] for i in junctions]

                        events[exon_tuple] = junctions
                    except:
                        pass
        events = self.event_dict_to_df(events,
                                       exon_names=['exon1', 'exon2', 'exon3',
                                                   'exon4'],
                                       junction_names=['junction13',
                                                       'junction34',
                                                       'junction12',
                                                       'junction24'])
        return events


