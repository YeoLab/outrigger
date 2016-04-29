import os

import pandas as pd
import pandas.util.testing as pdt
import pytest

class TestExonJunctionAdjacencies(object):

    @pytest.fixture
    def exon_id(self):
        return 'exon:chr7:130200983-130201173:-'

    @pytest.fixture
    def exon(self, db, exon_id):
        return db[exon_id]

    @pytest.fixture
    def adjacencies(self, metadata, db):
        from outrigger.index.adjacencies import ExonJunctionAdjacencies

        return ExonJunctionAdjacencies(metadata, db)

    @pytest.fixture
    def adjacent_in_genome_template(self, treutlein_adjacencies):
        return os.path.join(
            treutlein_adjacencies,
            'junctions_genome_adjacent_to_exon_{}stream.csv')

    @pytest.fixture
    def adjacent_in_genome_upstream(self, adjacent_in_genome_template):
        return pd.read_csv(adjacent_in_genome_template.format('up'),
                           squeeze=True, index_col=0)
    @pytest.fixture
    def adjacent_in_genome_downstream(self, adjacent_in_genome_template):
        return pd.read_csv(adjacent_in_genome_template.format('down'),
                           squeeze=True, index_col=0)

    @pytest.fixture
    def adjacent_in_genome(self, adjacent_in_genome_upstream,
                           adjacent_in_genome_downstream):
        """Dict of upstream and downstream boolean junctions"""
        return {'upstream': adjacent_in_genome_upstream,
                'downstream': adjacent_in_genome_downstream}

    def test___init(self, metadata, db):
        from outrigger.index.adjacencies import ExonJunctionAdjacencies
        from outrigger.io.common import (JUNCTION_ID, EXON_START, EXON_STOP,
                                         CHROM, STRAND)

        adjacencies = ExonJunctionAdjacencies(metadata, db)

        true_metadata = metadata.copy()
        true_metadata = true_metadata.set_index(JUNCTION_ID)
        true_metadata = true_metadata.sort_index()

        pdt.assert_frame_equal(adjacencies.metadata, true_metadata)

        assert adjacencies.junction_id == JUNCTION_ID
        assert adjacencies.exon_start == EXON_START
        assert adjacencies.exon_stop == EXON_STOP
        assert adjacencies.chrom == CHROM
        assert adjacencies.strand == STRAND

        assert adjacencies.db == db

    @pytest.mark.xfail
    def test___init_missing_required_column(self, metadata, db):
        from outrigger.index.adjacencies import ExonJunctionAdjacencies
        from outrigger.io.common import JUNCTION_ID

        test_metadata = metadata.copy()
        test_metadata = test_metadata.drop(JUNCTION_ID, axis=1)
        ExonJunctionAdjacencies(test_metadata, db)

    def test__single_junction_exon_triple(self, adjacencies, exon_id,
                                          treutlein_adjacencies,
                                          adjacent_in_genome_upstream):
        test = adjacencies._single_junction_exon_triple(
            adjacent_in_genome_upstream, 'downstream', exon_id)

        true = pd.read_csv(os.path.join(treutlein_adjacencies,
                                        'single_junction_exon_triple.csv'))
        pdt.assert_frame_equal(test, true)

    def test__to_stranded_transcript_adjacency(self, adjacencies, strand,
                                               adjacent_in_genome):
        test = adjacencies._to_stranded_transcript_adjacency(
            adjacent_in_genome, strand)

        if strand == '-':
            pdt.assert_series_equal(test['upstream'],
                                    adjacent_in_genome['downstream'])
            pdt.assert_series_equal(test['downstream'],
                                    adjacent_in_genome['upstream'])
        elif strand == '+':
            pdt.assert_series_equal(test['upstream'],
                                    adjacent_in_genome['upstream'])
            pdt.assert_series_equal(test['downstream'],
                                    adjacent_in_genome['downstream'])

    def test__junctions_genome_adjacent_to_exon(self, adjacencies, exon,
                                                adjacent_in_genome):
        test = adjacencies._junctions_genome_adjacent_to_exon(exon)



        pdt.assert_dict_equal(test, adjacent_in_genome)

    def test__adjacent_junctions_single_exon(self, adjacencies, exon,
                                             treutlein_adjacencies):
        test = adjacencies._adjacent_junctions_single_exon(exon)

        true = pd.read_csv(os.path.join(treutlein_adjacencies,
                                        'adjacent_junctions_single_exon.csv'))
        pdt.assert_frame_equal(test, true)

