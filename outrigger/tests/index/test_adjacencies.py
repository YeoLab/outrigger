import os

import pandas as pd
import pandas.util.testing as pdt
import pytest

class TestExonJunctionAdjacencies(object):

    exon_id = 'exon:chr7:130200983-130201173:-'

    @pytest.fixture
    def exon(self, db):
        return db[self.exon_id]

    @pytest.fixture
    def adjacencies(self, metadata, db):
        from outrigger.index.adjacencies import ExonJunctionAdjacencies

        return ExonJunctionAdjacencies(metadata, db)

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

    def test__single_junction_exon_triple(self):
        pass

    def test__to_stranded_transcript_adjacency(self):
        pass

    def test__junctions_genome_adjacent_to_exon(self, adjacencies, exon,
                                                treutlein_adjacencies):
        test = adjacencies._junctions_genome_adjacent_to_exon(exon)

        template = os.path.join(
            treutlein_adjacencies,
            'junctions_genome_adjacent_to_exon_{}stream.csv')
        true_upstream = pd.read_csv(template.format('up'), squeeze=True,
                                    index_col=0)
        true_downstream = pd.read_csv(template.format('down'), squeeze=True,
                                      index_col=0)

        true = {'upstream': true_upstream, 'downstream': true_downstream}

        pdt.assert_dict_equal(test, true)

    def test__adjacent_junctions_single_exon(self, adjacencies, exon):
        test = adjacencies._adjacent_junctions_single_exon(exon)
        assert False
