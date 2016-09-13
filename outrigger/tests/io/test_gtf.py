
from ..conftest import check_db_equal

def test_create_db(gtf_filename, db, snap25_exon_id):
    from outrigger.io import gtf

    true = db
    test = gtf.create_db(gtf_filename, 'test.gtf.db')

    check_db_equal(test, true)

    # SNAP25 should be in both the true and test databases
    assert true[snap25_exon_id] is not None
    assert test[snap25_exon_id] is not None
