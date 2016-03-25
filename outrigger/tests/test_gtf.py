def test_create_db(gtf_filename):
    from outrigger.index import gtf

    db = gtf.create_db(gtf_filename)

    assert False
