def test_create_db(gtf_filename):
    from outrigger.io import gtf

    db = gtf.create_db(gtf_filename)

    assert False
