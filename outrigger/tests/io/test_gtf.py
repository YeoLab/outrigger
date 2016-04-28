
import gffutils
import pytest

def test_create_db(gtf_filename, db):
    from outrigger.io import gtf

    true = db

    test = gtf.create_db(gtf_filename)

    # Check that all the true db features are in the test database
    for featuretype in true.featuretypes():
        for feature in true.features_of_type(featuretype=featuretype):
            try:
                test[feature.id]
            except gffutils.FeatureNotFoundError:
                pytest.fail('Feature in true database not found in test '
                            'database')

    # Check that all the test db features are in the true database
    for featuretype in test.featuretypes():
        for feature in test.features_of_type(featuretype=featuretype):
            try:
                true[feature.id]
            except gffutils.FeatureNotFoundError:
                pytest.fail('Feature in test database not found in true '
                            'database')
