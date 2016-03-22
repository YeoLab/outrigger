import numpy as np
import pytest
import pandas as pd
import pandas.util.testing as pdt
import six
import gffutils


def test_create_db(gtf_filename):
    from outrigger import gtf

    db = gtf.create_db(gtf_filename)

    assert False
