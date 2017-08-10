# -*- coding: utf-8 -*-

from __future__ import absolute_import

from os.path import abspath, dirname, join

try:
    import pytest
    import pytest_benchmark
except ImportError:
    pytest = None

thermodynamics_directory = abspath(join(dirname(abspath(__file__)), ".."))
thermodynamics_location = abspath(join(thermodynamics_directory, ".."))
data_dir = join(thermodynamics_directory, "thermodynamics_data", "")

def test_all(args=None):
    """ alias for running all unit-tests on installed thermodynamics
    """
    if pytest:
        args = args if args else []

        return pytest.main(
            ['--pyargs', 'thermodynamics', '--benchmark-skip', '-v', '-rs'] + args
        )
    else:
        raise ImportError('missing package pytest and pytest_benchmark'
                          ' required for testing')