#!/usr/bin/env python2
"""Sequential tests (num_threads == 0) for RazerS 3.

See run_tests.py for the actual tests.  We only call the code from there.
"""

import logging
import os.path
import sys

# Automagically add util/py_lib to PYTHONPATH environment variable.
path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..',
                                    '..', '..', 'util', 'py_lib'))
sys.path.insert(0, path)

import seqan.app_tests as app_tests

import run_tests

if __name__ == '__main__':
    sys.exit(app_tests.main(run_tests.main, num_threads=0))
