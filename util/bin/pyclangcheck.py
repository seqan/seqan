#!/usr/bin/env python2
"""llvm-clang based style checker."""

__author__ = 'Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>'

import os.path
import sys

def main():
    path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'py_lib'))
    sys.path.insert(0, path)
    import seqan.pyclangcheck
    return seqan.pyclangcheck.main()

if __name__ == '__main__':
    ## import cProfile
    ## cProfile.run('main()', 'pyclangcheck.prof')
    sys.exit(main())
    ## import pstats
    ## p = pstats.Stats('pyclangcheck.prof')
    ## print 'BY TIME'
    ## p.sort_stats('time').print_stats()
    ## print 'BY CUMULATIVE'
    ## p.sort_stats('cumulative').print_stats()
