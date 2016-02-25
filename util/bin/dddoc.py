#!/usr/bin/env python2
"""SeqAn Documentation System DDDoc."""

import os.path
import sys

def main():
    path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'py_lib'))
    sys.path.insert(0, path)
    import seqan.dddoc
    return seqan.dddoc.main()

if __name__ == '__main__':
    sys.exit(main())
