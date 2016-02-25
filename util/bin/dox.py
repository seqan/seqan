#!/usr/bin/env python2
"""SeqAn doxygen-style documentation system."""

import os.path
import sys

def main():
    path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'py_lib'))
    sys.path.insert(0, path)
    import seqan.dox.pure
    return seqan.dox.pure.main()

if __name__ == '__main__':
    sys.exit(main())
