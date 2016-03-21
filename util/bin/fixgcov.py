#!/usr/bin/env python2
"""Fix gcov output.

Fix gcov output with templates.  This is done by first parsing in the .cpp files
(compilation units) with libclang.  The AST is then parsed and all lines within
composite statements ({ stmt; stmt; ... }) are memoized as 'interesting' lines.
The resulting interesting lines are serialized to a location file with pickle.
Finally, gcov output files are read and updated.  If a line is interesting but
not marked as covered or uncovered (marker '-'), it is marked as uncovered
(marker '#####').

See API documentation of seqan.fixgcov.app for more information.
"""

__author__ = 'Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>'

import os.path
import sys


def main():
    path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'py_lib'))
    sys.path.insert(0, path)
    import seqan.fixgcov
    return seqan.fixgcov.main()


if __name__ == '__main__':
    sys.exit(main())
