#!/usr/bin/env python
"""Build SeqAn Release Page."""

import os.path
import sys

def main():
    path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'py_lib'))
    sys.path.insert(0, path)
    import seqan.release_page
    return seqan.release_page.main()

if __name__ == '__main__':
    sys.exit(main())
