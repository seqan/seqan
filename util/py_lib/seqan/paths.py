#!/usr/bin/env python2
"""SeqAn path utilities.

Code to get paths within the SeqAn repository; Useful for setting defaults in
options.

When called as a program, the paths are printed.

The main assumption is that this module lies in '/util/py_lib/seqan.'

Copyright: (c) 2010, Knut Reinert, FU Berlin
License:   3-clause BSD (see LICENSE)
"""

__author__ = 'Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>'

import os
import os.path
import sys

def repositoryRoot():
    """Return path to directory root."""
    abs_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../..'))
    return os.path.relpath(abs_path, os.getcwd())

def pathToSkeletons():
    """Return path to '${REPOSITORY}/util/skel'."""
    return os.path.join(repositoryRoot(), 'util/skel')

def pathToTemplate(template, filename):
    """Return path to file with given name in given template."""
    return os.path.join(pathToSkeletons(), template, filename)

def pathToApp(location, app):
    """Return path to the app in the given location/repository."""
    return os.path.join(repositoryRoot(), location, 'apps', app)

def pathToInclude(location):
    """Return path to the include dir in the given location/repository."""
    return os.path.join(repositoryRoot(), location, 'include')

def pathToDemo(location, demo):
    """Return path to the demo file in the given location/repository."""
    return os.path.join(repositoryRoot(), location, 'demos', '%s.cpp' % demo)

def pathToTest(location, test):
    """Return path to the test in the given location/repository."""
    return os.path.join(repositoryRoot(), location, 'tests', test)

def pathToRepository(location):
    """Return path to the given location/repository."""
    return os.path.join(repositoryRoot(), location)

def pathToHeader(location, filename):
    """Return path to the given header - just concatenate."""
    return os.path.join(repositoryRoot(), location, filename)

def main(args):
    print 'SeqAn paths'
    print
    print 'repositoryRoot()   ==', repositoryRoot()
    print 'pathToSkeletons() ==', pathToSkeletons()

if __name__ == '__main__':
   sys.exit(main(sys.argv))
