#!/usr/bin/env python2
"""Tests for the include manager."""

__author__ = 'Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>'


import os
import os.path
import unittest

import seqan.dox.inc_mgr as inc_mgr


class TestIncludeManager(unittest.TestCase):
    def setUp(self):
        base_dir = os.path.dirname(os.path.realpath(__file__))
        base_dir = os.path.join(base_dir, '../', 'test_src')
        self.mgr = inc_mgr.IncludeManager([base_dir])

    def testIncludeFile(self):
        txt = self.mgr.loadFile('example.cpp')
        self.assert_(txt.splitlines()[0].startswith('#include <iostream>'))
        self.assert_(txt.splitlines()[-1].endswith('}'))

    def testIncludeSnippet(self):
        txt = self.mgr.loadSnippet('example.cpp', 'Print to stdout')
        self.assertEqual(len(txt.splitlines()), 1)
        self.assertEqual(txt.splitlines()[0], r'    std::cout << "This is an example.\n";')


if __name__ == '__main__':
    unittest.main()
