#!/usr/bin/env python2
"""Tests for the path manager."""

__author__ = 'Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>'


import os
import os.path
import unittest

import seqan.dox.file_mgr as file_mgr


class Test(unittest.TestCase):
    def setUp(self):
        self.base_dir = os.path.dirname(os.path.realpath(__file__))
        self.mgr = file_mgr.FileManager()

    def testSimpleExample(self):
        p = os.path.join(self.base_dir, 'test.cpp')
        f = self.mgr.loadFile(p)
        self.assertEqual(f.path, p)
        self.assertEqual(f.start_markers, ['*!'])
        c = file_mgr.Comment(6, 0, 62, 8, 2, 87, 7, 3, 'Test comment...\n',
                             '/*!\n * Test comment...\n*/')
        self.assertEqual(f.comments, [c])


if __name__ == '__main__':
    unittest.main()
