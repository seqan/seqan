#!/usr/bin/env python2.5
"""Tests for the module dddoc."""

import os
import sys
import unittest

import dddoc


class DataTest(unittest.TestCase):
    """Tests for the dddoc.Data class.

    These tests were originally written to support the refactoring and
    optimization of Data.find().
    """

    def setUp(self):
        self.lines_a_path = dddoc.Line(['a', 'path'], 'Text')
        self.lines_a_path_also = dddoc.Line(['a', 'path', 'also'], 'Also')
        self.lines_another_path = dddoc.Line(['another', 'path'], 'Another')

    def testQueryDepthOne(self):
        lines = [self.lines_a_path,
                 self.lines_a_path_also,
                 self.lines_another_path]
        data = dddoc.Data(lines, 0)
        result = data.find('a')
        self.assertEqual(result.lines, [self.lines_a_path, self.lines_a_path_also])
        self.assertEqual(result.level, 1)

    def testQueryDepthTwo(self):
        lines = [self.lines_a_path,
                 self.lines_a_path_also,
                 self.lines_another_path]
        data = dddoc.Data(lines, 0)
        result = data.find('a.path.also')
        self.assertEqual(result.lines, [self.lines_a_path_also])
        self.assertEqual(result.level, 3)

    def testQueryDepthOneVariant2(self):
        lines = [self.lines_a_path,
                 self.lines_a_path_also,
                 self.lines_another_path]
        data = dddoc.Data(lines, 0)
        result = data.find('another')
        self.assertEqual(result.lines, [self.lines_another_path])
        self.assertEqual(result.level, 1)

    def testWithLevelOf1(self):
        lines = [self.lines_a_path,
                 self.lines_a_path_also]
        data = dddoc.Data(lines, 1)
        result = data.find('a')
        self.assertEqual(result.lines, [])
        self.assertEqual(result.level, 2)

        result = data.find('path')
        self.assertEqual(result.lines, [self.lines_a_path, self.lines_a_path_also])
        self.assertEqual(result.level, 2)

if __name__ == '__main__':
    unittest.main()
