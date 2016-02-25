#!/usr/bin/env python2
"""Tests for the lexer module.

We got the lexer itself from an external source so we do not write tests for
it now as long as it works nicely.  However, we need some extension for
storing the precise token locations.
"""

import unittest

import seqan.dox.lexer as lexer


class TestLexer(unittest.TestCase):
    def testWithoutOffset(self):
        self._testWithOffset(0, 0)

    def testWithOffset(self):
        self._testWithOffset(3, 4)

    def _testWithOffset(self, line_offset, col_offset):
        rules = (('number', r'[0-9]+'),
                 ('word', r'\w+'))
        txt = ('word number 0123\n'
               'first second\n'
               'third fourth')
        lex = lexer.Lexer(rules, line_offset=line_offset, col_offset=col_offset)
        lex.input(txt)
        tokens = [t for t in lex.tokens()]
        #for t in tokens: print t;

        self.assertEqual(tokens[0].val, 'word')
        self.assertEqual(tokens[0].type, 'word')
        self.assertEqual(tokens[0].pos, 0)
        self.assertEqual(tokens[0].lineno, 0 + line_offset)
        self.assertEqual(tokens[0].column, 0 + col_offset)

        self.assertEqual(tokens[1].val, 'number')
        self.assertEqual(tokens[1].type, 'word')
        self.assertEqual(tokens[1].pos, 5)
        self.assertEqual(tokens[1].lineno, 0 + line_offset)
        self.assertEqual(tokens[1].column, 5 + col_offset)

        self.assertEqual(tokens[2].val, '0123')
        self.assertEqual(tokens[2].type, 'number')
        self.assertEqual(tokens[2].pos, 12)
        self.assertEqual(tokens[2].lineno, 0 + line_offset)
        self.assertEqual(tokens[2].column, 12 + col_offset)

        self.assertEqual(tokens[3].val, 'first')
        self.assertEqual(tokens[3].type, 'word')
        self.assertEqual(tokens[3].pos, 17)
        self.assertEqual(tokens[3].lineno, 1 + line_offset)
        self.assertEqual(tokens[3].column, 0 + col_offset)

        self.assertEqual(tokens[4].val, 'second')
        self.assertEqual(tokens[4].type, 'word')
        self.assertEqual(tokens[4].pos, 23)
        self.assertEqual(tokens[4].lineno, 1 + line_offset)
        self.assertEqual(tokens[4].column, 6 + col_offset)

        self.assertEqual(tokens[5].val, 'third')
        self.assertEqual(tokens[5].type, 'word')
        self.assertEqual(tokens[5].pos, 30)
        self.assertEqual(tokens[5].lineno, 2 + line_offset)
        self.assertEqual(tokens[5].column, 0 + col_offset)

        self.assertEqual(tokens[6].val, 'fourth')
        self.assertEqual(tokens[6].type, 'word')
        self.assertEqual(tokens[6].pos, 36)
        self.assertEqual(tokens[6].lineno, 2 + line_offset)
        self.assertEqual(tokens[6].column, 6 + col_offset)


if __name__ == '__main__':
    unittest.main()
