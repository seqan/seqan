#!/usr/bin/env python2
"""Tests for the dox_parser module.

The parser results are very complex.  We rely on parsing Doxygen-style
documentation and then dumping it back into this format.

We write the "full" test with the smallest amount of empty lines possible to
check that the tokenization works well in these corner cases.  Having empty
lines actually is the easier case.
"""

__author__ = 'Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>'

# TODO(holtgrew): Add tests for @implements @extends

import os
import unittest
import sys

import seqan.dox.lexer as lexer
import seqan.dox.dox_tokens as dox_tokens
import seqan.dox.dox_parser as dox_parser


class TestDoxParserBase(unittest.TestCase):
    """Base class for all dox parser tests."""

    def __init__(self, *args, **kwargs):
        unittest.TestCase.__init__(self, *args, **kwargs)
        self.maxDiff = 1024*1024

    def createLexer(self, text):
        """Create a lexer.Lexer object with the given text."""
        lex = lexer.Lexer(dox_tokens.LEXER_TOKENS, skip_whitespace=False)
        lex.input(text)
        return lex

    def parseText(self, text):
        """Create a dox parser and let it parse the given text.

        Return the created parser.
        """
        parser = dox_parser.Parser()
        parser.parse(self.createLexer(text))
        return parser


class TestClassParsing(TestDoxParserBase):
    def testMinimal(self):
        txt = '@class Klass'
        parser = self.parseText(txt)
        doc = parser.documentation
        self.assertMultiLineEqual(doc.getFormatted(), '@class Klass\n\n')

    def testTwoMinimal(self):
        txt = ('@class A\n'
               '@brief Brief A\n'
               '@class B\n')
        parser = self.parseText(txt)
        doc = parser.documentation
        self.assertMultiLineEqual(doc.getFormatted(), '@class A\n\n@brief Brief A\n\n\n\n@class B\n\n\n\n')

    def testFull(self):
        txt = ('@class Klass\n'
               '@implements Concept\n'
               '@headerfile <seqan/base.h>\n'
               '@headerfile <seqan/sequence.h>\n'
               '@extends Other Klass\n'
               '@brief This is a brief text.\n'
               '@deprecated Deprecation message.\n'
               '@signature template <typename T1, typename T2>\n'
               '           class Klass;\n'
               '@tparam T1 The first value and a very very very very long \n'
               '           description\n'
               '@tparam T2 The second value and a very very very very long \n'
               '           description\n'
               '\n'
               '@section This is the first section.\n'
               '\n'
               'This is the first paragraph.\n'
               '\n'
               '@subsection This is the first subsection.\n'
               '\n'
               'This is the second paragraph.\n'
               '@see Other')
        parser = self.parseText(txt)
        doc = parser.documentation
        expected = ('@class Klass\n'
                    '\n'
                    '@implements Concept\n'
                    '\n'
                    '@extends Other Klass\n'
                    '\n'
                    '@headerfile <seqan/base.h>\n'
                    '@headerfile <seqan/sequence.h>\n'
                    '\n'
                    '@brief This is a brief text.\n'
                    '\n'
                    '@deprecated Deprecation message.\n'
                    '\n'
                    '@signature template <typename T1, typename T2>\n'
                    '           class Klass;\n'
                    '\n'
                    '@tparam T1 The first value and a very very very\n'
                    '           very long description\n'
                    '@tparam T2 The second value and a very very very\n'
                    '           very long description\n'
                    '\n'
                    '@section This is the first section.\n'
                    '\n'
                    'This is the first paragraph.\n'
                    '\n'
                    '@subsection This is the first subsection.\n'
                    '\n'
                    'This is the second paragraph.\n'
                    '\n'
                    '@see Other\n'
                    '\n')
        self.assertMultiLineEqual(doc.getFormatted(50), expected)


class TestTypedefParsing(TestDoxParserBase):
    def testGlobalMinimal(self):
        txt = '@typedef Typedef'
        parser = self.parseText(txt)
        doc = parser.documentation
        self.assertMultiLineEqual(doc.getFormatted(), '@typedef Typedef\n\n')

    def testGlobalFull(self):
        txt = ('@typedef Typedef\n'
               '@brief This is an example for a typedef.\n'
               '@deprecated Deprecation message.\n'
               '@signature typedef Foo<Bar> Typedef;\n'
               '\n'
               '@section This is the first section.\n'
               '\n'
               'This is the first paragraph.\n'
               '\n'
               '@subsection This is the first subsection.\n'
               '\n'
               'This is the second paragraph.\n'
               '@see Other')
        parser = self.parseText(txt)
        doc = parser.documentation
        expected = ('@typedef Typedef\n'
                    '\n'
                    '@brief This is an example for a typedef.\n'
                    '\n'
                    '@deprecated Deprecation message.\n'
                    '\n'
                    '@signature typedef Foo<Bar> Typedef;\n'
                    '\n'
                    '@section This is the first section.\n'
                    '\n'
                    'This is the first paragraph.\n'
                    '\n'
                    '@subsection This is the first subsection.\n'
                    '\n'
                    'This is the second paragraph.\n'
                    '\n'
                    '@see Other\n'
                    '\n')
        self.assertMultiLineEqual(doc.getFormatted(50), expected)


class TestAdaptionParsing(TestDoxParserBase):
    def testGlobalMinimal(self):
        txt = '@adaption Adaption'
        parser = self.parseText(txt)
        doc = parser.documentation
        self.assertMultiLineEqual(doc.getFormatted(), '@adaption Adaption\n\n')

    def testGlobalFull(self):
        txt = ('@adaption Adaption\n'
               '@brief This is an example for an adaption.\n'
               '@deprecated Deprecation message.\n'
               '@signature std::string;\n'
               '\n'
               '@section This is the first section.\n'
               '\n'
               'This is the first paragraph.\n'
               '\n'
               '@subsection This is the first subsection.\n'
               '\n'
               'This is the second paragraph.\n'
               '@see Other')
        parser = self.parseText(txt)
        doc = parser.documentation
        expected = ('@adaption Adaption\n'
                    '\n'
                    '@brief This is an example for an adaption.\n'
                    '\n'
                    '@deprecated Deprecation message.\n'
                    '\n'
                    '@signature std::string;\n'
                    '\n'
                    '@section This is the first section.\n'
                    '\n'
                    'This is the first paragraph.\n'
                    '\n'
                    '@subsection This is the first subsection.\n'
                    '\n'
                    'This is the second paragraph.\n'
                    '\n'
                    '@see Other\n'
                    '\n')
        self.assertMultiLineEqual(doc.getFormatted(50), expected)


class TestMacroParsing(TestDoxParserBase):
    def testGlobalMinimal(self):
        txt = '@macro MACRO'
        parser = self.parseText(txt)
        doc = parser.documentation
        self.assertMultiLineEqual(doc.getFormatted(), '@macro MACRO\n\n')

    def testGlobalFull(self):
        txt = ('@macro MACRO\n'
               '@brief This is an example for a macro.\n'
               '@deprecated Deprecation message.\n'
               '@signature MACRO(param)\n'
               '@param param The parameter.\n'
               '@return TString A path as <tt>char const *</tt>.\n'
               '@throw Exception The exception type.\n'
               '@datarace This macro is not thread safe.\n'
               '\n'
               '@section This is the first section.\n'
               '\n'
               'This is the first paragraph.\n'
               '\n'
               '@subsection This is the first subsection.\n'
               '\n'
               'This is the second paragraph.\n'
               '@see Other')
        parser = self.parseText(txt)
        doc = parser.documentation
        expected = ('@macro MACRO\n'
                    '\n'
                    '@brief This is an example for a macro.\n'
                    '\n'
                    '@deprecated Deprecation message.\n'
                    '\n'
                    '@signature MACRO(param)\n'
                    '\n'
                    '@param param The parameter.\n'
                    '\n'
                    '@return TString A path as <tt>char const *</tt>.\n'
                    '\n'
                    '@throw Exception The exception type.\n'
                    '\n'
                    '@datarace This macro is not thread safe.\n'
                    '\n'
                    '@section This is the first section.\n'
                    '\n'
                    'This is the first paragraph.\n'
                    '\n'
                    '@subsection This is the first subsection.\n'
                    '\n'
                    'This is the second paragraph.\n'
                    '\n'
                    '@see Other\n'
                    '\n')
        self.assertMultiLineEqual(doc.getFormatted(50), expected)


class TestFunctionParsing(TestDoxParserBase):
    def testGlobalMinimal(self):
        txt = '@fn funktion'
        parser = self.parseText(txt)
        doc = parser.documentation
        self.assertMultiLineEqual(doc.getFormatted(), '@fn funktion\n\n')

    def testReturnValue(self):
        txt = ('@fn funktion\n'
               '@return bool <tt>true</tt> if empty, <tt>false</tt> otherwise.')
        parser = self.parseText(txt)
        doc = parser.documentation
        fn = doc.entries[0]
        self.assertEqual(len(fn.returns), 1)
        self.assertEqual(fn.returns[0].name.text, 'bool')
        txt = '<tt>true</tt> if empty, <tt>false</tt> otherwise.'
        self.assertEqual(fn.returns[0].text.text, txt)

    def testThrow(self):
        txt = ('@fn funktion\n'
               '@throw Exception The thrown exception')
        parser = self.parseText(txt)
        doc = parser.documentation
        fn = doc.entries[0]
        self.assertEqual(len(fn.throws), 1)
        self.assertEqual(fn.throws[0].name.text, 'Exception')
        txt = 'The thrown exception'
        self.assertEqual(fn.throws[0].text.text, txt)
        
    def testDataRace(self):
        txt = ('@fn funktion\n'
               '@datarace This function is thread safe.')
        parser = self.parseText(txt)
        doc = parser.documentation
        fn = doc.entries[0]
        self.assertEqual(len(fn.dataraces), 1)
        txt = 'This function is thread safe.'
        self.assertEqual(fn.dataraces[0].text.text, txt)

    def testGlobalFull(self):
        txt = ('@fn funktion\n'
               '@brief This is a brief text.\n'
               '@deprecated Deprecation message.\n'
               '@signature TRes1 funktion<T1>(TParam p1);\n'
               '@signature TRes2 funktion<T2>(TParam p2);\n'
               '@tparam T1 The first value and a very very very very long \n'
               '           description\n'
               '@tparam T2 The second value and a very very very very long \n'
               '           description\n'
               '@param[in] p1 The first parameter.\n'
               '@param     p2 The second parameter.\n'
               '@return TRes1 The first return type.\n'
               '@return TRes2 The second return type.\n'
               '\n'
               '@throw Exception The thrown exception.\n'
               '\n'
               '@datarace This function is thread safe.\n'
               '\n'
               '@section This is the first section.\n'
               '\n'
               'This is the first paragraph.\n'
               '\n'
               '@subsection This is the first subsection.\n'
               '\n'
               'This is the second paragraph.\n'
               '@see Other')
        parser = self.parseText(txt)
        doc = parser.documentation
        expected = ('@fn funktion\n'
                    '\n'
                    '@brief This is a brief text.\n'
                    '\n'
                    '@deprecated Deprecation message.\n'
                    '\n'
                    '@signature TRes1 funktion<T1>(TParam p1);\n'
                    '@signature TRes2 funktion<T2>(TParam p2);\n'
                    '\n'
                    '@tparam T1 The first value and a very very very\n'
                    '           very long description\n'
                    '@tparam T2 The second value and a very very very\n'
                    '           very long description\n'
                    '\n'
                    '@param[in] p1 The first parameter.\n'
                    '@param p2 The second parameter.\n'
                    '\n'
                    '@return TRes1 The first return type.\n'
                    '@return TRes2 The second return type.\n'
                    '\n'
                    '@throw Exception The thrown exception.\n'
                    '\n'
                    '@datarace This function is thread safe.\n'
                    '\n'
                    '@section This is the first section.\n'
                    '\n'
                    'This is the first paragraph.\n'
                    '\n'
                    '@subsection This is the first subsection.\n'
                    '\n'
                    'This is the second paragraph.\n'
                    '\n'
                    '@see Other\n'
                    '\n')
        self.assertMultiLineEqual(doc.getFormatted(50), expected)

    def testInterfaceMinimal(self):
        txt = '@fn Klass#funktion'
        parser = self.parseText(txt)
        doc = parser.documentation
        self.assertMultiLineEqual(doc.getFormatted(), '@fn Klass#funktion\n\n')

    def testInterfaceFull(self):
        txt = ('@fn Klass#funktion\n'
               '@brief This is a brief text.\n'
               '@deprecated Deprecation message.\n'
               '@signature TRes1 funktion<T1>(TParam p1);\n'
               '@signature TRes2 funktion<T2>(TParam p2);\n'
               '@tparam T1 The first value and a very very very very long \n'
               '           description\n'
               '@tparam T2 The second value and a very very very very long \n'
               '           description\n'
               '@param[in] p1 The first parameter.\n'
               '@param     p2 The second parameter.\n'
               '@return TRes1 The first return type.\n'
               '@return TRes2 The second return type.\n'
               '@datarace This function is thread safe.\n'
               '\n'
               '@section This is the first section.\n'
               '\n'
               'This is the first paragraph.\n'
               '\n'
               '@subsection This is the first subsection.\n'
               '\n'
               'This is the second paragraph.\n'
               '@see Other')
        parser = self.parseText(txt)
        doc = parser.documentation
        expected = ('@fn Klass#funktion\n'
                    '\n'
                    '@brief This is a brief text.\n'
                    '\n'
                    '@deprecated Deprecation message.\n'
                    '\n'
                    '@signature TRes1 funktion<T1>(TParam p1);\n'
                    '@signature TRes2 funktion<T2>(TParam p2);\n'
                    '\n'
                    '@tparam T1 The first value and a very very very\n'
                    '           very long description\n'
                    '@tparam T2 The second value and a very very very\n'
                    '           very long description\n'
                    '\n'
                    '@param[in] p1 The first parameter.\n'
                    '@param p2 The second parameter.\n'
                    '\n'
                    '@return TRes1 The first return type.\n'
                    '@return TRes2 The second return type.\n'
                    '\n'
                    '@datarace This function is thread safe.\n'
                    '\n'
                    '@section This is the first section.\n'
                    '\n'
                    'This is the first paragraph.\n'
                    '\n'
                    '@subsection This is the first subsection.\n'
                    '\n'
                    'This is the second paragraph.\n'
                    '\n'
                    '@see Other\n'
                    '\n')
        self.assertMultiLineEqual(doc.getFormatted(50), expected)


class TestConceptParsing(TestDoxParserBase):
    def testMinimal(self):
        txt = '@concept Konzept'
        parser = self.parseText(txt)
        doc = parser.documentation
        self.assertMultiLineEqual(doc.getFormatted(), '@concept Konzept\n\n')

    def testFull(self):
        txt = ('@concept Konzept\n'
               '@extends Konzert\n'
               '@headerfile <seqan/base.h>\n'
               '@headerfile <seqan/sequence.h>\n'
               '@brief This is a brief text.\n'
               '@deprecated Deprecation message.\n'
               '@signature concept Konzept;\n'
               '\n'
               '@section This is the first section.\n'
               '\n'
               'This is the first paragraph.\n'
               '\n'
               '@subsection This is the first subsection.\n'
               '\n'
               'This is the second paragraph.\n'
               '@see Other')
        parser = self.parseText(txt)
        doc = parser.documentation
        expected = ('@concept Konzept\n'
                    '\n'
                    '@extends Konzert\n'
                    '\n'
                    '@headerfile <seqan/base.h>\n'
                    '@headerfile <seqan/sequence.h>\n'
                    '\n'
                    '@brief This is a brief text.\n'
                    '\n'
                    '@deprecated Deprecation message.\n'
                    '\n'
                    '@signature concept Konzept;\n'
                    '\n'
                    '@section This is the first section.\n'
                    '\n'
                    'This is the first paragraph.\n'
                    '\n'
                    '@subsection This is the first subsection.\n'
                    '\n'
                    'This is the second paragraph.\n'
                    '\n'
                    '@see Other\n'
                    '\n')
        self.assertMultiLineEqual(doc.getFormatted(50), expected)


class TestMetafunctionParsing(TestDoxParserBase):
    def testGlobalMinimal(self):
        txt = '@mfn Metafunktion'
        parser = self.parseText(txt)
        doc = parser.documentation
        self.assertMultiLineEqual(doc.getFormatted(), '@mfn Metafunktion\n\n')

    def testGlobalFull(self):
        txt = ('@mfn Metafunktion\n'
               '@brief This is a brief text.\n'
               '@deprecated Deprecation message.\n'
               '@signature Metafunktion<T1>::Type;\n'
               '@signature Metafunktion<T2>::VALUE;\n'
               '@tparam T1 The first value and a very very very very long \n'
               '           description\n'
               '@tparam T2 The second value and a very very very very long \n'
               '           description\n'
               '@return Type The return type.\n'
               '@return VALUE The return value.\n'
               '\n'
               '@section This is the first section.\n'
               '\n'
               'This is the first paragraph.\n'
               '\n'
               '@subsection This is the first subsection.\n'
               '\n'
               'This is the second paragraph.\n'
               '@see Other')
        parser = self.parseText(txt)
        doc = parser.documentation
        expected = ('@mfn Metafunktion\n'
                    '\n'
                    '@brief This is a brief text.\n'
                    '\n'
                    '@deprecated Deprecation message.\n'
                    '\n'
                    '@signature Metafunktion<T1>::Type;\n'
                    '@signature Metafunktion<T2>::VALUE;\n'
                    '\n'
                    '@tparam T1 The first value and a very very very\n'
                    '           very long description\n'
                    '@tparam T2 The second value and a very very very\n'
                    '           very long description\n'
                    '\n'
                    '@return Type The return type.\n'
                    '@return VALUE The return value.\n'
                    '\n'
                    '@section This is the first section.\n'
                    '\n'
                    'This is the first paragraph.\n'
                    '\n'
                    '@subsection This is the first subsection.\n'
                    '\n'
                    'This is the second paragraph.\n'
                    '\n'
                    '@see Other\n'
                    '\n')
        self.assertMultiLineEqual(doc.getFormatted(50), expected)

    def testInterfaceMinimal(self):
        txt = '@mfn Klass#Metafunktion'
        parser = self.parseText(txt)
        doc = parser.documentation
        self.assertMultiLineEqual(doc.getFormatted(), '@mfn Klass#Metafunktion\n\n')

    def testInterfaceFull(self):
        txt = ('@mfn Klass#Metafunktion\n'
               '@brief This is a brief text.\n'
               '@deprecated Deprecation message.\n'
               '@signature Metafunktion<T1>::Type;\n'
               '@signature Metafunktion<T2>::VALUE;\n'
               '@tparam T1 The first value and a very very very very long \n'
               '           description\n'
               '@tparam T2 The second value and a very very very very long \n'
               '           description\n'
               '@return Type The return type.\n'
               '@return VALUE The return value.\n'
               '\n'
               '@section This is the first section.\n'
               '\n'
               'This is the first paragraph.\n'
               '\n'
               '@subsection This is the first subsection.\n'
               '\n'
               'This is the second paragraph.\n'
               '@see Other')
        parser = self.parseText(txt)
        doc = parser.documentation
        expected = ('@mfn Klass#Metafunktion\n'
                    '\n'
                    '@brief This is a brief text.\n'
                    '\n'
                    '@deprecated Deprecation message.\n'
                    '\n'
                    '@signature Metafunktion<T1>::Type;\n'
                    '@signature Metafunktion<T2>::VALUE;\n'
                    '\n'
                    '@tparam T1 The first value and a very very very\n'
                    '           very long description\n'
                    '@tparam T2 The second value and a very very very\n'
                    '           very long description\n'
                    '\n'
                    '@return Type The return type.\n'
                    '@return VALUE The return value.\n'
                    '\n'
                    '@section This is the first section.\n'
                    '\n'
                    'This is the first paragraph.\n'
                    '\n'
                    '@subsection This is the first subsection.\n'
                    '\n'
                    'This is the second paragraph.\n'
                    '\n'
                    '@see Other\n'
                    '\n')
        self.assertMultiLineEqual(doc.getFormatted(50), expected)


class TestConceptParsing(TestDoxParserBase):
    def testMinimal(self):
        txt = '@concept Konzept'
        parser = self.parseText(txt)
        doc = parser.documentation
        self.assertMultiLineEqual(doc.getFormatted(), '@concept Konzept\n\n')

    def testFull(self):
        txt = ('@concept Konzept\n'
               '@brief This is a brief text.\n'
               '@deprecated Deprecation message.\n'
               '@headerfile <seqan/base.h>\n'
               '@headerfile <seqan/sequence.h>\n'
               '@signature concept Konzept;\n'
               '\n'
               '@section This is the first section.\n'
               '\n'
               'This is the first paragraph.\n'
               '\n'
               '@subsection This is the first subsection.\n'
               '\n'
               'This is the second paragraph.\n'
               '@see Other')
        parser = self.parseText(txt)
        doc = parser.documentation
        expected = ('@concept Konzept\n'
                    '\n'
                    '@headerfile <seqan/base.h>\n'
                    '@headerfile <seqan/sequence.h>\n'
                    '\n'
                    '@brief This is a brief text.\n'
                    '\n'
                    '@deprecated Deprecation message.\n'
                    '\n'
                    '@signature concept Konzept;\n'
                    '\n'
                    '@section This is the first section.\n'
                    '\n'
                    'This is the first paragraph.\n'
                    '\n'
                    '@subsection This is the first subsection.\n'
                    '\n'
                    'This is the second paragraph.\n'
                    '\n'
                    '@see Other\n'
                    '\n')
        self.assertMultiLineEqual(doc.getFormatted(50), expected)


class TestVariableParsing(TestDoxParserBase):
    def testFull(self):
        txt = ('@var int var\n'
               '@brief This is a brief text.\n'
               '@deprecated Deprecation message.\n'
               '@headerfile <seqan/base.h>\n'
               '@headerfile <seqan/sequence.h>\n'
               '@signature int var;\n'
               '\n'
               '@section This is the first section.\n'
               '\n'
               'This is the first paragraph.\n'
               '\n'
               '@subsection This is the first subsection.\n'
               '\n'
               'This is the second paragraph.\n'
               '@see Other')
        parser = self.parseText(txt)
        doc = parser.documentation
        expected = ('@var int var;\n'
                    '\n'
                    '@headerfile <seqan/base.h>\n'
                    '@headerfile <seqan/sequence.h>\n'
                    '\n'
                    '@brief This is a brief text.\n'
                    '\n'
                    '@deprecated Deprecation message.\n'
                    '\n'
                    '@signature int var;\n'
                    '\n'
                    '@section This is the first section.\n'
                    '\n'
                    'This is the first paragraph.\n'
                    '\n'
                    '@subsection This is the first subsection.\n'
                    '\n'
                    'This is the second paragraph.\n'
                    '\n'
                    '@see Other\n'
                    '\n')
        self.assertMultiLineEqual(doc.getFormatted(50), expected)

    def testMemberMinimal(self):
        txt = '@var Klass::Type Klass::var'
        parser = self.parseText(txt)
        doc = parser.documentation
        self.assertMultiLineEqual(doc.getFormatted(), '@var Klass::Type Klass::var;\n\n')

    def testMemberFull(self):
        txt = ('@var Klass::Type Klass::var\n'
               '@brief This is a brief text.\n'
               '@deprecated Deprecation message.\n'
               '@signature Klass::Type Klass::var;\n'
               '\n'
               '@section This is the first section.\n'
               '\n'
               'This is the first paragraph.\n'
               '\n'
               '@subsection This is the first subsection.\n'
               '\n'
               'This is the second paragraph.\n'
               '@see Other')
        parser = self.parseText(txt)
        doc = parser.documentation
        expected = ('@var Klass::Type Klass::var;\n'
                    '\n'
                    '@brief This is a brief text.\n'
                    '\n'
                    '@deprecated Deprecation message.\n'
                    '\n'
                    '@signature Klass::Type Klass::var;\n'
                    '\n'
                    '@section This is the first section.\n'
                    '\n'
                    'This is the first paragraph.\n'
                    '\n'
                    '@subsection This is the first subsection.\n'
                    '\n'
                    'This is the second paragraph.\n'
                    '\n'
                    '@see Other\n'
                    '\n')
        self.assertMultiLineEqual(doc.getFormatted(50), expected)


class TestTagParsing(TestDoxParserBase):
    def testFull(self):
        txt = ('@tag TagName\n'
               '@brief This is a brief text.\n'
               '@deprecated Deprecation message.\n'
               '@headerfile <seqan/base.h>\n'
               '@headerfile <seqan/sequence.h>\n'
               '@signature typedef Tag<TagName_> TagName;\n'
               '\n'
               '@section This is the first section.\n'
               '\n'
               'This is the first paragraph.\n'
               '\n'
               '@subsection This is the first subsection.\n'
               '\n'
               'This is the second paragraph.\n'
               '@see Other')
        parser = self.parseText(txt)
        doc = parser.documentation
        expected = ('@tag TagName\n'
                    '\n'
                    '@headerfile <seqan/base.h>\n'
                    '@headerfile <seqan/sequence.h>\n'
                    '\n'
                    '@brief This is a brief text.\n'
                    '\n'
                    '@deprecated Deprecation message.\n'
                    '\n'
                    '@signature typedef Tag<TagName_> TagName;\n'
                    '\n'
                    '@section This is the first section.\n'
                    '\n'
                    'This is the first paragraph.\n'
                    '\n'
                    '@subsection This is the first subsection.\n'
                    '\n'
                    'This is the second paragraph.\n'
                    '\n'
                    '@see Other\n'
                    '\n')
        self.assertMultiLineEqual(doc.getFormatted(50), expected)


class TestEnumParsing(TestDoxParserBase):
    def testMinimal(self):
        txt = '@enum Enum'
        parser = self.parseText(txt)
        doc = parser.documentation
        self.assertMultiLineEqual(doc.getFormatted(), '@enum Enum\n\n')

    def testFull(self):
        txt = ('@enum Enum\n'
               '@brief This is a brief text.\n'
               '@deprecated Deprecation message.\n'
               '@headerfile <seqan/base.h>\n'
               '@headerfile <seqan/sequence.h>\n'
               '@signature enum Enum;\n'
               '\n'
               '@section This is the first section.\n'
               '\n'
               'This is the first paragraph.\n'
               '\n'
               '@subsection This is the first subsection.\n'
               '\n'
               'This is the second paragraph.\n'
               '@see Other')
        parser = self.parseText(txt)
        doc = parser.documentation
        expected = ('@enum Enum\n'
                    '\n'
                    '@headerfile <seqan/base.h>\n'
                    '@headerfile <seqan/sequence.h>\n'
                    '\n'
                    '@brief This is a brief text.\n'
                    '\n'
                    '@deprecated Deprecation message.\n'
                    '\n'
                    '@signature enum Enum;\n'
                    '\n'
                    '@section This is the first section.\n'
                    '\n'
                    'This is the first paragraph.\n'
                    '\n'
                    '@subsection This is the first subsection.\n'
                    '\n'
                    'This is the second paragraph.\n'
                    '\n'
                    '@see Other\n'
                    '\n')
        self.assertMultiLineEqual(doc.getFormatted(50), expected)


class TestPageParsingWithInclude(TestDoxParserBase):
    """Tests for the @include and @snippet command.

    We use a @page for simplicity.
    """

    def testInclude(self):
        txt = ('@page Page Title\n'
               '@include example.cpp')
        parser = self.parseText(txt)
        doc = parser.documentation
        txt = ('@page Page Title\n'
               '\n'
               '@include example.cpp\n\n')
        self.assertMultiLineEqual(doc.getFormatted(), txt)

    def testSnippet(self):
        txt = ('@page Page Title\n'
               '@snippet example.cpp Snippet Name')
        parser = self.parseText(txt)
        doc = parser.documentation
        txt = ('@page Page Title\n'
               '\n'
               '@snippet example.cpp Snippet Name\n\n')
        self.assertMultiLineEqual(doc.getFormatted(), txt)

        
class TestPageParsing(TestDoxParserBase):
    def testMinimal(self):
        txt = '@page Page Title'
        parser = self.parseText(txt)
        doc = parser.documentation
        self.assertMultiLineEqual(doc.getFormatted(), '@page Page Title\n\n')

    def testSmallBody(self):
        txt = ('@page Page Title\n'
               'This is the body.')
        parser = self.parseText(txt)
        doc = parser.documentation
        txt = ('@page Page Title\n'
               '\n'
               'This is the body.\n\n')
        self.assertMultiLineEqual(doc.getFormatted(), txt)

    def testFull(self):
        txt = ('@page Page Title\n'
               '@section This is the first section.\n'
               '\n'
               'This is the first paragraph.\n'
               '@include path/to/file#hash\n'
               'This is another paragraph.\n'
               '\n'
               '@subsection This is the first subsection.\n'
               '\n'
               'This is the second paragraph.\n'
               '@see Other')
        parser = self.parseText(txt)
        doc = parser.documentation
        expected = ('@page Page Title\n'
                    '\n'
                    '@section This is the first section.\n'
                    '\n'
                    'This is the first paragraph.\n'
                    '\n'
                    '@include path/to/file#hash\n'
                    '\n'
                    'This is another paragraph.\n'
                    '\n'
                    '@subsection This is the first subsection.\n'
                    '\n'
                    'This is the second paragraph.\n'
                    '\n'
                    '@see Other\n'
                    '\n')
        self.assertMultiLineEqual(doc.getFormatted(50), expected)

        
class TestPageParsing(TestDoxParserBase):
    def testMinimal(self):
        txt = '@defgroup GroupName Group Title'
        parser = self.parseText(txt)
        doc = parser.documentation
        self.assertMultiLineEqual(doc.getFormatted(), '@defgroup GroupName Group Title\n\n')

    def testSmallBody(self):
        txt = ('@defgroup GroupName Group Title\n'
               'This is the body.')
        parser = self.parseText(txt)
        doc = parser.documentation
        txt = ('@defgroup GroupName Group Title\n'
               '\n'
               'This is the body.\n\n')
        self.assertMultiLineEqual(doc.getFormatted(), txt)

    def testFull(self):
        txt = ('@defgroup GroupName Group Title\n'
               '@section This is the first section.\n'
               '\n'
               'This is the first paragraph.\n'
               '@include path/to/file#hash\n'
               'This is another paragraph.\n'
               '\n'
               '@subsection This is the first subsection.\n'
               '\n'
               'This is the second paragraph.\n'
               '@see Other')
        parser = self.parseText(txt)
        doc = parser.documentation
        expected = ('@defgroup GroupName Group Title\n'
                    '\n'
                    '@section This is the first section.\n'
                    '\n'
                    'This is the first paragraph.\n'
                    '\n'
                    '@include path/to/file#hash\n'
                    '\n'
                    'This is another paragraph.\n'
                    '\n'
                    '@subsection This is the first subsection.\n'
                    '\n'
                    'This is the second paragraph.\n'
                    '\n'
                    '@see Other\n'
                    '\n')
        self.assertMultiLineEqual(doc.getFormatted(50), expected)


class TestLinkParsing(TestDoxParserBase):
    def testMinimal(self):
        txt = ('@page PageTitle Page Title\n'
               '\n'
               '@link PageTitle the page title @endlink.\n')
        parser = self.parseText(txt)
        doc = parser.documentation
        expected = ('@page PageTitle Page Title\n'
                    '\n'
                    '@link PageTitle the page title @endlink.\n\n')
        self.assertMultiLineEqual(doc.getFormatted(), expected)
        

if __name__ == '__main__':
    unittest.main()
