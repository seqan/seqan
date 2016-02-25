#!/usr/bin/env python2
"""Tests for the raw_doc module."""

__author__ = 'Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>'

# TODO(holtgrew): Implement test for Class, Function, Metafunction, Macro, Page, Group, Tag.

import unittest

import seqan.dox.lexer as lexer
import seqan.dox.raw_doc as raw_doc


class TestDoxFormatter(unittest.TestCase):
    def setUp(self):
        self.fmt = raw_doc.DoxFormatter(50)
        self.txt = ('This is a quite long string that is used to determine '
                    'whether the formatter wraps correctly.\n')

    def testFormatParagraph(self):
        res = self.fmt.formatParagraph(self.txt)
        txt = ('This is a quite long string that is used to\n'
               'determine whether the formatter wraps correctly.\n')
        self.assertMultiLineEqual(res, txt)

    def testFormatCommandNoLeading(self):
        res = self.fmt.formatCommand('cmd', self.txt)
        txt = ('@cmd This is a quite long string that is used to\n'
               '     determine whether the formatter wraps\n'
               '     correctly.\n')
        self.assertMultiLineEqual(res, txt)

    def testFormatCommandLeading(self):
        res = self.fmt.formatCommand('cmd', self.txt, 'leading')
        txt = ('@cmd leading This is a quite long string that is\n'
               '             used to determine whether the\n'
               '             formatter wraps correctly.\n')
        self.assertMultiLineEqual(res, txt)


class TestText(unittest.TestCase):
    def testConstructionWithTokens(self):
        tokens = [lexer.Token('WORD', 'test', 0, 0, 0),
                  lexer.Token('SPACE', ' ', 0, 0, 0),
                  lexer.Token('WORD', 'foo', 0, 0, 0)]
        text = raw_doc.RawText(tokens)
        self.assertEqual(text.tokens, tokens)
        self.failIf(text.empty)

    def testConstructionWithoutTokens(self):
        text = raw_doc.RawText()
        self.assertEqual(text.tokens, [])
        self.assert_(text.empty)


class TestDocumentation(unittest.TestCase):
    # TODO(holtgrew): Write more tests!

    def testMerge(self):
        self.page_token = lexer.Token('COMMAND_PAGE', '@page', 0, 0, 0)
        doc_left = raw_doc.RawDoc()
        page_left = raw_doc.RawPage(self.page_token)
        doc_left.entries.append(page_left)
        doc_right = raw_doc.RawDoc()
        page_right = raw_doc.RawPage(self.page_token)
        doc_right.entries.append(page_right)
        doc_left.merge(doc_right)

        self.assertEqual(len(doc_left.entries), 2)
        self.assertEqual(len(doc_right.entries), 1)
        self.assertEqual(doc_left.entries[0], page_left)
        self.assertEqual(doc_left.entries[1], page_right)


class TestEntry(unittest.TestCase):
    def setUp(self):
        self.tok = lexer.Token('WORD', 'entry', 0, 0, 0)
        self.brief_tok = lexer.Token('WORD', 'This is brief.', 0, 0, 0)
        self.name_tok = lexer.Token('WORD', 'Concept', 0, 0, 0)
        self.title_tok = lexer.Token('WORD', 'Concept Title', 0, 0, 0)
        self.tok_see = lexer.Token('WORD', 'See', 0, 0, 0)

    def testInitialization(self):
        entry = raw_doc.RawEntry(self.tok)
        self.assertEqual(entry.name.text, raw_doc.RawText().text)
        self.assertEqual(len(entry.briefs), 0)
        self.assertEqual(entry.body, raw_doc.RawBody())
        self.assertEqual(entry.sees, [])
        self.assertEqual(entry.command, '<entry>')

    def testAddBrief(self):
        entry = raw_doc.RawEntry(self.tok)
        self.assertEqual(entry.briefs, [])
        b = raw_doc.RawBrief(self.brief_tok, raw_doc.RawText())
        entry.addBrief(b)
        self.assertEqual(len(entry.briefs), 1)
        self.assertEqual(entry.briefs, [b])

    def testAddSee(self):
        entry = raw_doc.RawEntry(self.tok_see)
        self.assertEqual(entry.sees, [])
        s = raw_doc.RawSee(self.tok_see, raw_doc.RawText())
        entry.addSee(s)
        self.assertEqual(len(entry.sees), 1)
        self.assertEqual(entry.sees, [s])

    def testEntryTypes(self):
        expected = ('concept', 'class', 'function', 'metafunction', 'page',
                    'enum', 'var', 'tag', 'defgroup', 'macro', 'enum_value')
        self.assertEqual(raw_doc.RawEntry.entryTypes(), expected)

    def testAddParagraph(self):
        entry = raw_doc.RawEntry(self.tok)
        self.assertEqual(entry.body, raw_doc.RawBody())
        p = raw_doc.RawParagraph(self.tok)
        entry.addParagraph(p)
        b = raw_doc.RawBody()
        b.addParagraph(p)
        self.assertEqual(entry.body, b)

    def testGetFormatted(self):
        b = raw_doc.RawBrief(self.brief_tok, raw_doc.RawText([self.brief_tok]))
        entry = raw_doc.RawEntry(self.brief_tok, [b])
        entry.name = raw_doc.RawText([self.name_tok])
        entry.title = raw_doc.RawText([self.title_tok])
        entry.sees = [raw_doc.RawSee(self.tok_see, raw_doc.RawText([self.tok_see]))]
        formatter = raw_doc.DoxFormatter()
        msg = ('@<entry> Concept Concept Title\n\n'
               '@brief This is brief.\n\n'
               '@see See\n\n')
        self.assertMultiLineEqual(entry.getFormatted(formatter), msg)


class CodeEntryTest(unittest.TestCase):
    def setUp(self):
        self.code_tok = lexer.Token('COMMAND_CODE', '@code', 0, 0, 0)
        self.sig_tok = lexer.Token('COMMAND_SIGNATURE', '@signature', 0, 0, 0)
        self.brief_tok = lexer.Token('WORD', 'This is brief.', 0, 0, 0)
        self.name_tok = lexer.Token('WORD', 'Concept', 0, 0, 0)
        self.title_tok = lexer.Token('WORD', 'Concept Title', 0, 0, 0)
        self.tok_see = lexer.Token('WORD', 'See', 0, 0, 0)
        self.tok_sig = lexer.Token('WORD', 'payload', 0, 0, 0)

    def testCreation(self):
        code_entry = raw_doc.RawCodeEntry(self.code_tok)
        self.assertEqual(code_entry.signatures, [])
        self.assertEqual(code_entry.command, '<code entry>')

    def testAddSignature(self):
        code_entry = raw_doc.RawCodeEntry(self.code_tok)
        s = raw_doc.RawSignature(self.sig_tok, raw_doc.RawText(
            [lexer.Token('WORD', 'payload', 0, 0, 0)]))
        code_entry.addSignature(s)

    def testGetFormatted(self):
        code_entry = raw_doc.RawCodeEntry(self.code_tok)
        b = raw_doc.RawBrief(self.brief_tok, raw_doc.RawText([self.brief_tok]))
        code_entry = raw_doc.RawCodeEntry(self.code_tok, [b])
        code_entry.name = raw_doc.RawText([self.name_tok])
        code_entry.title = raw_doc.RawText([self.title_tok])
        code_entry.sees = [raw_doc.RawSee(self.tok_see, raw_doc.RawText([self.tok_see]))]
        s = raw_doc.RawSignature(self.tok_sig, raw_doc.RawText([self.tok_sig]))
        code_entry.addSignature(s)
        formatter = raw_doc.DoxFormatter()
        txt = ('@<code entry> Concept Concept Title\n\n'
               '@brief This is brief.\n\n'
               '@signature payload\n\n'
               '@see See\n\n')
        self.assertMultiLineEqual(code_entry.getFormatted(formatter), txt)


class ConceptTest(unittest.TestCase):
    def setUp(self):
        self.brief_tok = lexer.Token('WORD', 'This is brief.', 0, 0, 0)
        self.name_tok = lexer.Token('WORD', 'Concept', 0, 0, 0)
        self.title_tok = lexer.Token('WORD', 'Concept Title', 0, 0, 0)
        self.tok_see = lexer.Token('WORD', 'See', 0, 0, 0)
        self.tok_sig = lexer.Token('WORD', 'payload', 0, 0, 0)
        self.formatter = raw_doc.DoxFormatter()

    def testGetFormatted(self):
        b = raw_doc.RawBrief(self.brief_tok, raw_doc.RawText([self.brief_tok]))
        code_entry = raw_doc.RawConcept(self.brief_tok, [b])
        code_entry.name = raw_doc.RawText([self.name_tok])
        code_entry.title = raw_doc.RawText([self.title_tok])
        code_entry.sees = [raw_doc.RawSee(self.tok_see, raw_doc.RawText([self.tok_see]))]
        s = raw_doc.RawSignature(self.tok_sig, raw_doc.RawText([self.tok_sig]))
        code_entry.addSignature(s)
        txt = ('@concept Concept Concept Title\n\n'
               '@brief This is brief.\n\n'
               '@signature payload\n\n'
               '@see See\n\n')
        self.assertMultiLineEqual(code_entry.getFormatted(self.formatter), txt)


class EnumTest(unittest.TestCase):
    def setUp(self):
        self.brief_tok = lexer.Token('WORD', 'This is brief.', 0, 0, 0)
        self.name_tok = lexer.Token('WORD', 'Enum', 0, 0, 0)
        self.title_tok = lexer.Token('WORD', 'Enum Title', 0, 0, 0)
        self.tok_see = lexer.Token('WORD', 'See', 0, 0, 0)
        self.tok_sig = lexer.Token('WORD', 'payload', 0, 0, 0)
        self.formatter = raw_doc.DoxFormatter()

    def testGetFormatted(self):
        b = raw_doc.RawBrief(self.brief_tok, raw_doc.RawText([self.brief_tok]))
        code_entry = raw_doc.RawEnum(self.brief_tok, [b])
        code_entry.name = raw_doc.RawText([self.name_tok])
        code_entry.title = raw_doc.RawText([self.title_tok])
        code_entry.sees = [raw_doc.RawSee(self.tok_see, raw_doc.RawText([self.tok_see]))]
        s = raw_doc.RawSignature(self.tok_sig, raw_doc.RawText([self.tok_sig]))
        code_entry.addSignature(s)
        txt = ('@enum Enum Enum Title\n\n'
               '@brief This is brief.\n\n'
               '@signature payload\n\n'
               '@see See\n\n')
        self.assertMultiLineEqual(code_entry.getFormatted(self.formatter), txt)


class TypedefTest(unittest.TestCase):
    def setUp(self):
        self.brief_tok = lexer.Token('WORD', 'This is brief.', 0, 0, 0)
        self.name_tok = lexer.Token('WORD', 'TypeDef', 0, 0, 0)
        self.title_tok = lexer.Token('WORD', 'Typedef Title', 0, 0, 0)
        self.tok_see = lexer.Token('WORD', 'See', 0, 0, 0)
        self.tok_sig = lexer.Token('WORD', 'payload', 0, 0, 0)
        self.formatter = raw_doc.DoxFormatter()

    def testGetFormatted(self):
        b = raw_doc.RawBrief(self.brief_tok, raw_doc.RawText([self.brief_tok]))
        code_entry = raw_doc.RawTypedef(self.brief_tok, [b])
        code_entry.name = raw_doc.RawText([self.name_tok])
        code_entry.title = raw_doc.RawText([self.title_tok])
        code_entry.sees = [raw_doc.RawSee(self.tok_see, raw_doc.RawText([self.tok_see]))]
        s = raw_doc.RawSignature(self.tok_sig, raw_doc.RawText([self.tok_sig]))
        code_entry.addSignature(s)
        txt = ('@typedef TypeDef Typedef Title\n\n'
               '@brief This is brief.\n\n'
               '@signature payload\n\n'
               '@see See\n\n')
        self.assertMultiLineEqual(code_entry.getFormatted(self.formatter), txt)


class AdaptionTest(unittest.TestCase):
    def setUp(self):
        self.brief_tok = lexer.Token('WORD', 'This is brief.', 0, 0, 0)
        self.name_tok = lexer.Token('WORD', 'Adaption', 0, 0, 0)
        self.title_tok = lexer.Token('WORD', 'Adaption Title', 0, 0, 0)
        self.tok_see = lexer.Token('WORD', 'See', 0, 0, 0)
        self.tok_sig = lexer.Token('WORD', 'payload', 0, 0, 0)
        self.formatter = raw_doc.DoxFormatter()

    def testGetFormatted(self):
        b = raw_doc.RawBrief(self.brief_tok, raw_doc.RawText([self.brief_tok]))
        code_entry = raw_doc.RawAdaption(self.brief_tok, [b])
        code_entry.name = raw_doc.RawText([self.name_tok])
        code_entry.title = raw_doc.RawText([self.title_tok])
        code_entry.sees = [raw_doc.RawSee(self.tok_see, raw_doc.RawText([self.tok_see]))]
        s = raw_doc.RawSignature(self.tok_sig, raw_doc.RawText([self.tok_sig]))
        code_entry.addSignature(s)
        txt = ('@adaption Adaption Adaption Title\n\n'
               '@brief This is brief.\n\n'
               '@signature payload\n\n'
               '@see See\n\n')
        self.assertMultiLineEqual(code_entry.getFormatted(self.formatter), txt)


class VariableTest(unittest.TestCase):
    def setUp(self):
        self.brief_tok = lexer.Token('WORD', 'This is brief.', 0, 0, 0)
        self.name_tok = lexer.Token('WORD', 'var', 0, 0, 0)
        self.type_tok = lexer.Token('WORD', 'int', 0, 0, 0)
        self.tok_see = lexer.Token('WORD', 'See', 0, 0, 0)
        self.tok_sig = lexer.Token('WORD', 'payload', 0, 0, 0)
        self.formatter = raw_doc.DoxFormatter()

    def testGetFormatted(self):
        b = raw_doc.RawBrief(self.brief_tok, raw_doc.RawText([self.brief_tok]))
        code_entry = raw_doc.RawVariable(self.brief_tok, [b])
        code_entry.name = raw_doc.RawText([self.name_tok])
        code_entry.type = raw_doc.RawText([self.type_tok])
        code_entry.sees = [raw_doc.RawSee(self.tok_see, raw_doc.RawText([self.tok_see]))]
        s = raw_doc.RawSignature(self.tok_sig, raw_doc.RawText([self.tok_sig]))
        code_entry.addSignature(s)
        txt = ('@var int var;\n\n'
               '@brief This is brief.\n\n'
               '@signature payload\n\n'
               '@see See\n\n')
        self.assertMultiLineEqual(code_entry.getFormatted(self.formatter), txt)


class BodyTest(unittest.TestCase):
    def setUp(self):
        self.t = lexer.Token('WORD', 'aword', 0, 0, 0)
        self.p = raw_doc.RawParagraph(self.t, raw_doc.RawText([self.t]))
        self.formatter = raw_doc.DoxFormatter()

    def testCreation(self):
        body = raw_doc.RawBody()
        self.assertEqual(body.paragraphs, [])
        body.addParagraph(self.p)
        self.assertEqual(body.paragraphs, [self.p])

    def testGetFormatted(self):
        body = raw_doc.RawBody()
        body.addParagraph(self.p)
        self.assertEqual(body.getFormatted(self.formatter), 'aword\n')


class SectionTest(unittest.TestCase):
    def setUp(self):
        self.t = lexer.Token('WORD', 'aword', 0, 0, 0)
        self.formatter = raw_doc.DoxFormatter()
        self.txt = raw_doc.RawText([self.t])

    def testCreation(self):
        section = raw_doc.RawSection(self.t, self.txt)
        self.assertEqual(section.heading, self.txt)
        self.assertEqual(section.level, 0)

    def testGetType(self):
        section = raw_doc.RawSection(self.t, self.txt, 1)
        self.assertEqual(section.getType(), 'section')

    def testCreationWithLevel(self):
        section = raw_doc.RawSection(self.t, self.txt, 1)
        self.assertEqual(section.heading, self.txt)
        self.assertEqual(section.level, 1)

    def testGetCommand(self):
        section = raw_doc.RawSection(self.t, self.txt, 0)
        self.assertEqual(section.getCommand(), 'section')
        section = raw_doc.RawSection(self.t, self.txt, 1)
        self.assertEqual(section.getCommand(), 'subsection')

    def testGetFormatted(self):
        section = raw_doc.RawSection(self.t, self.txt, 1)
        self.assertEqual(section.getFormatted(self.formatter), '@subsection aword\n')


class IncludeTest(unittest.TestCase):
    def setUp(self):
        self.path_t = lexer.Token('WORD', 'apath', 0, 0, 0)
        self.path = raw_doc.RawText([self.path_t])
        self.formatter = raw_doc.DoxFormatter()

    def testCreation(self):
        include = raw_doc.RawInclude(self.path_t, [self.path_t])
        self.assertEqual(include.path, self.path)
        self.assertEqual(include.getFormatted(self.formatter), '@include apath\n')


class SnippetTest(unittest.TestCase):
    def setUp(self):
        self.path_t = lexer.Token('WORD', 'apath', 0, 0, 0)
        self.path = raw_doc.RawText([self.path_t])
        self.snippet_t0 = lexer.Token('WORD', 'The', 0, 0, 0)
        self.snippet_t1 = lexer.Token('SPACE', ' ', 0, 0, 0)
        self.snippet_t2 = lexer.Token('WORD', 'snippet', 0, 0, 0)
        self.snippet = raw_doc.RawText([self.snippet_t0, self.snippet_t1,
                                        self.snippet_t2])
        self.formatter = raw_doc.DoxFormatter()

    def testCreation(self):
        snippet = raw_doc.RawSnippet(self.path_t, self.path.tokens, self.snippet.tokens)
        self.assertEqual(snippet.getFormatted(self.formatter),
                         '@snippet apath The snippet\n')


class CodeTest(unittest.TestCase):
    def setUp(self):
        self.t = lexer.Token('WORD', 'aword', 0, 0, 0)
        self.txt = raw_doc.RawText([self.t])
        self.formatter = raw_doc.DoxFormatter()

    def testCreation(self):
        paragraph = raw_doc.RawCode(self.t, self.txt)
        self.assertEqual(paragraph.text, self.txt)
        self.assertEqual(paragraph.extension, '.txt')

    def testCreationWithExtension(self):
        paragraph = raw_doc.RawCode(self.t, self.txt, '.cpp')
        self.assertEqual(paragraph.text, self.txt)
        self.assertEqual(paragraph.extension, '.cpp')

    def testGetType(self):
        pass

    def testGetFormatted(self):
        pass


class BriefTest(unittest.TestCase):
    def setUp(self):
        self.t = lexer.Token('WORD', 'aword', 0, 0, 0)
        self.txt = raw_doc.RawText([self.t])
        self.formatter = raw_doc.DoxFormatter()

    def testCreation(self):
        brief = raw_doc.RawBrief(self.t, self.txt)
        self.assertEqual(brief.text, self.txt)

    def testGetType(self):
        brief = raw_doc.RawBrief(self.t, self.txt)
        self.assertEqual(brief.getType(), 'brief')

    def testGetFormatted(self):
        brief = raw_doc.RawBrief(self.t, self.txt)
        self.assertEqual(brief.getFormatted(self.formatter), '@brief aword\n')


class ExtendsTest(unittest.TestCase):
    def setUp(self):
        self.t = lexer.Token('WORD', 'aword', 0, 0, 0)
        self.txt = raw_doc.RawText([self.t])
        self.formatter = raw_doc.DoxFormatter()

    def testCreation(self):
        extends = raw_doc.RawExtends(self.t, self.txt)
        self.assertEqual(extends.text, self.txt)

    def testGetType(self):
        extends = raw_doc.RawExtends(self.t, self.txt)
        self.assertEqual(extends.getType(), 'extends')

    def testGetFormatted(self):
        extends = raw_doc.RawExtends(self.t, self.txt)
        self.assertEqual(extends.getFormatted(self.formatter), '@extends aword\n')


class ImplementsTest(unittest.TestCase):
    def setUp(self):
        self.t = lexer.Token('WORD', 'aword', 0, 0, 0)
        self.txt = raw_doc.RawText([self.t])
        self.formatter = raw_doc.DoxFormatter()

    def testCreation(self):
        implements = raw_doc.RawImplements(self.t, self.txt)
        self.assertEqual(implements.text, self.txt)

    def testGetType(self):
        implements = raw_doc.RawImplements(self.t, self.txt)
        self.assertEqual(implements.getType(), 'implements')

    def testGetFormatted(self):
        implements = raw_doc.RawImplements(self.t, self.txt)
        self.assertEqual(implements.getFormatted(self.formatter), '@implements aword\n')


class SeeTest(unittest.TestCase):
    def setUp(self):
        self.t = lexer.Token('WORD', 'aword', 0, 0, 0)
        self.txt = raw_doc.RawText([self.t])
        self.formatter = raw_doc.DoxFormatter()

    def testCreation(self):
        see = raw_doc.RawSee(self.t, self.txt)
        self.assertEqual(see.text, self.txt)

    def testGetType(self):
        see = raw_doc.RawSee(self.t, self.txt)
        self.assertEqual(see.getType(), 'see')

    def testGetFormatted(self):
        see = raw_doc.RawSee(self.t, self.txt)
        self.assertEqual(see.getFormatted(self.formatter), '@see aword\n')


class ParamTest(unittest.TestCase):
    def setUp(self):
        self.tok = lexer.Token('COMMAND_PARAM', '@param', 0, 0, 0)
        self.tok_name = lexer.Token('WORD', 'name', 0, 0, 0)
        self.tok_text = lexer.Token('WORD', 'text', 0, 0, 0)
        self.tok_inout = lexer.Token('PARAM_IN_OUT', '[in,out]', 0, 0, 0)
        self.txt_name = raw_doc.RawText([self.tok_name])
        self.txt_text = raw_doc.RawText([self.tok_text])
        self.formatter = raw_doc.DoxFormatter()

    def testCreation(self):
        param = raw_doc.RawParam(self.tok, self.txt_name, self.txt_text)
        self.assertEqual(param.name, self.txt_name)
        self.assertEqual(param.text, self.txt_text)
        self.assertEqual(param.inout, None)

    def testCreationInOut(self):
        param = raw_doc.RawParam(self.tok, self.txt_name, self.txt_text, self.tok_inout)
        self.assertEqual(param.name, self.txt_name)
        self.assertEqual(param.text, self.txt_text)
        self.assertEqual(param.inout, self.tok_inout)

    def testGetType(self):
        param = raw_doc.RawParam(self.tok, self.txt_name, self.txt_text)
        self.assertEqual(param.getType(), 'param')

    def testGetFormatted(self):
        param = raw_doc.RawParam(self.tok, self.txt_name, self.txt_text)
        self.assertEqual(param.getFormatted(self.formatter),
                         '@param name text\n')
        param = raw_doc.RawParam(self.tok, self.txt_name, self.txt_text, self.tok_inout)
        self.assertEqual(param.getFormatted(self.formatter),
                         '@param[in,out] name text\n')


class TParamTest(unittest.TestCase):
    def setUp(self):
        self.tok = lexer.Token('COMMAND_TPARAM', '@tparam', 0, 0, 0)
        self.tok_name = lexer.Token('WORD', 'name', 0, 0, 0)
        self.tok_text = lexer.Token('WORD', 'text', 0, 0, 0)
        self.txt_name = raw_doc.RawText([self.tok_name])
        self.txt_text = raw_doc.RawText([self.tok_text])
        self.formatter = raw_doc.DoxFormatter()

    def testCreation(self):
        param = raw_doc.RawTParam(self.tok, self.txt_name, self.txt_text)
        self.assertEqual(param.name, self.txt_name)
        self.assertEqual(param.text, self.txt_text)

    def testGetType(self):
        param = raw_doc.RawTParam(self.tok, self.txt_name, self.txt_text)
        self.assertEqual(param.getType(), 'tparam')

    def testGetFormatted(self):
        param = raw_doc.RawTParam(self.tok, self.txt_name, self.txt_text)
        self.assertEqual(param.getFormatted(self.formatter),
                         '@tparam name text\n')


class ReturnTest(unittest.TestCase):
    def setUp(self):
        self.tok = lexer.Token('COMMAND_RETURN', '@return', 0, 0, 0)
        self.tok_name = lexer.Token('WORD', 'name', 0, 0, 0)
        self.tok_text = lexer.Token('WORD', 'text', 0, 0, 0)
        self.txt_name = raw_doc.RawText([self.tok_name])
        self.txt_text = raw_doc.RawText([self.tok_text])
        self.formatter = raw_doc.DoxFormatter()

    def testCreation(self):
        ret = raw_doc.RawReturn(self.tok, self.txt_name, self.txt_text)
        self.assertEqual(ret.name, self.txt_name)
        self.assertEqual(ret.text, self.txt_text)
        self.assertEqual(ret.inout, None)

    def testGetType(self):
        ret = raw_doc.RawReturn(self.tok, self.txt_name, self.txt_text)
        self.assertEqual(ret.getType(), 'return')

    def testGetFormatted(self):
        ret = raw_doc.RawReturn(self.tok, self.txt_name, self.txt_text)
        self.assertEqual(ret.getFormatted(self.formatter),
                         '@return name text\n')


class ThrowTest(unittest.TestCase):
    """Test for the RawThrow class, mostly tests instance variables."""

    def setUp(self):
        self.tok = lexer.Token('COMMAND_THROW', '@throw', 0, 0, 0)
        self.tok_type = lexer.Token('WORD', 'type', 0, 0, 0)
        self.tok_text = lexer.Token('WORD', 'text', 0, 0, 0)
        self.txt_type = raw_doc.RawText([self.tok_type])
        self.txt_text = raw_doc.RawText([self.tok_text])
        self.formatter = raw_doc.DoxFormatter()

    def testCreation(self):
        ret = raw_doc.RawThrow(self.tok, self.txt_type, self.txt_text)
        self.assertEqual(ret.name, self.txt_type)
        self.assertEqual(ret.text, self.txt_text)
        self.assertEqual(ret.inout, None)

    def testGetType(self):
        ret = raw_doc.RawThrow(self.tok, self.txt_type, self.txt_text)
        self.assertEqual(ret.getType(), 'throw')

    def testGetFormatted(self):
        ret = raw_doc.RawThrow(self.tok, self.txt_type, self.txt_text)
        self.assertEqual(ret.getFormatted(self.formatter),
                         '@throw type text\n')
        

class DataRaceTest(unittest.TestCase):
    """Test for the RawDataRace class, mostly tests instance variables."""

    def setUp(self):
        self.tok = lexer.Token('COMMAND_DATARACE', '@datarace', 0, 0, 0)
        self.tok_text = lexer.Token('WORD', 'text', 0, 0, 0)
        self.txt_text = raw_doc.RawText([self.tok_text])
        self.formatter = raw_doc.DoxFormatter()

    def testCreation(self):
        dataRace = raw_doc.RawDataRace(self.tok, self.txt_text)
        self.assertEqual(dataRace.text, self.txt_text)

    def testGetType(self):
        dataRace = raw_doc.RawDataRace(self.tok, self.txt_text)
        self.assertEqual(dataRace.getType(), 'datarace')

    def testGetFormatted(self):
        dataRace = raw_doc.RawDataRace(self.tok, self.txt_text)
        self.assertEqual(dataRace.getFormatted(self.formatter), '@datarace text\n')


class SignatureTest(unittest.TestCase):
    def setUp(self):
        self.tok = lexer.Token('COMMAND_SIGNATURE', '@signature', 0, 0, 0)
        self.tok_text = lexer.Token('WORD', 'text', 0, 0, 0)
        self.txt_text = raw_doc.RawText([self.tok_text])
        self.formatter = raw_doc.DoxFormatter()

    def testCreation(self):
        signature = raw_doc.RawSignature(self.tok, self.txt_text)
        self.assertEqual(signature.text, self.txt_text)

    def testGetType(self):
        signature = raw_doc.RawSignature(self.tok, self.txt_text)
        self.assertEqual(signature.getType(), 'signature')

    def testGetFormatted(self):
        signature = raw_doc.RawSignature(self.tok, self.txt_text)
        self.assertEqual(signature.getFormatted(self.formatter),
                         '@signature text\n')


if __name__ == '__main__':
    unittest.main()
