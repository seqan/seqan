#!/usr/bin/env python
"""Tests for the raw_doc module."""

__author__ = 'Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>'

# TODO(holtgrew): Implement test for Class, Function, Metafunction, Macro, Page, Group, Tag.

import unittest

import lexer
import raw_doc


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
        doc_left = raw_doc.RawDoc()
        page_left = raw_doc.RawPage()
        doc_left.entries.append(page_left)
        doc_right = raw_doc.RawDoc()
        page_right = raw_doc.RawPage()
        doc_right.entries.append(page_right)
        doc_left.merge(doc_right)

        self.assertEqual(len(doc_left.entries), 2)
        self.assertEqual(len(doc_right.entries), 1)
        self.assertEqual(doc_left.entries[0], page_left)
        self.assertEqual(doc_left.entries[1], page_right)


class TestEntry(unittest.TestCase):
    def setUp(self):
        self.brief_tok = lexer.Token('WORD', 'This is brief.', 0, 0, 0)
        self.name_tok = lexer.Token('WORD', 'Concept', 0, 0, 0)
        self.tok_see = lexer.Token('WORD', 'See', 0, 0, 0)

    def testInitialization(self):
        entry = raw_doc.RawEntry()
        self.assertEqual(entry.name.text, raw_doc.RawText().text)
        self.assertEqual(len(entry.briefs), 0)
        self.assertEqual(entry.body, raw_doc.RawBody())
        self.assertEqual(entry.sees, [])
        self.assertEqual(entry.command, '<entry>')

    def testAddBrief(self):
        entry = raw_doc.RawEntry()
        self.assertEqual(entry.briefs, [])
        b = raw_doc.RawBrief(raw_doc.RawText())
        entry.addBrief(b)
        self.assertEqual(len(entry.briefs), 1)
        self.assertEqual(entry.briefs, [b])

    def testAddSee(self):
        entry = raw_doc.RawEntry()
        self.assertEqual(entry.sees, [])
        s = raw_doc.RawSee(raw_doc.RawText())
        entry.addSee(s)
        self.assertEqual(len(entry.sees), 1)
        self.assertEqual(entry.sees, [s])

    def testEntryTypes(self):
        expected = ('concept', 'class', 'function', 'metafunction', 'page',
                    'enum', 'var', 'tag', 'defgroup', 'macro')
        self.assertEqual(raw_doc.RawEntry.entryTypes(), expected)

    def testAddParagraph(self):
        entry = raw_doc.RawEntry()
        self.assertEqual(entry.body, raw_doc.RawBody())
        p = raw_doc.RawParagraph()
        entry.addParagraph(p)
        b = raw_doc.RawBody()
        b.addParagraph(p)
        self.assertEqual(entry.body, b)

    def testGetFormatted(self):
        b = raw_doc.RawBrief(raw_doc.RawText([self.brief_tok]))
        entry = raw_doc.RawEntry([b])
        entry.name = raw_doc.RawText([self.name_tok])
        entry.sees = [raw_doc.RawSee(raw_doc.RawText([self.tok_see]))]
        formatter = raw_doc.DoxFormatter()
        msg = ('@<entry> Concept\n\n'
               '@brief This is brief.\n\n'
               '@see See\n\n')
        self.assertMultiLineEqual(entry.getFormatted(formatter), msg)


class CodeEntryTest(unittest.TestCase):
    def setUp(self):
        self.brief_tok = lexer.Token('WORD', 'This is brief.', 0, 0, 0)
        self.name_tok = lexer.Token('WORD', 'Concept', 0, 0, 0)
        self.tok_see = lexer.Token('WORD', 'See', 0, 0, 0)
        self.tok_sig = lexer.Token('WORD', 'payload', 0, 0, 0)

    def testCreation(self):
        code_entry = raw_doc.RawCodeEntry()
        self.assertEqual(code_entry.signatures, [])
        self.assertEqual(code_entry.command, '<code entry>')

    def testAddSignature(self):
        code_entry = raw_doc.RawCodeEntry()
        s = raw_doc.RawSignature(raw_doc.RawText([lexer.Token('WORD', 'payload', 0, 0, 0)]))
        code_entry.addSignature(s)

    def testGetFormatted(self):
        code_entry = raw_doc.RawCodeEntry()
        b = raw_doc.RawBrief(raw_doc.RawText([self.brief_tok]))
        code_entry = raw_doc.RawCodeEntry([b])
        code_entry.name = raw_doc.RawText([self.name_tok])
        code_entry.sees = [raw_doc.RawSee(raw_doc.RawText([self.tok_see]))]
        s = raw_doc.RawSignature(raw_doc.RawText([self.tok_sig]))
        code_entry.addSignature(s)
        formatter = raw_doc.DoxFormatter()
        txt = ('@<code entry> Concept\n\n'
               '@brief This is brief.\n\n'
               '@signature payload\n\n'
               '@see See\n\n')
        self.assertMultiLineEqual(code_entry.getFormatted(formatter), txt)


class ConceptTest(unittest.TestCase):
    def setUp(self):
        self.brief_tok = lexer.Token('WORD', 'This is brief.', 0, 0, 0)
        self.name_tok = lexer.Token('WORD', 'Concept', 0, 0, 0)
        self.tok_see = lexer.Token('WORD', 'See', 0, 0, 0)
        self.tok_sig = lexer.Token('WORD', 'payload', 0, 0, 0)
        self.formatter = raw_doc.DoxFormatter()

    def testGetFormatted(self):
        b = raw_doc.RawBrief(raw_doc.RawText([self.brief_tok]))
        code_entry = raw_doc.RawConcept([b])
        code_entry.name = raw_doc.RawText([self.name_tok])
        code_entry.sees = [raw_doc.RawSee(raw_doc.RawText([self.tok_see]))]
        s = raw_doc.RawSignature(raw_doc.RawText([self.tok_sig]))
        code_entry.addSignature(s)
        txt = ('@concept Concept\n\n'
               '@brief This is brief.\n\n'
               '@signature payload\n\n'
               '@see See\n\n')
        self.assertMultiLineEqual(code_entry.getFormatted(self.formatter), txt)


class EnumTest(unittest.TestCase):
    def setUp(self):
        self.brief_tok = lexer.Token('WORD', 'This is brief.', 0, 0, 0)
        self.name_tok = lexer.Token('WORD', 'Enum', 0, 0, 0)
        self.tok_see = lexer.Token('WORD', 'See', 0, 0, 0)
        self.tok_sig = lexer.Token('WORD', 'payload', 0, 0, 0)
        self.formatter = raw_doc.DoxFormatter()

    def testGetFormatted(self):
        b = raw_doc.RawBrief(raw_doc.RawText([self.brief_tok]))
        code_entry = raw_doc.RawEnum([b])
        code_entry.name = raw_doc.RawText([self.name_tok])
        code_entry.sees = [raw_doc.RawSee(raw_doc.RawText([self.tok_see]))]
        s = raw_doc.RawSignature(raw_doc.RawText([self.tok_sig]))
        code_entry.addSignature(s)
        txt = ('@enum Enum\n\n'
               '@brief This is brief.\n\n'
               '@signature payload\n\n'
               '@see See\n\n')
        self.assertMultiLineEqual(code_entry.getFormatted(self.formatter), txt)


class TypedefTest(unittest.TestCase):
    def setUp(self):
        self.brief_tok = lexer.Token('WORD', 'This is brief.', 0, 0, 0)
        self.name_tok = lexer.Token('WORD', 'TypeDef', 0, 0, 0)
        self.tok_see = lexer.Token('WORD', 'See', 0, 0, 0)
        self.tok_sig = lexer.Token('WORD', 'payload', 0, 0, 0)
        self.formatter = raw_doc.DoxFormatter()

    def testGetFormatted(self):
        b = raw_doc.RawBrief(raw_doc.RawText([self.brief_tok]))
        code_entry = raw_doc.RawTypedef([b])
        code_entry.name = raw_doc.RawText([self.name_tok])
        code_entry.sees = [raw_doc.RawSee(raw_doc.RawText([self.tok_see]))]
        s = raw_doc.RawSignature(raw_doc.RawText([self.tok_sig]))
        code_entry.addSignature(s)
        txt = ('@typedef TypeDef\n\n'
               '@brief This is brief.\n\n'
               '@signature payload\n\n'
               '@see See\n\n')
        self.assertMultiLineEqual(code_entry.getFormatted(self.formatter), txt)


class AdaptionTest(unittest.TestCase):
    def setUp(self):
        self.brief_tok = lexer.Token('WORD', 'This is brief.', 0, 0, 0)
        self.name_tok = lexer.Token('WORD', 'Adaption', 0, 0, 0)
        self.tok_see = lexer.Token('WORD', 'See', 0, 0, 0)
        self.tok_sig = lexer.Token('WORD', 'payload', 0, 0, 0)
        self.formatter = raw_doc.DoxFormatter()

    def testGetFormatted(self):
        b = raw_doc.RawBrief(raw_doc.RawText([self.brief_tok]))
        code_entry = raw_doc.RawAdaption([b])
        code_entry.name = raw_doc.RawText([self.name_tok])
        code_entry.sees = [raw_doc.RawSee(raw_doc.RawText([self.tok_see]))]
        s = raw_doc.RawSignature(raw_doc.RawText([self.tok_sig]))
        code_entry.addSignature(s)
        txt = ('@adaption Adaption\n\n'
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
        b = raw_doc.RawBrief(raw_doc.RawText([self.brief_tok]))
        code_entry = raw_doc.RawVariable([b])
        code_entry.name = raw_doc.RawText([self.name_tok])
        code_entry.type = raw_doc.RawText([self.type_tok])
        code_entry.sees = [raw_doc.RawSee(raw_doc.RawText([self.tok_see]))]
        s = raw_doc.RawSignature(raw_doc.RawText([self.tok_sig]))
        code_entry.addSignature(s)
        txt = ('@var int var\n\n'
               '@brief This is brief.\n\n'
               '@signature payload\n\n'
               '@see See\n\n')
        self.assertMultiLineEqual(code_entry.getFormatted(self.formatter), txt)


class BodyTest(unittest.TestCase):
    def setUp(self):
        self.t = lexer.Token('WORD', 'aword', 0, 0, 0)
        self.p = raw_doc.RawParagraph(raw_doc.RawText([self.t]))
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
        section = raw_doc.RawSection(self.txt)
        self.assertEqual(section.heading, self.txt)
        self.assertEqual(section.level, 0)

    def testGetType(self):
        section = raw_doc.RawSection(self.txt, 1)
        self.assertEqual(section.getType(), 'section')

    def testCreationWithLevel(self):
        section = raw_doc.RawSection(self.txt, 1)
        self.assertEqual(section.heading, self.txt)
        self.assertEqual(section.level, 1)

    def testGetCommand(self):
        section = raw_doc.RawSection(self.txt, 0)
        self.assertEqual(section.getCommand(), 'section')
        section = raw_doc.RawSection(self.txt, 1)
        self.assertEqual(section.getCommand(), 'subsection')

    def testGetFormatted(self):
        section = raw_doc.RawSection(self.txt, 1)
        self.assertEqual(section.getFormatted(self.formatter), '@subsection aword\n')


class IncludeTest(unittest.TestCase):
    def setUp(self):
        self.path_t = lexer.Token('WORD', 'apath', 0, 0, 0)
        self.path = raw_doc.RawText([self.path_t])
        self.formatter = raw_doc.DoxFormatter()

    def testCreation(self):
        include = raw_doc.RawInclude(self.path)
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
        snippet = raw_doc.RawSnippet(self.path, self.snippet)
        self.assertEqual(snippet.getFormatted(self.formatter),
                         '@snippet apath The snippet\n')


class CodeTest(unittest.TestCase):
    def setUp(self):
        self.t = lexer.Token('WORD', 'aword', 0, 0, 0)
        self.txt = raw_doc.RawText([self.t])
        self.formatter = raw_doc.DoxFormatter()

    def testCreation(self):
        paragraph = raw_doc.RawCode(self.txt)
        self.assertEqual(paragraph.text, self.txt)
        self.assertEqual(paragraph.extension, '.txt')

    def testCreationWithExtension(self):
        paragraph = raw_doc.RawCode(self.txt, '.cpp')
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
        brief = raw_doc.RawBrief(self.txt)
        self.assertEqual(brief.text, self.txt)

    def testGetType(self):
        brief = raw_doc.RawBrief(self.txt)
        self.assertEqual(brief.getType(), 'brief')

    def testGetFormatted(self):
        brief = raw_doc.RawBrief(self.txt)
        self.assertEqual(brief.getFormatted(self.formatter), '@brief aword\n')


class ExtendsTest(unittest.TestCase):
    def setUp(self):
        self.t = lexer.Token('WORD', 'aword', 0, 0, 0)
        self.txt = raw_doc.RawText([self.t])
        self.formatter = raw_doc.DoxFormatter()

    def testCreation(self):
        extends = raw_doc.RawExtends(self.txt)
        self.assertEqual(extends.text, self.txt)

    def testGetType(self):
        extends = raw_doc.RawExtends(self.txt)
        self.assertEqual(extends.getType(), 'extends')

    def testGetFormatted(self):
        extends = raw_doc.RawExtends(self.txt)
        self.assertEqual(extends.getFormatted(self.formatter), '@extends aword\n')


class ImplementsTest(unittest.TestCase):
    def setUp(self):
        self.t = lexer.Token('WORD', 'aword', 0, 0, 0)
        self.txt = raw_doc.RawText([self.t])
        self.formatter = raw_doc.DoxFormatter()

    def testCreation(self):
        implements = raw_doc.RawImplements(self.txt)
        self.assertEqual(implements.text, self.txt)

    def testGetType(self):
        implements = raw_doc.RawImplements(self.txt)
        self.assertEqual(implements.getType(), 'implements')

    def testGetFormatted(self):
        implements = raw_doc.RawImplements(self.txt)
        self.assertEqual(implements.getFormatted(self.formatter), '@implements aword\n')


class SeeTest(unittest.TestCase):
    def setUp(self):
        self.t = lexer.Token('WORD', 'aword', 0, 0, 0)
        self.txt = raw_doc.RawText([self.t])
        self.formatter = raw_doc.DoxFormatter()

    def testCreation(self):
        see = raw_doc.RawSee(self.txt)
        self.assertEqual(see.text, self.txt)

    def testGetType(self):
        see = raw_doc.RawSee(self.txt)
        self.assertEqual(see.getType(), 'see')

    def testGetFormatted(self):
        see = raw_doc.RawSee(self.txt)
        self.assertEqual(see.getFormatted(self.formatter), '@see aword\n')


class ParamTest(unittest.TestCase):
    def setUp(self):
        self.tok_name = lexer.Token('WORD', 'name', 0, 0, 0)
        self.tok_text = lexer.Token('WORD', 'text', 0, 0, 0)
        self.tok_inout = lexer.Token('PARAM_IN_OUT', '[in,out]', 0, 0, 0)
        self.txt_name = raw_doc.RawText([self.tok_name])
        self.txt_text = raw_doc.RawText([self.tok_text])
        self.formatter = raw_doc.DoxFormatter()

    def testCreation(self):
        param = raw_doc.RawParam(self.txt_name, self.txt_text)
        self.assertEqual(param.name, self.txt_name)
        self.assertEqual(param.text, self.txt_text)
        self.assertEqual(param.inout, None)

    def testCreationInOut(self):
        param = raw_doc.RawParam(self.txt_name, self.txt_text, self.tok_inout)
        self.assertEqual(param.name, self.txt_name)
        self.assertEqual(param.text, self.txt_text)
        self.assertEqual(param.inout, self.tok_inout)

    def testGetType(self):
        param = raw_doc.RawParam(self.txt_name, self.txt_text)
        self.assertEqual(param.getType(), 'param')

    def testGetFormatted(self):
        param = raw_doc.RawParam(self.txt_name, self.txt_text)
        self.assertEqual(param.getFormatted(self.formatter),
                         '@param name text\n')
        param = raw_doc.RawParam(self.txt_name, self.txt_text, self.tok_inout)
        self.assertEqual(param.getFormatted(self.formatter),
                         '@param[in,out] name text\n')


class TParamTest(unittest.TestCase):
    def setUp(self):
        self.tok_name = lexer.Token('WORD', 'name', 0, 0, 0)
        self.tok_text = lexer.Token('WORD', 'text', 0, 0, 0)
        self.txt_name = raw_doc.RawText([self.tok_name])
        self.txt_text = raw_doc.RawText([self.tok_text])
        self.formatter = raw_doc.DoxFormatter()

    def testCreation(self):
        param = raw_doc.RawTParam(self.txt_name, self.txt_text)
        self.assertEqual(param.name, self.txt_name)
        self.assertEqual(param.text, self.txt_text)

    def testGetType(self):
        param = raw_doc.RawTParam(self.txt_name, self.txt_text)
        self.assertEqual(param.getType(), 'tparam')

    def testGetFormatted(self):
        param = raw_doc.RawTParam(self.txt_name, self.txt_text)
        self.assertEqual(param.getFormatted(self.formatter),
                         '@tparam name text\n')


class ReturnTest(unittest.TestCase):
    def setUp(self):
        self.tok_name = lexer.Token('WORD', 'name', 0, 0, 0)
        self.tok_text = lexer.Token('WORD', 'text', 0, 0, 0)
        self.txt_name = raw_doc.RawText([self.tok_name])
        self.txt_text = raw_doc.RawText([self.tok_text])
        self.formatter = raw_doc.DoxFormatter()

    def testCreation(self):
        ret = raw_doc.RawReturn(self.txt_name, self.txt_text)
        self.assertEqual(ret.name, self.txt_name)
        self.assertEqual(ret.text, self.txt_text)
        self.assertEqual(ret.inout, None)

    def testGetType(self):
        ret = raw_doc.RawReturn(self.txt_name, self.txt_text)
        self.assertEqual(ret.getType(), 'return')

    def testGetFormatted(self):
        ret = raw_doc.RawReturn(self.txt_name, self.txt_text)
        self.assertEqual(ret.getFormatted(self.formatter),
                         '@return name text\n')


class SignatureTest(unittest.TestCase):
    def setUp(self):
        self.tok_text = lexer.Token('WORD', 'text', 0, 0, 0)
        self.txt_text = raw_doc.RawText([self.tok_text])
        self.formatter = raw_doc.DoxFormatter()

    def testCreation(self):
        signature = raw_doc.RawSignature(self.txt_text)
        self.assertEqual(signature.text, self.txt_text)

    def testGetType(self):
        signature = raw_doc.RawSignature(self.txt_text)
        self.assertEqual(signature.getType(), 'signature')

    def testGetFormatted(self):
        signature = raw_doc.RawSignature(self.txt_text)
        self.assertEqual(signature.getFormatted(self.formatter),
                         '@signature text\n')


if __name__ == '__main__':
    unittest.main()
