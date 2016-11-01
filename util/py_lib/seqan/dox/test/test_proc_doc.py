#!/usr/bin/env python2
"""Tests for the proc_doc module."""

import sys
import os.path
import unittest

import seqan.dox.lexer as lexer
import seqan.dox.dox_tokens as dox_tokens
import seqan.dox.dox_parser as dox_parser
import seqan.dox.proc_doc as proc_doc
import seqan.dox.raw_doc as raw_doc
import seqan.dox.pure as pure


class TextNodeTest(unittest.TestCase):
    def testRenderSimple(self):
        parent = proc_doc.TextNode(text='This is some text.')
        self.assertEqual(parent.toHtmlLike(), 'This is some text.')

    def testRenderNested(self):
        parent = proc_doc.TextNode(
            type='a', attrs={'href': 'http://www.example.com'})
        parent.addChild(proc_doc.TextNode(text='A word'))
        parent.addChild(proc_doc.TextNode(text=' does not make a '))
        parent.addChild(proc_doc.TextNode(text='sentence'))
        self.assertEqual(parent.toHtmlLike(),
                         '<a href="http://www.example.com">A word does '
                         'not make a sentence</a>')


class TestTextNodeConversion(unittest.TestCase):
    def setUp(self):
        base_dir = os.path.dirname(os.path.realpath(__file__))
        self.lexer = lexer.Lexer(dox_tokens.LEXER_TOKENS, skip_whitespace=False)
        self.doc_proc = proc_doc.DocProcessor(include_dirs=[base_dir],
                                              expected_tags=pure.EXPECTED_TAGS)
        self.conv = proc_doc.RawTextToTextNodeConverter(doc_proc=self.doc_proc,
                                                        expected_tags=pure.EXPECTED_TAGS)

    def strToTokens(self, s):
        self.lexer.input(s)
        tokens = [t for t in self.lexer.tokens()]
        return tokens[:-1]

    def testConversionPlain(self):
        r = raw_doc.RawText(self.strToTokens('This is some example.'))
        n = self.conv.run(r)
        self.assertEqual(n.toHtmlLike(), '<div>This is some example.</div>')

    def testConversionOneLevel(self):
        r = raw_doc.RawText(self.strToTokens('This <b>is</b> some example.'))
        n = self.conv.run(r)
        self.assertEqual(n.toHtmlLike(), '<div>This <b>is</b> some example.</div>')

    def testConversionOneLevel(self):
        txt = 'This is a list: <ul><li>foo</li><li>ba</li></ul>'
        r = raw_doc.RawText(self.strToTokens(txt))
        n = self.conv.run(r)
        self.assertEqual(n.toHtmlLike(), '<div>This is a list: <ul><li>foo</li><li>ba</li></ul></div>')

    def testConversionList(self):
        txt = '<ul><li>Lists</li><li>Lists again!</li></ul>'
        r = raw_doc.RawText(self.strToTokens(txt))
        n = self.conv.run(r)
        self.assertEqual(n.toHtmlLike(), '<div>%s</div>' % txt)


class TestConverterBase(unittest.TestCase):
    def createDocProcessor(self):
        base_dir = os.path.dirname(os.path.realpath(__file__))
        return proc_doc.DocProcessor(include_dirs=[base_dir],
                                     expected_tags=pure.EXPECTED_TAGS)

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
        for e in parser.documentation.entries:
            parser.documentation.filenames.append('dummy' + e.name.text)
        return parser.documentation


class TestConvertPageWithIncludes(TestConverterBase):
    """Tests for @page with @include and @snippet commands."""

    def setUp(self):
        self.proc = self.createDocProcessor()
        self.conv = self.proc.converters['page']

    def testInclude(self):
        txt = ('@page Page Page Title\n'
               '@section Example\n'
               '@include example.cpp')
        raw_page = self.parseText(txt).entries[0]
        proc_page = self.conv.process(raw_page)
        txt = ('<div><h1>Example</h1><dox:code source="include" type=".cpp" path="example.cpp">#include <iostream>\n'
               '\n'
               'int main(int arg, char const ** argv)\n'
               '{\n'
               '    //![Print to stdout]\n'
               '    std::cout << "This is an example.\\n";\n'
               '    //![Print to stdout]\n'
               '    return 0;\n'
               '}\n'
               '</dox:code></div>')
        self.assertMultiLineEqual(proc_page.body.toHtmlLike(), txt)

    def testSnippet(self):
        txt = ('@page Page Page Title\n'
               '@section Example\n'
               '@snippet example.cpp Print to stdout')
        raw_page = self.parseText(txt).entries[0]
        proc_page = self.conv.process(raw_page)
        txt = ('<div><h1>Example</h1><dox:code source="snippet" type=".cpp" path="example.cpp">'
               '    std::cout << "This is an example.\\n";'
               '</dox:code></div>')
        self.assertMultiLineEqual(proc_page.body.toHtmlLike(), txt)


class TestConvertPageWithLink(TestConverterBase):
    """Tests for @page with @link command."""

    def setUp(self):
        base_dir = os.path.dirname(os.path.realpath(__file__))
        self.proc = self.createDocProcessor()
        self.conv = self.proc.converters['page']

    def testLinkWithTitle(self):
        txt = ('@page Page Page Title\n'
               '\n'
               'A link with @link OtherPage a title @endlink.\n')
        raw_page = self.parseText(txt).entries[0]
        proc_page = self.conv.process(raw_page)
        txt = ('<div><p>A link with <a href="seqan:OtherPage">a title</a>.'
               '</p></div>')
        self.assertMultiLineEqual(proc_page.body.toHtmlLike(), txt)

    def testLinkWithoutTitle(self):
        txt = ('@page Page Page Title\n'
               '\n'
               'And a link without a title: @link OtherPage @endlink.\n')
        raw_page = self.parseText(txt).entries[0]
        proc_page = self.conv.process(raw_page)
        txt = ('<div><p>And a link without a title: <a href="seqan:OtherPage">OtherPage</a>'
               '.</p></div>')
        self.assertMultiLineEqual(proc_page.body.toHtmlLike(), txt)


class TestConvertPageWithImage(TestConverterBase):
    """Tests for @page with <img> tag."""

    def setUp(self):
        self.proc = self.createDocProcessor()
        self.conv = self.proc.converters['page']

    def testLinkWithTitle(self):
        txt = ('@page Page Page Title\n'
               '\n'
               'Here is an image: <img src="img.png" title="My image" />.\n')
        raw_page = self.parseText(txt).entries[0]
        proc_page = self.conv.process(raw_page)
        txt = ('<div><p>Here is an image: <img src="img.png" '
               'title="My image" />.</p></div>')
        self.assertMultiLineEqual(proc_page.body.toHtmlLike(), txt)

    def testLinkWithoutTitle(self):
        txt = ('@page Page Page Title\n'
               '\n'
               'And a link without a title: @link OtherPage @endlink.\n')
        raw_page = self.parseText(txt).entries[0]
        proc_page = self.conv.process(raw_page)
        txt = ('<div><p>And a link without a title: <a href="seqan:OtherPage">'
               'OtherPage</a>.</p></div>')
        self.assertMultiLineEqual(proc_page.body.toHtmlLike(), txt)


class TestConvertPage(TestConverterBase):
    def setUp(self):
        self.proc = self.createDocProcessor()
        self.conv = self.proc.converters['page']

    def testConvertMinimal(self):
        txt = '@page Page Page Title <i>italic</i>.'
        raw_page = self.parseText(txt).entries[0]
        proc_page = self.conv.process(raw_page)
        self.assertEqual(proc_page.name, 'Page')
        self.assertEqual(proc_page.kind, 'page')
        self.assertEqual(proc_page.title, 'Page Title <i>italic</i>.')

    def testConvertFull(self):
        txt = ('@page Page Page Title\n'
               '@brief This is the <i>very important</i> page brief.\n'
               '\n'
               'This is the first paragraph.\n'
               '@section First <em>heading</em>\n'
               '\n'
               'Second paragraph\n'
               '@see Link Target\n')
        raw_page = self.parseText(txt).entries[0]
        proc_page = self.conv.process(raw_page)
        self.assertEqual(proc_page.name, 'Page')
        self.assertEqual(proc_page.kind, 'page')
        self.assertEqual(proc_page.title, 'Page Title')
        txt = '<div>This is the <i>very important</i> page brief.</div>'
        self.assertEqual(proc_page.brief.toHtmlLike(), txt)
        txt = ('<div>'
               '<p>This is the first paragraph.</p>'
               '<h1>First <em>heading</em></h1>'
               '<p>Second paragraph</p>'
               '</div>'
               )
        self.assertEqual(proc_page.body.toHtmlLike(), txt)
        txt = '<a href="seqan:Link Target">Link Target</a>'
        self.assertEqual(len(proc_page.sees), 1)
        self.assertEqual(proc_page.sees[0].toHtmlLike(), txt)

    def testConvertWithCode(self):
        txt = ('@page Page Page Title\n'
               '@code{.cpp}\n'
               'int main(int argc, char const ** argv) {\n    return 0;\n}\n'
               '@endcode\n')
        raw_page = self.parseText(txt).entries[0]
        proc_page = self.conv.process(raw_page)
        txt = ('<div>'
               '<dox:code type=".cpp">'
               'int main(int argc, char const ** argv) {\n    return 0;\n}'
               '</dox:code>'
               '</div>')
        self.assertEqual(proc_page.body.toHtmlLike(), txt)


class TestConvertGroup(TestConverterBase):
    def setUp(self):
        self.proc = self.createDocProcessor()
        self.conv = self.proc.converters['defgroup']

    def testConvertMinimal(self):
        txt = '@defgroup Group Group Title <i>italic</i>.'
        raw_group = self.parseText(txt).entries[0]
        proc_group = self.conv.process(raw_group)
        self.assertEqual(proc_group.name, 'Group')
        self.assertEqual(proc_group.kind, 'group')
        self.assertEqual(proc_group.title, 'Group Title <i>italic</i>.')

    def testConvertFull(self):
        txt = ('@defgroup Group Group Title\n'
               '@brief This is the <i>very important</i> group brief.\n'
               '\n'
               'This is the first paragraph.\n'
               '@section First <em>heading</em>\n'
               '\n'
               'Second paragraph\n'
               '@see Link Target\n')
        raw_group = self.parseText(txt).entries[0]
        proc_group = self.conv.process(raw_group)
        self.assertEqual(proc_group.name, 'Group')
        self.assertEqual(proc_group.kind, 'group')
        self.assertEqual(proc_group.title, 'Group Title')
        txt = '<div>This is the <i>very important</i> group brief.</div>'
        self.assertEqual(proc_group.brief.toHtmlLike(), txt)
        txt = ('<div>'
               '<p>This is the first paragraph.</p>'
               '<h1>First <em>heading</em></h1>'
               '<p>Second paragraph</p>'
               '</div>'
               )
        self.assertEqual(proc_group.body.toHtmlLike(), txt)
        txt = '<a href="seqan:Link Target">Link Target</a>'
        self.assertEqual(len(proc_group.sees), 1)
        self.assertEqual(proc_group.sees[0].toHtmlLike(), txt)

    def testConvertWithCode(self):
        txt = ('@defgroup Group Group Title\n'
               '@code{.cpp}\n'
               'int main(int argc, char const ** argv) {\n    return 0;\n}\n'
               '@endcode\n')
        raw_group = self.parseText(txt).entries[0]
        proc_group = self.conv.process(raw_group)
        txt = ('<div>'
               '<dox:code type=".cpp">'
               'int main(int argc, char const ** argv) {\n    return 0;\n}'
               '</dox:code>'
               '</div>')
        self.assertEqual(proc_group.body.toHtmlLike(), txt)


class TestConvertEnum(TestConverterBase):
    def setUp(self):
        self.proc = self.createDocProcessor()
        self.conv = self.proc.converters['enum']

    def testConvertMinimal(self):
        txt = '@enum MyEnum My Enum'
        raw_enum = self.parseText(txt).entries[0]
        proc_enum = self.conv.process(raw_enum)
        self.assertEqual(proc_enum.name, 'MyEnum')
        self.assertEqual(proc_enum.title, 'My Enum')

    def testConvertFull(self):
        txt = ('@enum EnumName Enum Name\n'
               '@headerfile <seqan/header.h>\n'
               '@headerfile <seqan/header2.h>\n'
               '@brief This is the <i>very important</i> enum brief.\n'
               '@deprecated Deprecation msg.\n'
               '@signature Enum Name;\n'
               '\n'
               'This is the first paragraph.\n'
               '@section First <em>heading</em>\n'
               '\n'
               'Second paragraph\n'
               '@see Link Target\n')
        raw_enum = self.parseText(txt).entries[0]
        proc_enum = self.conv.process(raw_enum)
        self.assertEqual(proc_enum.name, 'EnumName')
        self.assertEqual(proc_enum.title, 'Enum Name')
        self.assertEqual(proc_enum.kind, 'enum')
        self.assertEqual(proc_enum.headerfiles, ['<seqan/header.h>', '<seqan/header2.h>'])
        self.assertEqual(len(proc_enum.signatures), 1)
        self.assertEqual(proc_enum.signatures[0].toHtmlLike(), '<div>Enum Name;</div>')
        txt = '<div>This is the <i>very important</i> enum brief.</div>'
        self.assertEqual(proc_enum.brief.toHtmlLike(), txt)
        self.assertEqual(len(proc_enum.deprecation_msgs), 1)
        self.assertEqual(proc_enum.deprecation_msgs[0].toHtmlLike(), '<div>Deprecation msg.</div>')
        txt = ('<div>'
               '<p>This is the first paragraph.</p>'
               '<h1>First <em>heading</em></h1>'
               '<p>Second paragraph</p>'
               '</div>'
               )
        self.assertEqual(proc_enum.body.toHtmlLike(), txt)
        txt = '<a href="seqan:Link Target">Link Target</a>'
        self.assertEqual(len(proc_enum.sees), 1)
        self.assertEqual(proc_enum.sees[0].toHtmlLike(), txt)


class TestConvertAdaption(TestConverterBase):
    def setUp(self):
        self.proc = self.createDocProcessor()
        self.conv = self.proc.converters['adaption']

    def testConvertMinimal(self):
        txt = '@adaption MyAdaption My Adaption'
        raw_adaption = self.parseText(txt).entries[0]
        proc_adaption = self.conv.process(raw_adaption)
        self.assertEqual(proc_adaption.name, 'MyAdaption')
        self.assertEqual(proc_adaption.title, 'My Adaption')

    def testConvertFull(self):
        txt = ('@adaption AdaptionName Adaption Name\n'
               '@headerfile <seqan/header.h>\n'
               '@headerfile <seqan/header2.h>\n'
               '@brief This is the <i>very important</i> adaption brief.\n'
               '@deprecated Deprecation msg.\n'
               '@signature Adaption Name;\n'
               '\n'
               'This is the first paragraph.\n'
               '@section First <em>heading</em>\n'
               '\n'
               'Second paragraph\n'
               '@see Link Target\n')
        raw_adaption = self.parseText(txt).entries[0]
        proc_adaption = self.conv.process(raw_adaption)
        self.assertEqual(proc_adaption.name, 'AdaptionName')
        self.assertEqual(proc_adaption.title, 'Adaption Name')
        self.assertEqual(proc_adaption.kind, 'adaption')
        self.assertEqual(proc_adaption.headerfiles, ['<seqan/header.h>', '<seqan/header2.h>'])
        self.assertEqual(len(proc_adaption.signatures), 1)
        self.assertEqual(proc_adaption.signatures[0].toHtmlLike(), '<div>Adaption Name;</div>')
        txt = '<div>This is the <i>very important</i> adaption brief.</div>'
        self.assertEqual(proc_adaption.brief.toHtmlLike(), txt)
        self.assertEqual(len(proc_adaption.deprecation_msgs), 1)
        self.assertEqual(proc_adaption.deprecation_msgs[0].toHtmlLike(), '<div>Deprecation msg.</div>')
        txt = ('<div>'
               '<p>This is the first paragraph.</p>'
               '<h1>First <em>heading</em></h1>'
               '<p>Second paragraph</p>'
               '</div>'
               )
        self.assertEqual(proc_adaption.body.toHtmlLike(), txt)
        txt = '<a href="seqan:Link Target">Link Target</a>'
        self.assertEqual(len(proc_adaption.sees), 1)
        self.assertEqual(proc_adaption.sees[0].toHtmlLike(), txt)


class TestConvertTypedef(TestConverterBase):
    def setUp(self):
        self.proc = self.createDocProcessor()
        self.conv = self.proc.converters['global_typedef']

    def testConvertMinimal(self):
        txt = '@typedef MyTypedef My Typedef'
        raw_typedef = self.parseText(txt).entries[0]
        proc_typedef = self.conv.process(raw_typedef)
        self.assertEqual(proc_typedef.name, 'MyTypedef')
        self.assertEqual(proc_typedef.title, 'My Typedef')

    def testConvertFull(self):
        txt = ('@typedef TypedefName Typedef Name\n'
               '@headerfile <seqan/header.h>\n'
               '@headerfile <seqan/header2.h>\n'
               '@brief This is the <i>very important</i> typedef brief.\n'
               '@deprecated Deprecation msg.\n'
               '@signature typedef int Name;\n'
               '\n'
               'This is the first paragraph.\n'
               '@section First <em>heading</em>\n'
               '\n'
               'Second paragraph\n'
               '@see Link Target\n')
        raw_typedef = self.parseText(txt).entries[0]
        proc_typedef = self.conv.process(raw_typedef)
        self.assertEqual(proc_typedef.name, 'TypedefName')
        self.assertEqual(proc_typedef.title, 'Typedef Name')
        self.assertEqual(proc_typedef.kind, 'global_typedef')
        self.assertEqual(proc_typedef.headerfiles, ['<seqan/header.h>', '<seqan/header2.h>'])
        self.assertEqual(len(proc_typedef.signatures), 1)
        self.assertEqual(proc_typedef.signatures[0].toHtmlLike(), '<div>typedef int Name;</div>')
        txt = '<div>This is the <i>very important</i> typedef brief.</div>'
        self.assertEqual(proc_typedef.brief.toHtmlLike(), txt)
        self.assertEqual(len(proc_typedef.deprecation_msgs), 1)
        self.assertEqual(proc_typedef.deprecation_msgs[0].toHtmlLike(), '<div>Deprecation msg.</div>')
        txt = ('<div>'
               '<p>This is the first paragraph.</p>'
               '<h1>First <em>heading</em></h1>'
               '<p>Second paragraph</p>'
               '</div>'
               )
        self.assertEqual(proc_typedef.body.toHtmlLike(), txt)
        txt = '<a href="seqan:Link Target">Link Target</a>'
        self.assertEqual(len(proc_typedef.sees), 1)
        self.assertEqual(proc_typedef.sees[0].toHtmlLike(), txt)


class TestConvertConcept(TestConverterBase):
    def setUp(self):
        self.proc = self.createDocProcessor()
        self.conv = self.proc.converters['concept']

    def testConvertMinimal(self):
        txt = '@concept MyConcept My Concept'
        raw_concept = self.parseText(txt).entries[0]
        proc_concept = self.conv.process(raw_concept)
        self.assertEqual(proc_concept.name, 'MyConcept')
        self.assertEqual(proc_concept.title, 'My Concept')
        self.assertEqual(proc_concept.kind, 'concept')

    def testConvertFull(self):
        txt = ('@concept ConceptName Concept Name\n'
               '@headerfile <seqan/header.h>\n'
               '@headerfile <seqan/header2.h>\n'
               '@extends Other Concept\n'
               '@brief This is the <i>very important</i> concept brief.\n'
               '@deprecated Deprecation msg.\n'
               '@signature concept Name;\n'
               '\n'
               'This is the first paragraph.\n'
               '@section First <em>heading</em>\n'
               '\n'
               'Second paragraph\n'
               '@see Link Target\n')
        raw_concept = self.parseText(txt).entries[0]
        proc_concept = self.conv.process(raw_concept)
        self.assertEqual(proc_concept.name, 'ConceptName')
        self.assertEqual(proc_concept.title, 'Concept Name')
        self.assertEqual(proc_concept.kind, 'concept')
        self.assertEqual(proc_concept.headerfiles, ['<seqan/header.h>', '<seqan/header2.h>'])
        self.assertEqual(len(proc_concept.signatures), 1)
        self.assertEqual(proc_concept.signatures[0].toHtmlLike(), '<div>concept Name;</div>')
        self.assertEqual(len(proc_concept.extends), 1)
        self.assertEqual(proc_concept.extends[0], 'Other Concept')
        txt = '<div>This is the <i>very important</i> concept brief.</div>'
        self.assertEqual(proc_concept.brief.toHtmlLike(), txt)
        self.assertEqual(len(proc_concept.deprecation_msgs), 1)
        self.assertEqual(proc_concept.deprecation_msgs[0].toHtmlLike(), '<div>Deprecation msg.</div>')
        txt = ('<div>'
               '<p>This is the first paragraph.</p>'
               '<h1>First <em>heading</em></h1>'
               '<p>Second paragraph</p>'
               '</div>'
               )
        self.assertEqual(proc_concept.body.toHtmlLike(), txt)
        txt = '<a href="seqan:Link Target">Link Target</a>'
        self.assertEqual(len(proc_concept.sees), 1)
        self.assertEqual(proc_concept.sees[0].toHtmlLike(), txt)


class TestConvertClass(TestConverterBase):
    def setUp(self):
        self.proc = self.createDocProcessor()
        self.conv = self.proc.converters['class']

    def testConvertMinimal(self):
        txt = '@class MyClass My Class'
        raw_class = self.parseText(txt).entries[0]
        proc_class = self.conv.process(raw_class)
        self.assertEqual(proc_class.name, 'MyClass')
        self.assertEqual(proc_class.title, 'My Class')
        self.assertEqual(proc_class.kind, 'class')

    def testConvertFull(self):
        txt = ('@class ClassName Class Name\n'
               '@headerfile <seqan/header.h>\n'
               '@headerfile <seqan/header2.h>\n'
               '@implements A Class\n'
               '@extends Other Class\n'
               '@brief This is the <i>very important</i> class brief.\n'
               '@deprecated Deprecation msg.\n'
               '@signature template <typename T>\n'
               '           class Name;\n'
               '\n'
               'This is the first paragraph.\n'
               '@section First <em>heading</em>\n'
               '\n'
               'Second paragraph\n'
               '@see Link Target\n')
        raw_class = self.parseText(txt).entries[0]
        proc_class = self.conv.process(raw_class)
        self.assertEqual(proc_class.name, 'ClassName')
        self.assertEqual(proc_class.title, 'Class Name')
        self.assertEqual(proc_class.kind, 'specialization')
        self.assertEqual(proc_class.headerfiles, ['<seqan/header.h>', '<seqan/header2.h>'])
        self.assertEqual(len(proc_class.signatures), 1)
        self.assertEqual(proc_class.signatures[0].toHtmlLike(), '<div>template &lt;typename T&gt;\nclass Name;</div>')
        self.assertEqual(len(proc_class.extends), 1)
        self.assertEqual(proc_class.extends[0], 'Other Class')
        self.assertEqual(len(proc_class.implements), 1)
        self.assertEqual(proc_class.implements[0], 'A Class')
        txt = '<div>This is the <i>very important</i> class brief.</div>'
        self.assertEqual(proc_class.brief.toHtmlLike(), txt)
        self.assertEqual(len(proc_class.deprecation_msgs), 1)
        self.assertEqual(proc_class.deprecation_msgs[0].toHtmlLike(), '<div>Deprecation msg.</div>')
        txt = ('<div>'
               '<p>This is the first paragraph.</p>'
               '<h1>First <em>heading</em></h1>'
               '<p>Second paragraph</p>'
               '</div>'
               )
        self.assertEqual(proc_class.body.toHtmlLike(), txt)
        txt = '<a href="seqan:Link Target">Link Target</a>'
        self.assertEqual(len(proc_class.sees), 1)
        self.assertEqual(proc_class.sees[0].toHtmlLike(), txt)


class TestConvertTag(TestConverterBase):
    def setUp(self):
        self.proc = self.createDocProcessor()
        self.conv = self.proc.converters['tag']

    def testConvertMinimal(self):
        txt = '@tag MyTag My Tag'
        raw_tag = self.parseText(txt).entries[0]
        proc_tag = self.conv.process(raw_tag)
        self.assertEqual(proc_tag.name, 'MyTag')
        self.assertEqual(proc_tag.title, 'My Tag')
        self.assertEqual(proc_tag.kind, 'tag')

    def testConvertFull(self):
        txt = ('@tag TagName Tag Name\n'
               '@headerfile <seqan/header.h>\n'
               '@headerfile <seqan/header2.h>\n'
               '@brief This is the <i>very important</i> tag brief.\n'
               '@deprecated Deprecation msg.\n'
               '@signature typedef Tag<TagName_> TagName;\n'
               '\n'
               'This is the first paragraph.\n'
               '@section First <em>heading</em>\n'
               '\n'
               'Second paragraph\n'
               '@see Link Target\n')
        raw_tag = self.parseText(txt).entries[0]
        proc_tag = self.conv.process(raw_tag)
        self.assertEqual(proc_tag.name, 'TagName')
        self.assertEqual(proc_tag.title, 'Tag Name')
        self.assertEqual(proc_tag.kind, 'tag')
        self.assertEqual(proc_tag.headerfiles, ['<seqan/header.h>', '<seqan/header2.h>'])
        self.assertEqual(len(proc_tag.signatures), 1)
        self.assertEqual(proc_tag.signatures[0].toHtmlLike(), '<div>typedef Tag&lt;TagName_&gt; TagName;</div>')
        txt = '<div>This is the <i>very important</i> tag brief.</div>'
        self.assertEqual(proc_tag.brief.toHtmlLike(), txt)
        self.assertEqual(len(proc_tag.deprecation_msgs), 1)
        self.assertEqual(proc_tag.deprecation_msgs[0].toHtmlLike(), '<div>Deprecation msg.</div>')
        txt = ('<div>'
               '<p>This is the first paragraph.</p>'
               '<h1>First <em>heading</em></h1>'
               '<p>Second paragraph</p>'
               '</div>'
               )
        self.assertEqual(proc_tag.body.toHtmlLike(), txt)
        txt = '<a href="seqan:Link Target">Link Target</a>'
        self.assertEqual(len(proc_tag.sees), 1)
        self.assertEqual(proc_tag.sees[0].toHtmlLike(), txt)

    def testConvertFullGrouped(self):
        txt = ('@tag Group#TagName Tag Name\n'
               '@headerfile <seqan/header.h>\n'
               '@headerfile <seqan/header2.h>\n'
               '@brief This is the <i>very important</i> tag brief.\n'
               '@deprecated Deprecation msg.\n'
               '@signature typedef Tag<TagName_> TagName;\n'
               '\n'
               'This is the first paragraph.\n'
               '@section First <em>heading</em>\n'
               '\n'
               'Second paragraph\n'
               '@see Link Target\n')
        raw_tag = self.parseText(txt).entries[0]
        proc_tag = self.conv.process(raw_tag)
        self.assertEqual(proc_tag.name, 'Group#TagName')
        self.assertEqual(proc_tag.title, 'Tag Name')
        self.assertEqual(proc_tag.kind, 'grouped_tag')
        self.assertEqual(proc_tag.headerfiles, ['<seqan/header.h>', '<seqan/header2.h>'])
        self.assertEqual(len(proc_tag.signatures), 1)
        self.assertEqual(proc_tag.signatures[0].toHtmlLike(), '<div>typedef Tag&lt;TagName_&gt; TagName;</div>')
        txt = '<div>This is the <i>very important</i> tag brief.</div>'
        self.assertEqual(proc_tag.brief.toHtmlLike(), txt)
        self.assertEqual(len(proc_tag.deprecation_msgs), 1)
        self.assertEqual(proc_tag.deprecation_msgs[0].toHtmlLike(), '<div>Deprecation msg.</div>')
        txt = ('<div>'
               '<p>This is the first paragraph.</p>'
               '<h1>First <em>heading</em></h1>'
               '<p>Second paragraph</p>'
               '</div>'
               )
        self.assertEqual(proc_tag.body.toHtmlLike(), txt)
        txt = '<a href="seqan:Link Target">Link Target</a>'
        self.assertEqual(len(proc_tag.sees), 1)
        self.assertEqual(proc_tag.sees[0].toHtmlLike(), txt)


class TestConvertFunction(TestConverterBase):
    def setUp(self):
        self.proc = self.createDocProcessor()
        self.conv = self.proc.converters['global_function']

    def testConvertMinimalGlobal(self):
        txt = '@fn myFunction my Function'
        raw_function = self.parseText(txt).entries[0]
        proc_function = self.conv.process(raw_function)
        self.assertEqual(proc_function.name, 'myFunction')
        self.assertEqual(proc_function.title, 'my Function')
        self.assertEqual(proc_function.kind, 'global_function')

    def testConvertMinimalInterface(self):
        txt = '@fn Klass#myFunction my Function'
        raw_function = self.parseText(txt).entries[0]
        proc_function = self.conv.process(raw_function)
        self.assertEqual(proc_function.name, 'Klass#myFunction')
        self.assertEqual(proc_function.title, 'my Function')
        self.assertEqual(proc_function.kind, 'interface_function')

    def testConvertMinimalMember(self):
        txt = '@fn Klass::myFunction my Function'
        raw_function = self.parseText(txt).entries[0]
        proc_function = self.conv.process(raw_function)
        self.assertEqual(proc_function.name, 'Klass::myFunction')
        self.assertEqual(proc_function.title, 'my Function')
        self.assertEqual(proc_function.kind, 'member_function')

    def testConvertFullGlobal(self):
        txt = ('@fn myFunction my Function\n'
               '@headerfile <seqan/header.h>\n'
               '@headerfile <seqan/header2.h>\n'
               '@brief This is the <i>very important</i> class brief.\n'
               '@deprecated Deprecation msg.\n'
               '@signature template <typename T1>\n'
               '           TReturn foo(T1 x);\n'
               '@param[in] x       The parameter\n'
               '@tparam    T1      The type of the first template parameter.\n'
               '@return    TReturn The return value.\n'
               '@throw     Exception The exception type.\n'
               '@datarace  This function is thread safe.\n'
               '\n'
               'This is the first paragraph.\n'
               '@section First <em>heading</em>\n'
               '\n'
               'Second paragraph\n'
               '@see Link Target\n')
        raw_function = self.parseText(txt).entries[0]
        proc_function = self.conv.process(raw_function)
        self.assertEqual(proc_function.name, 'myFunction')
        self.assertEqual(proc_function.title, 'my Function')
        self.assertEqual(proc_function.kind, 'global_function')
        self.assertEqual(proc_function.headerfiles, ['<seqan/header.h>', '<seqan/header2.h>'])
        # params
        self.assertEqual(len(proc_function.params), 1)
        self.assertEqual(proc_function.params[0].name, 'x')
        txt = '<div>The parameter</div>'
        self.assertEqual(proc_function.params[0].desc.toHtmlLike(), txt)
        self.assertEqual(proc_function.params[0].in_out, proc_doc.ProcParam.IN)
        # tparams
        self.assertEqual(len(proc_function.tparams), 1)
        self.assertEqual(proc_function.tparams[0].type, 'T1')
        txt = '<div>The type of the first template parameter.</div>'
        self.assertEqual(proc_function.tparams[0].desc.toHtmlLike(), txt)
        # returns
        self.assertEqual(len(proc_function.returns), 1)
        txt = '<div>The return value.</div>'
        self.assertEqual(proc_function.returns[0].desc.toHtmlLike(), txt)
        # throws
        self.assertEqual(len(proc_function.throws), 1)
        txt = '<div>The exception type.</div>'
        self.assertEqual(proc_function.throws[0].desc.toHtmlLike(), txt)
        # dataraces
        self.assertEqual(len(proc_function.dataraces), 1)
        txt = '<div>This function is thread safe.</div>'
        self.assertEqual(proc_function.dataraces[0].desc.toHtmlLike(), txt)
        # brief
        txt = '<div>This is the <i>very important</i> class brief.</div>'
        self.assertEqual(proc_function.brief.toHtmlLike(), txt)
        self.assertEqual(len(proc_function.deprecation_msgs), 1)
        self.assertEqual(proc_function.deprecation_msgs[0].toHtmlLike(), '<div>Deprecation msg.</div>')
        txt = ('<div>'
               '<p>This is the first paragraph.</p>'
               '<h1>First <em>heading</em></h1>'
               '<p>Second paragraph</p>'
               '</div>'
               )
        self.assertEqual(proc_function.body.toHtmlLike(), txt)
        txt = '<a href="seqan:Link Target">Link Target</a>'
        self.assertEqual(len(proc_function.sees), 1)
        self.assertEqual(proc_function.sees[0].toHtmlLike(), txt)

    def testConvertFullInterface(self):
        txt = ('@fn Klass#myFunction my Function\n'
               '@headerfile <seqan/header.h>\n'
               '@headerfile <seqan/header2.h>\n'
               '@brief This is the <i>very important</i> class brief.\n'
               '@deprecated Deprecation msg.\n'
               '@signature template <typename T1>\n'
               '           TReturn foo(T1 x);\n'
               '@param[in] x       The parameter\n'
               '@tparam    T1      The type of the first template parameter.\n'
               '@return    TReturn The return value.\n'
               '@throw     Excpetion The exception type.\n'
               '@datarace  This function is thread safe.\n'
               '\n'
               'This is the first paragraph.\n'
               '@section First <em>heading</em>\n'
               '\n'
               'Second paragraph\n'
               '@see Link Target\n')
        raw_function = self.parseText(txt).entries[0]
        proc_function = self.conv.process(raw_function)
        self.assertEqual(proc_function.name, 'Klass#myFunction')
        self.assertEqual(proc_function.title, 'my Function')
        self.assertEqual(proc_function.kind, 'interface_function')
        self.assertEqual(proc_function.headerfiles, ['<seqan/header.h>', '<seqan/header2.h>'])
        # params
        self.assertEqual(len(proc_function.params), 1)
        self.assertEqual(proc_function.params[0].name, 'x')
        txt = '<div>The parameter</div>'
        self.assertEqual(proc_function.params[0].desc.toHtmlLike(), txt)
        self.assertEqual(proc_function.params[0].in_out, proc_doc.ProcParam.IN)
        # tparams
        self.assertEqual(len(proc_function.tparams), 1)
        self.assertEqual(proc_function.tparams[0].type, 'T1')
        txt = '<div>The type of the first template parameter.</div>'
        self.assertEqual(proc_function.tparams[0].desc.toHtmlLike(), txt)
        # returns
        self.assertEqual(len(proc_function.returns), 1)
        txt = '<div>The return value.</div>'
        self.assertEqual(proc_function.returns[0].desc.toHtmlLike(), txt)
        # throws
        self.assertEqual(len(proc_function.throws), 1)
        txt = '<div>The exception type.</div>'
        self.assertEqual(proc_function.throws[0].desc.toHtmlLike(), txt)
        # dataraces
        self.assertEqual(len(proc_function.dataraces), 1)
        txt = '<div>This function is thread safe.</div>'
        self.assertEqual(proc_function.dataraces[0].desc.toHtmlLike(), txt)
        # brief
        txt = '<div>This is the <i>very important</i> class brief.</div>'
        self.assertEqual(proc_function.brief.toHtmlLike(), txt)
        self.assertEqual(len(proc_function.deprecation_msgs), 1)
        self.assertEqual(proc_function.deprecation_msgs[0].toHtmlLike(), '<div>Deprecation msg.</div>')
        txt = ('<div>'
               '<p>This is the first paragraph.</p>'
               '<h1>First <em>heading</em></h1>'
               '<p>Second paragraph</p>'
               '</div>'
               )
        self.assertEqual(proc_function.body.toHtmlLike(), txt)
        txt = '<a href="seqan:Link Target">Link Target</a>'
        self.assertEqual(len(proc_function.sees), 1)
        self.assertEqual(proc_function.sees[0].toHtmlLike(), txt)

    def testConvertFullMember(self):
        txt = ('@fn Klass::myFunction my Function\n'
               '@headerfile <seqan/header.h>\n'
               '@headerfile <seqan/header2.h>\n'
               '@brief This is the <i>very important</i> class brief.\n'
               '@deprecated Deprecation msg.\n'
               '@signature template <typename T1>\n'
               '           TReturn foo(T1 x);\n'
               '@param[in] x       The parameter\n'
               '@tparam    T1      The type of the first template parameter.\n'
               '@return    TReturn The return value.\n'
               '@throw     Excpetion The exception type.\n'
               '@datarace This function is thread safe.\n'
               '\n'
               'This is the first paragraph.\n'
               '@section First <em>heading</em>\n'
               '\n'
               'Second paragraph\n'
               '@see Link Target\n')
        raw_function = self.parseText(txt).entries[0]
        proc_function = self.conv.process(raw_function)
        self.assertEqual(proc_function.name, 'Klass::myFunction')
        self.assertEqual(proc_function.title, 'my Function')
        self.assertEqual(proc_function.kind, 'member_function')
        self.assertEqual(proc_function.headerfiles, ['<seqan/header.h>', '<seqan/header2.h>'])
        # params
        self.assertEqual(len(proc_function.params), 1)
        self.assertEqual(proc_function.params[0].name, 'x')
        txt = '<div>The parameter</div>'
        self.assertEqual(proc_function.params[0].desc.toHtmlLike(), txt)
        self.assertEqual(proc_function.params[0].in_out, proc_doc.ProcParam.IN)
        # tparams
        self.assertEqual(len(proc_function.tparams), 1)
        self.assertEqual(proc_function.tparams[0].type, 'T1')
        txt = '<div>The type of the first template parameter.</div>'
        self.assertEqual(proc_function.tparams[0].desc.toHtmlLike(), txt)
        # returns
        self.assertEqual(len(proc_function.returns), 1)
        txt = '<div>The return value.</div>'
        self.assertEqual(proc_function.returns[0].desc.toHtmlLike(), txt)
        # throws
        self.assertEqual(len(proc_function.throws), 1)
        txt = '<div>The exception type.</div>'
        self.assertEqual(proc_function.throws[0].desc.toHtmlLike(), txt)
        # dataraces
        self.assertEqual(len(proc_function.dataraces), 1)
        txt = '<div>This function is thread safe.</div>'
        self.assertEqual(proc_function.dataraces[0].desc.toHtmlLike(), txt)
        # brief
        txt = '<div>This is the <i>very important</i> class brief.</div>'
        self.assertEqual(proc_function.brief.toHtmlLike(), txt)
        self.assertEqual(len(proc_function.deprecation_msgs), 1)
        self.assertEqual(proc_function.deprecation_msgs[0].toHtmlLike(), '<div>Deprecation msg.</div>')
        txt = ('<div>'
               '<p>This is the first paragraph.</p>'
               '<h1>First <em>heading</em></h1>'
               '<p>Second paragraph</p>'
               '</div>'
               )
        self.assertEqual(proc_function.body.toHtmlLike(), txt)
        txt = '<a href="seqan:Link Target">Link Target</a>'
        self.assertEqual(len(proc_function.sees), 1)
        self.assertEqual(proc_function.sees[0].toHtmlLike(), txt)


class TestConvertMacro(TestConverterBase):
    def setUp(self):
        self.proc = self.createDocProcessor()
        self.conv = self.proc.converters['macro']

    def testConvertMinimalGlobal(self):
        txt = '@macro MACRO macro title'
        raw_macro = self.parseText(txt).entries[0]
        proc_macro = self.conv.process(raw_macro)
        self.assertEqual(proc_macro.name, 'MACRO')
        self.assertEqual(proc_macro.title, 'macro title')
        self.assertEqual(proc_macro.kind, 'macro')

    def testConvertMinimalGrouped(self):
        txt = '@fn Group#MACRO macro title'
        raw_macro = self.parseText(txt).entries[0]
        proc_macro = self.conv.process(raw_macro)
        self.assertEqual(proc_macro.name, 'Group#MACRO')
        self.assertEqual(proc_macro.title, 'macro title')
        self.assertEqual(proc_macro.kind, 'grouped_macro')

    def testConvertFullGlobal(self):
        txt = ('@macro MACRO macro title\n'
               '@headerfile <seqan/header.h>\n'
               '@headerfile <seqan/header2.h>\n'
               '@brief This is the <i>very important</i> macro brief.\n'
               '@deprecated Deprecation msg.\n'
               '@signature MACRO(param, param2)\n'
               '@param param  The parameter.\n'
               '@param param2 The second parameter.\n'
               '@return    TReturn The return value.\n'
               '@throw     Exception The exception type.\n'
               '@datarace This function is not thread safe.\n'
               '\n'
               'This is the first paragraph.\n'
               '@section First <em>heading</em>\n'
               '\n'
               'Second paragraph\n'
               '@see Link Target\n')
        raw_macro = self.parseText(txt).entries[0]
        proc_macro = self.conv.process(raw_macro)
        self.assertEqual(proc_macro.name, 'MACRO')
        self.assertEqual(proc_macro.title, 'macro title')
        self.assertEqual(proc_macro.kind, 'macro')
        self.assertEqual(proc_macro.headerfiles, ['<seqan/header.h>', '<seqan/header2.h>'])
        # params
        self.assertEqual(len(proc_macro.params), 2)
        self.assertEqual(proc_macro.params[0].name, 'param')
        txt = '<div>The parameter.</div>'
        self.assertEqual(proc_macro.params[0].desc.toHtmlLike(), txt)
        self.assertEqual(proc_macro.params[0].in_out, None)
        self.assertEqual(proc_macro.params[1].name, 'param2')
        txt = '<div>The second parameter.</div>'
        self.assertEqual(proc_macro.params[1].desc.toHtmlLike(), txt)
        self.assertEqual(proc_macro.params[1].in_out, None)
        # returns
        self.assertEqual(len(proc_macro.returns), 1)
        txt = '<div>The return value.</div>'
        self.assertEqual(proc_macro.returns[0].desc.toHtmlLike(), txt)
        # throws
        self.assertEqual(len(proc_macro.throws), 1)
        txt = '<div>The exception type.</div>'
        self.assertEqual(proc_macro.throws[0].desc.toHtmlLike(), txt)
        # dataraces
        self.assertEqual(len(proc_macro.dataraces), 1)
        txt = '<div>This function is not thread safe.</div>'
        self.assertEqual(proc_macro.dataraces[0].desc.toHtmlLike(), txt)
        # brief
        txt = '<div>This is the <i>very important</i> macro brief.</div>'
        self.assertEqual(proc_macro.brief.toHtmlLike(), txt)
        self.assertEqual(len(proc_macro.deprecation_msgs), 1)
        self.assertEqual(proc_macro.deprecation_msgs[0].toHtmlLike(), '<div>Deprecation msg.</div>')
        txt = ('<div>'
               '<p>This is the first paragraph.</p>'
               '<h1>First <em>heading</em></h1>'
               '<p>Second paragraph</p>'
               '</div>'
               )
        self.assertEqual(proc_macro.body.toHtmlLike(), txt)
        txt = '<a href="seqan:Link Target">Link Target</a>'
        self.assertEqual(len(proc_macro.sees), 1)
        self.assertEqual(proc_macro.sees[0].toHtmlLike(), txt)

    def testConvertFullGrouped(self):
        txt = ('@macro Group#MACRO macro title\n'
               '@headerfile <seqan/header.h>\n'
               '@headerfile <seqan/header2.h>\n'
               '@brief This is the <i>very important</i> class brief.\n'
               '@deprecated Deprecation msg.\n'
               '@signature MACRO(param)\n'
               '@param param The parameter\n'
               '@return    TReturn The return value.\n'
               '@throw     Exception The exception type.\n'
               '@datarace This function is not thread safe.\n'
               '\n'
               'This is the first paragraph.\n'
               '@section First <em>heading</em>\n'
               '\n'
               'Second paragraph\n'
               '@see Link Target\n')
        raw_macro = self.parseText(txt).entries[0]
        proc_macro = self.conv.process(raw_macro)
        self.assertEqual(proc_macro.name, 'Group#MACRO')
        self.assertEqual(proc_macro.title, 'macro title')
        self.assertEqual(proc_macro.kind, 'grouped_macro')
        self.assertEqual(proc_macro.headerfiles, ['<seqan/header.h>', '<seqan/header2.h>'])
        # params
        self.assertEqual(len(proc_macro.params), 1)
        self.assertEqual(proc_macro.params[0].name, 'param')
        txt = '<div>The parameter</div>'
        self.assertEqual(proc_macro.params[0].desc.toHtmlLike(), txt)
        self.assertEqual(proc_macro.params[0].in_out, None)
        # returns
        self.assertEqual(len(proc_macro.returns), 1)
        txt = '<div>The return value.</div>'
        self.assertEqual(proc_macro.returns[0].desc.toHtmlLike(), txt)
        # throws
        self.assertEqual(len(proc_macro.throws), 1)
        txt = '<div>The exception type.</div>'
        self.assertEqual(proc_macro.throws[0].desc.toHtmlLike(), txt)
        # dataraces
        self.assertEqual(len(proc_macro.dataraces), 1)
        txt = '<div>This function is not thread safe.</div>'
        self.assertEqual(proc_macro.dataraces[0].desc.toHtmlLike(), txt)
        # brief
        txt = '<div>This is the <i>very important</i> class brief.</div>'
        self.assertEqual(proc_macro.brief.toHtmlLike(), txt)
        self.assertEqual(len(proc_macro.deprecation_msgs), 1)
        self.assertEqual(proc_macro.deprecation_msgs[0].toHtmlLike(), '<div>Deprecation msg.</div>')
        txt = ('<div>'
               '<p>This is the first paragraph.</p>'
               '<h1>First <em>heading</em></h1>'
               '<p>Second paragraph</p>'
               '</div>'
               )
        self.assertEqual(proc_macro.body.toHtmlLike(), txt)
        txt = '<a href="seqan:Link Target">Link Target</a>'
        self.assertEqual(len(proc_macro.sees), 1)
        self.assertEqual(proc_macro.sees[0].toHtmlLike(), txt)


class TestConvertMetafunction(TestConverterBase):
    def setUp(self):
        self.proc = self.createDocProcessor()
        self.conv = self.proc.converters['global_metafunction']

    def testConvertMinimalGlobal(self):
        txt = '@mfn Metafunktion metafunction title'
        raw_function = self.parseText(txt).entries[0]
        proc_function = self.conv.process(raw_function)
        self.assertEqual(proc_function.name, 'Metafunktion')
        self.assertEqual(proc_function.title, 'metafunction title')
        self.assertEqual(proc_function.kind, 'global_metafunction')

    def testConvertMinimalInterface(self):
        txt = '@mfn Klass#Metafunktion metafunction title'
        raw_function = self.parseText(txt).entries[0]
        proc_function = self.conv.process(raw_function)
        self.assertEqual(proc_function.name, 'Klass#Metafunktion')
        self.assertEqual(proc_function.title, 'metafunction title')
        self.assertEqual(proc_function.kind, 'interface_metafunction')

    def testConvertFullGlobal(self):
        txt = ('@mfn Metafunktion metafunction title\n'
               '@headerfile <seqan/header.h>\n'
               '@headerfile <seqan/header2.h>\n'
               '@brief This is the <i>very important</i> class brief.\n'
               '@deprecated Deprecation msg.\n'
               '@signature template <typename T1>\n'
               '           Metafunktion<T1>::Type;\n'
               '@tparam    T1   The type of the first template parameter.\n'
               '@return    Type The return type.\n'
               '\n'
               'This is the first paragraph.\n'
               '@section First <em>heading</em>\n'
               '\n'
               'Second paragraph\n'
               '@see Link Target\n')
        raw_metafunction = self.parseText(txt).entries[0]
        proc_metafunction = self.conv.process(raw_metafunction)
        self.assertEqual(proc_metafunction.name, 'Metafunktion')
        self.assertEqual(proc_metafunction.title, 'metafunction title')
        self.assertEqual(proc_metafunction.kind, 'global_metafunction')
        self.assertEqual(proc_metafunction.headerfiles, ['<seqan/header.h>', '<seqan/header2.h>'])
        # tparams
        self.assertEqual(len(proc_metafunction.tparams), 1)
        self.assertEqual(proc_metafunction.tparams[0].type, 'T1')
        txt = '<div>The type of the first template parameter.</div>'
        self.assertEqual(proc_metafunction.tparams[0].desc.toHtmlLike(), txt)
        # returns
        self.assertEqual(len(proc_metafunction.returns), 1)
        txt = '<div>The return type.</div>'
        self.assertEqual(proc_metafunction.returns[0].desc.toHtmlLike(), txt)
        # brief
        txt = '<div>This is the <i>very important</i> class brief.</div>'
        self.assertEqual(proc_metafunction.brief.toHtmlLike(), txt)
        self.assertEqual(len(proc_metafunction.deprecation_msgs), 1)
        self.assertEqual(proc_metafunction.deprecation_msgs[0].toHtmlLike(), '<div>Deprecation msg.</div>')
        txt = ('<div>'
               '<p>This is the first paragraph.</p>'
               '<h1>First <em>heading</em></h1>'
               '<p>Second paragraph</p>'
               '</div>'
               )
        self.assertEqual(proc_metafunction.body.toHtmlLike(), txt)
        txt = '<a href="seqan:Link Target">Link Target</a>'
        self.assertEqual(len(proc_metafunction.sees), 1)
        self.assertEqual(proc_metafunction.sees[0].toHtmlLike(), txt)

    def testConvertFullInterface(self):
        txt = ('@fn Klass#Metafunktion metafunction title\n'
               '@headerfile <seqan/header.h>\n'
               '@headerfile <seqan/header2.h>\n'
               '@brief This is the <i>very important</i> class brief.\n'
               '@deprecated Deprecation msg.\n'
               '@signature template <typename T1>\n'
               '           Metafunktion<T1>::Type;\n'
               '@tparam T1   The type of the first template parameter.\n'
               '@return Type The return type.\n'
               '\n'
               'This is the first paragraph.\n'
               '@section First <em>heading</em>\n'
               '\n'
               'Second paragraph\n'
               '@see Link Target\n')
        raw_metafunction = self.parseText(txt).entries[0]
        proc_metafunction = self.conv.process(raw_metafunction)
        self.assertEqual(proc_metafunction.name, 'Klass#Metafunktion')
        self.assertEqual(proc_metafunction.title, 'metafunction title')
        self.assertEqual(proc_metafunction.kind, 'interface_metafunction')
        # tparams
        self.assertEqual(len(proc_metafunction.tparams), 1)
        self.assertEqual(proc_metafunction.tparams[0].type, 'T1')
        txt = '<div>The type of the first template parameter.</div>'
        self.assertEqual(proc_metafunction.tparams[0].desc.toHtmlLike(), txt)
        # returns
        self.assertEqual(len(proc_metafunction.returns), 1)
        txt = '<div>The return type.</div>'
        self.assertEqual(proc_metafunction.returns[0].desc.toHtmlLike(), txt)
        # brief
        txt = '<div>This is the <i>very important</i> class brief.</div>'
        self.assertEqual(proc_metafunction.brief.toHtmlLike(), txt)
        self.assertEqual(len(proc_metafunction.deprecation_msgs), 1)
        self.assertEqual(proc_metafunction.deprecation_msgs[0].toHtmlLike(), '<div>Deprecation msg.</div>')
        txt = ('<div>'
               '<p>This is the first paragraph.</p>'
               '<h1>First <em>heading</em></h1>'
               '<p>Second paragraph</p>'
               '</div>'
               )
        self.assertEqual(proc_metafunction.body.toHtmlLike(), txt)
        txt = '<a href="seqan:Link Target">Link Target</a>'
        self.assertEqual(len(proc_metafunction.sees), 1)
        self.assertEqual(proc_metafunction.sees[0].toHtmlLike(), txt)


class TestConvertVariable(TestConverterBase):
    def setUp(self):
        self.proc = self.createDocProcessor()
        self.conv = self.proc.converters['variable']

    def testConvertMinimalGlobal(self):
        txt = '@var Type myVar'
        raw_function = self.parseText(txt).entries[0]
        proc_function = self.conv.process(raw_function)
        self.assertEqual(proc_function.name, 'myVar')
        self.assertEqual(proc_function.title, 'myVar')
        self.assertEqual(proc_function.type, 'Type')
        self.assertEqual(proc_function.kind, 'variable')

    def testConvertMinimalMember(self):
        txt = '@var Type Klass::myVar'
        raw_function = self.parseText(txt).entries[0]
        doc = proc_doc.ProcDoc
        proc_function = self.conv.process(raw_function)
        self.assertEqual(proc_function.name, 'Klass::myVar')
        self.assertEqual(proc_function.title, 'Klass::myVar')
        self.assertEqual(proc_function.type, 'Type')
        self.assertEqual(proc_function.kind, 'member_variable')

    def testConvertFullGlobal(self):
        txt = ('@var Type myVar\n'
               '@headerfile <seqan/header.h>\n'
               '@headerfile <seqan/header2.h>\n'
               '@brief This is the <i>very important</i> variable brief.\n'
               '@signature Type myVar;\n'
               '\n'
               'This is the first paragraph.\n'
               '@section First <em>heading</em>\n'
               '\n'
               'Second paragraph\n'
               '@see Link Target\n')
        raw_variable = self.parseText(txt).entries[0]
        proc_variable = self.conv.process(raw_variable)
        self.assertEqual(proc_variable.name, 'myVar')
        self.assertEqual(proc_variable.title, 'myVar')
        self.assertEqual(proc_variable.type, 'Type')
        self.assertEqual(proc_variable.kind, 'variable')
        # brief
        txt = '<div>This is the <i>very important</i> variable brief.</div>'
        self.assertEqual(proc_variable.brief.toHtmlLike(), txt)
        txt = ('<div>'
               '<p>This is the first paragraph.</p>'
               '<h1>First <em>heading</em></h1>'
               '<p>Second paragraph</p>'
               '</div>'
               )
        self.assertEqual(proc_variable.body.toHtmlLike(), txt)
        txt = '<a href="seqan:Link Target">Link Target</a>'
        self.assertEqual(len(proc_variable.sees), 1)
        self.assertEqual(proc_variable.sees[0].toHtmlLike(), txt)


    def testConvertFullMember(self):
        txt = ('@var Type Klass::myVar\n'
               '@headerfile <seqan/header.h>\n'
               '@headerfile <seqan/header2.h>\n'
               '@brief This is the <i>very important</i> variable brief.\n'
               '@deprecated Deprecation msg.\n'
               '@signature Type myVar;\n'
               '\n'
               'This is the first paragraph.\n'
               '@section First <em>heading</em>\n'
               '\n'
               'Second paragraph\n'
               '@see Link Target\n')
        raw_variable = self.parseText(txt).entries[0]
        proc_variable = self.conv.process(raw_variable)
        self.assertEqual(proc_variable.name, 'Klass::myVar')
        self.assertEqual(proc_variable.title, 'Klass::myVar')
        self.assertEqual(proc_variable.type, 'Type')
        self.assertEqual(proc_variable.kind, 'member_variable')
        self.assertEqual(proc_variable.headerfiles, ['<seqan/header.h>', '<seqan/header2.h>'])
        # brief
        txt = '<div>This is the <i>very important</i> variable brief.</div>'
        self.assertEqual(proc_variable.brief.toHtmlLike(), txt)
        self.assertEqual(len(proc_variable.deprecation_msgs), 1)
        self.assertEqual(proc_variable.deprecation_msgs[0].toHtmlLike(), '<div>Deprecation msg.</div>')
        txt = ('<div>'
               '<p>This is the first paragraph.</p>'
               '<h1>First <em>heading</em></h1>'
               '<p>Second paragraph</p>'
               '</div>'
               )
        self.assertEqual(proc_variable.body.toHtmlLike(), txt)
        txt = '<a href="seqan:Link Target">Link Target</a>'
        self.assertEqual(len(proc_variable.sees), 1)
        self.assertEqual(proc_variable.sees[0].toHtmlLike(), txt)

    def testConvertFullEnumValue(self):
        txt = ('@var Enum CONSTANT\n'
               '@headerfile <seqan/header.h>\n'
               '@headerfile <seqan/header2.h>\n'
               '@brief This is the <i>very important</i> value brief.\n'
               '@deprecated Deprecation msg.\n'
               '\n'
               'This is the first paragraph.\n'
               '@section First <em>heading</em>\n'
               '\n'
               'Second paragraph\n'
               '@see Link Target\n')
        raw_variable = self.parseText(txt).entries[0]
        proc_variable = self.conv.process(raw_variable)
        self.assertEqual(proc_variable.name, 'CONSTANT')
        self.assertEqual(proc_variable.title, 'CONSTANT')
        self.assertEqual(proc_variable.type, 'Enum')
        self.assertEqual(proc_variable.kind, 'variable')
        self.assertEqual(proc_variable.headerfiles, ['<seqan/header.h>', '<seqan/header2.h>'])
        # brief
        txt = '<div>This is the <i>very important</i> value brief.</div>'
        self.assertEqual(proc_variable.brief.toHtmlLike(), txt)
        self.assertEqual(len(proc_variable.deprecation_msgs), 1)
        self.assertEqual(proc_variable.deprecation_msgs[0].toHtmlLike(), '<div>Deprecation msg.</div>')
        txt = ('<div>'
               '<p>This is the first paragraph.</p>'
               '<h1>First <em>heading</em></h1>'
               '<p>Second paragraph</p>'
               '</div>'
               )
        self.assertEqual(proc_variable.body.toHtmlLike(), txt)
        txt = '<a href="seqan:Link Target">Link Target</a>'
        self.assertEqual(len(proc_variable.sees), 1)
        self.assertEqual(proc_variable.sees[0].toHtmlLike(), txt)

    def testConvertFullMemberEnumValue(self):
        txt = ('@var Klass::Enum Klass::CONSTANT\n'
               '@headerfile <seqan/header.h>\n'
               '@headerfile <seqan/header2.h>\n'
               '@brief This is the <i>very important</i> value brief.\n'
               '@deprecated Deprecation msg.\n'
               '\n'
               'This is the first paragraph.\n'
               '@section First <em>heading</em>\n'
               '\n'
               'Second paragraph\n'
               '@see Link Target\n')
        raw_variable = self.parseText(txt).entries[0]
        proc_variable = self.conv.process(raw_variable)
        self.assertEqual(proc_variable.name, 'Klass::CONSTANT')
        self.assertEqual(proc_variable.title, 'Klass::CONSTANT')
        self.assertEqual(proc_variable.type, 'Klass::Enum')
        self.assertEqual(proc_variable.kind, 'member_variable')
        self.assertEqual(proc_variable.headerfiles, ['<seqan/header.h>', '<seqan/header2.h>'])
        # brief
        txt = '<div>This is the <i>very important</i> value brief.</div>'
        self.assertEqual(proc_variable.brief.toHtmlLike(), txt)
        self.assertEqual(len(proc_variable.deprecation_msgs), 1)
        self.assertEqual(proc_variable.deprecation_msgs[0].toHtmlLike(), '<div>Deprecation msg.</div>')
        txt = ('<div>'
               '<p>This is the first paragraph.</p>'
               '<h1>First <em>heading</em></h1>'
               '<p>Second paragraph</p>'
               '</div>'
               )
        self.assertEqual(proc_variable.body.toHtmlLike(), txt)
        txt = '<a href="seqan:Link Target">Link Target</a>'
        self.assertEqual(len(proc_variable.sees), 1)
        self.assertEqual(proc_variable.sees[0].toHtmlLike(), txt)


class TestDocProcessorInheritance(TestConverterBase):
    def setUp(self):
        self.proc = self.createDocProcessor()

    def testConceptInheritance(self):
        txt = ('@concept ConceptA1\n'
               '@brief Concept A1\n'
               '\n'
               '@concept ConceptA2\n'
               '@brief Concept A2\n'
               '\n'
               '@concept ConceptB\n'
               '@brief Concept B\n'
               '@extends ConceptA1\n'
               '@extends ConceptA2\n'
               '\n'
               '@concept ConceptC\n'
               '@brief Concept C\n'
               '@extends ConceptB\n')
        raw_doc = self.parseText(txt)
        proc_doc = self.proc.run(raw_doc)
        concept_a1 = proc_doc.top_level_entries['ConceptA1']
        self.assertEqual(concept_a1.all_extended, set())
        self.assertEqual(concept_a1.all_extending, set(['ConceptB', 'ConceptC']))
        concept_a2 = proc_doc.top_level_entries['ConceptA2']
        self.assertEqual(concept_a2.all_extended, set())
        self.assertEqual(concept_a2.all_extending, set(['ConceptB', 'ConceptC']))
        concept_b = proc_doc.top_level_entries['ConceptB']
        self.assertEqual(concept_b.all_extended, set(['ConceptA1', 'ConceptA2']))
        self.assertEqual(concept_b.all_extending, set(['ConceptC']))
        concept_c = proc_doc.top_level_entries['ConceptC']
        self.assertEqual(concept_c.all_extended, set(['ConceptA1', 'ConceptA2', 'ConceptB']))
        self.assertEqual(concept_c.all_extending, set([]))

    def testClassInheritance(self):
        txt = ('@class ClassA\n'
               '@brief Brief A\n'
               '@signature class A\n'
               '\n'
               '@class ClassB\n'
               '@brief Brief B\n'
               '@signature class B\n'
               '@extends ClassA\n'
               '\n'
               '@class ClassC\n'
               '@brief Brief C\n'
               '@signature class C\n'
               '@extends ClassB\n')
        raw_doc = self.parseText(txt)
        proc_doc = self.proc.run(raw_doc)
        class_a = proc_doc.top_level_entries['ClassA']
        self.assertEqual(class_a.all_extended, set())
        self.assertEqual(class_a.all_extending, set(['ClassB', 'ClassC']))
        class_b = proc_doc.top_level_entries['ClassB']
        self.assertEqual(class_b.all_extended, set(['ClassA']))
        self.assertEqual(class_b.all_extending, set(['ClassC']))
        class_c = proc_doc.top_level_entries['ClassC']
        self.assertEqual(class_c.all_extended, set(['ClassA', 'ClassB']))

    def testConceptClassInheritance(self):
        txt = ('@concept ConceptA\n'
               '@brief Concept A\n'
               '@signature concept A;\n'
               '\n'
               '@concept ConceptB\n'
               '@brief Concept B\n'
               '@signature concept B;\n'
               '@extends ConceptA\n'
               '\n'
               '@class ClassA\n'
               '@brief Class A\n'
               '@signature class A\n'
               '@implements ConceptB\n'
               '\n'
               '@class ClassB\n'
               '@brief Class B\n'
               '@signature class B\n'
               '@extends ClassA\n'
               '\n'
               '@class ClassC\n'
               '@brief Class C\n'
               '@signature class C\n'
               '@extends ClassB\n')
        raw_doc = self.parseText(txt)
        proc_doc = self.proc.run(raw_doc)
        concept_a = proc_doc.top_level_entries['ConceptA']
        self.assertEqual(concept_a.all_extended, set())
        self.assertEqual(concept_a.all_extending, set(['ConceptB']))
        self.assertEqual(concept_a.all_implementing, set(['ClassA', 'ClassB', 'ClassC']))
        concept_b = proc_doc.top_level_entries['ConceptB']
        self.assertEqual(concept_b.all_extended, set(['ConceptA']))
        self.assertEqual(concept_b.all_extending, set([]))
        self.assertEqual(concept_b.all_implementing, set(['ClassA', 'ClassB', 'ClassC']))
        class_a = proc_doc.top_level_entries['ClassA']
        self.assertEqual(class_a.all_extended, set())
        self.assertEqual(class_a.all_extending, set(['ClassB', 'ClassC']))
        self.assertEqual(class_a.all_implemented, set(['ConceptA', 'ConceptB']))
        class_b = proc_doc.top_level_entries['ClassB']
        self.assertEqual(class_b.all_extended, set(['ClassA']))
        self.assertEqual(class_b.all_extending, set(['ClassC']))
        self.assertEqual(class_b.all_implemented, set(['ConceptA', 'ConceptB']))
        class_c = proc_doc.top_level_entries['ClassC']
        self.assertEqual(class_c.all_extended, set(['ClassA', 'ClassB']))
        self.assertEqual(class_c.all_extending, set([]))
        self.assertEqual(class_c.all_implemented, set(['ConceptA', 'ConceptB']))


if __name__ == '__main__':
    unittest.main()
