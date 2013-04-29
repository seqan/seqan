#!/usr/bin/env python
"""Tests for the proc_doc module."""

import sys
import unittest

import lexer
import os.path
import dox_tokens
import dox_parser
import proc_doc
import raw_doc

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

        
class TestTextNodeConverstion(unittest.TestCase):
    def setUp(self):
        self.lexer = lexer.Lexer(dox_tokens.LEXER_TOKENS, skip_whitespace=False)
        self.conv = proc_doc.RawTextToTextNodeConverter()

    def strToTokens(self, s):
        self.lexer.input(s)
        tokens = [t for t in self.lexer.tokens() ]
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
        return parser.documentation


class TestConvertPageWithIncludes(TestConverterBase):
    """Tests for @page with @include and @snippet commands."""

    def setUp(self):
        base_dir = os.path.dirname(os.path.realpath(__file__))
        self.proc = proc_doc.DocProcessor(include_dir=base_dir)
        self.conv = self.proc.converters['page']

    def testInclude(self):
        txt = ('@page Page Page Title\n'
               '@section Example\n'
               '@include example.cpp')
        raw_page = self.parseText(txt).entries[0]
        proc_page = self.conv.process(raw_page)
        txt = ('<div><h0>Example</h0><code type=".cpp">#include <iostream>\n'
               '\n'
               'int main(int arg, char const ** argv)\n'
               '{\n'
               '    //![Print to stdout]\n'
               '    std::cout << "This is an example.\\n";\n'
               '    //![Print to stdout]\n'
               '    return 0;\n'
               '}\n'
               '</code></div>')
        self.assertMultiLineEqual(proc_page.body.toHtmlLike(), txt)
        
    def testSnippet(self):
        txt = ('@page Page Page Title\n'
               '@section Example\n'
               '@snippet example.cpp Print to stdout')
        raw_page = self.parseText(txt).entries[0]
        proc_page = self.conv.process(raw_page)
        txt = ('<div><h0>Example</h0><code type=".cpp">'
               '    std::cout << "This is an example.\\n";'
               '</code></div>')
        self.assertMultiLineEqual(proc_page.body.toHtmlLike(), txt)
        

class TestConvertPage(TestConverterBase):
    def setUp(self):
        self.proc = proc_doc.DocProcessor()
        self.conv = self.proc.converters['page']

    def testConvertMinimal(self):
        txt = '@page Page Page Title <i>italic</i>.'
        raw_page = self.parseText(txt).entries[0]
        proc_page = self.conv.process(raw_page)
        self.assertEqual(proc_page.name, 'Page')
        self.assertEqual(proc_page.kind, 'page')
        self.assertEqual(proc_page.title.toHtmlLike(),
                         '<div>Page Title <i>italic</i>.</div>')

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
        self.assertEqual(proc_page.title.toHtmlLike(), '<div>Page Title</div>')
        txt = '<div>This is the <i>very important</i> page brief.</div>'
        self.assertEqual(proc_page.brief.toHtmlLike(), txt)
        txt = ('<div>'
               '<p>This is the first paragraph.</p>'
               '<h0>First <em>heading</em></h0>'
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
               '<code type=".cpp">'
               'int main(int argc, char const ** argv) {\n    return 0;\n}'
               '</code>'
               '</div>')
        self.assertEqual(proc_page.body.toHtmlLike(), txt)


class TestConvertGroup(TestConverterBase):
    def setUp(self):
        self.proc = proc_doc.DocProcessor()
        self.conv = self.proc.converters['defgroup']

    def testConvertMinimal(self):
        txt = '@defgroup Group Group Title <i>italic</i>.'
        raw_group = self.parseText(txt).entries[0]
        proc_group = self.conv.process(raw_group)
        self.assertEqual(proc_group.name, 'Group')
        self.assertEqual(proc_group.kind, 'group')
        self.assertEqual(proc_group.title.toHtmlLike(),
                         '<div>Group Title <i>italic</i>.</div>')

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
        self.assertEqual(proc_group.title.toHtmlLike(), '<div>Group Title</div>')
        txt = '<div>This is the <i>very important</i> group brief.</div>'
        self.assertEqual(proc_group.brief.toHtmlLike(), txt)
        txt = ('<div>'
               '<p>This is the first paragraph.</p>'
               '<h0>First <em>heading</em></h0>'
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
               '<code type=".cpp">'
               'int main(int argc, char const ** argv) {\n    return 0;\n}'
               '</code>'
               '</div>')
        self.assertEqual(proc_group.body.toHtmlLike(), txt)


class TestConvertEnum(TestConverterBase):
    def setUp(self):
        self.proc = proc_doc.DocProcessor()
        self.conv = self.proc.converters['enum']

    def testConvertMinimal(self):
        txt = '@enum My Enum'
        raw_enum = self.parseText(txt).entries[0]
        proc_enum = self.conv.process(raw_enum)
        self.assertEqual(proc_enum.name, 'My Enum')

    def testConvertFull(self):
        txt = ('@enum Enum Name\n'
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
        self.assertEqual(proc_enum.name, 'Enum Name')
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
               '<h0>First <em>heading</em></h0>'
               '<p>Second paragraph</p>'
               '</div>'
               )
        self.assertEqual(proc_enum.body.toHtmlLike(), txt)
        txt = '<a href="seqan:Link Target">Link Target</a>'
        self.assertEqual(len(proc_enum.sees), 1)
        self.assertEqual(proc_enum.sees[0].toHtmlLike(), txt)


class TestConvertAdaption(TestConverterBase):
    def setUp(self):
        self.proc = proc_doc.DocProcessor()
        self.conv = self.proc.converters['adaption']

    def testConvertMinimal(self):
        txt = '@adaption My Adaption'
        raw_adaption = self.parseText(txt).entries[0]
        proc_adaption = self.conv.process(raw_adaption)
        self.assertEqual(proc_adaption.name, 'My Adaption')

    def testConvertFull(self):
        txt = ('@adaption Adaption Name\n'
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
        self.assertEqual(proc_adaption.name, 'Adaption Name')
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
               '<h0>First <em>heading</em></h0>'
               '<p>Second paragraph</p>'
               '</div>'
               )
        self.assertEqual(proc_adaption.body.toHtmlLike(), txt)
        txt = '<a href="seqan:Link Target">Link Target</a>'
        self.assertEqual(len(proc_adaption.sees), 1)
        self.assertEqual(proc_adaption.sees[0].toHtmlLike(), txt)


class TestConvertTypedef(TestConverterBase):
    def setUp(self):
        self.proc = proc_doc.DocProcessor()
        self.conv = self.proc.converters['global_typedef']

    def testConvertMinimal(self):
        txt = '@typedef My Typedef'
        raw_typedef = self.parseText(txt).entries[0]
        proc_typedef = self.conv.process(raw_typedef)
        self.assertEqual(proc_typedef.name, 'My Typedef')

    def testConvertFull(self):
        txt = ('@typedef Typedef Name\n'
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
        self.assertEqual(proc_typedef.name, 'Typedef Name')
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
               '<h0>First <em>heading</em></h0>'
               '<p>Second paragraph</p>'
               '</div>'
               )
        self.assertEqual(proc_typedef.body.toHtmlLike(), txt)
        txt = '<a href="seqan:Link Target">Link Target</a>'
        self.assertEqual(len(proc_typedef.sees), 1)
        self.assertEqual(proc_typedef.sees[0].toHtmlLike(), txt)


class TestConvertConcept(TestConverterBase):
    def setUp(self):
        self.proc = proc_doc.DocProcessor()
        self.conv = self.proc.converters['concept']

    def testConvertMinimal(self):
        txt = '@concept My Concept'
        raw_concept = self.parseText(txt).entries[0]
        proc_concept = self.conv.process(raw_concept)
        self.assertEqual(proc_concept.name, 'My Concept')
        self.assertEqual(proc_concept.kind, 'concept')

    def testConvertFull(self):
        txt = ('@concept Concept Name\n'
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
        self.assertEqual(proc_concept.name, 'Concept Name')
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
               '<h0>First <em>heading</em></h0>'
               '<p>Second paragraph</p>'
               '</div>'
               )
        self.assertEqual(proc_concept.body.toHtmlLike(), txt)
        txt = '<a href="seqan:Link Target">Link Target</a>'
        self.assertEqual(len(proc_concept.sees), 1)
        self.assertEqual(proc_concept.sees[0].toHtmlLike(), txt)


class TestConvertClass(TestConverterBase):
    def setUp(self):
        self.proc = proc_doc.DocProcessor()
        self.conv = self.proc.converters['class']

    def testConvertMinimal(self):
        txt = '@class My Class'
        raw_class = self.parseText(txt).entries[0]
        proc_class = self.conv.process(raw_class)
        self.assertEqual(proc_class.name, 'My Class')
        self.assertEqual(proc_class.kind, 'class')

    def testConvertFull(self):
        txt = ('@class Class Name\n'
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
        self.assertEqual(proc_class.name, 'Class Name')
        self.assertEqual(proc_class.kind, 'class')
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
               '<h0>First <em>heading</em></h0>'
               '<p>Second paragraph</p>'
               '</div>'
               )
        self.assertEqual(proc_class.body.toHtmlLike(), txt)
        txt = '<a href="seqan:Link Target">Link Target</a>'
        self.assertEqual(len(proc_class.sees), 1)
        self.assertEqual(proc_class.sees[0].toHtmlLike(), txt)


class TestConvertTag(TestConverterBase):
    def setUp(self):
        self.proc = proc_doc.DocProcessor()
        self.conv = self.proc.converters['tag']

    def testConvertMinimal(self):
        txt = '@tag MyTag'
        raw_tag = self.parseText(txt).entries[0]
        proc_tag = self.conv.process(raw_tag)
        self.assertEqual(proc_tag.name, 'MyTag')
        self.assertEqual(proc_tag.kind, 'tag')

    def testConvertFull(self):
        txt = ('@tag TagName\n'
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
               '<h0>First <em>heading</em></h0>'
               '<p>Second paragraph</p>'
               '</div>'
               )
        self.assertEqual(proc_tag.body.toHtmlLike(), txt)
        txt = '<a href="seqan:Link Target">Link Target</a>'
        self.assertEqual(len(proc_tag.sees), 1)
        self.assertEqual(proc_tag.sees[0].toHtmlLike(), txt)

    def testConvertFullGrouped(self):
        txt = ('@tag Group#TagName\n'
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
               '<h0>First <em>heading</em></h0>'
               '<p>Second paragraph</p>'
               '</div>'
               )
        self.assertEqual(proc_tag.body.toHtmlLike(), txt)
        txt = '<a href="seqan:Link Target">Link Target</a>'
        self.assertEqual(len(proc_tag.sees), 1)
        self.assertEqual(proc_tag.sees[0].toHtmlLike(), txt)


class TestConvertFunction(TestConverterBase):
    def setUp(self):
        self.proc = proc_doc.DocProcessor()
        self.conv = self.proc.converters['global_function']

    def testConvertMinimalGlobal(self):
        txt = '@fn myFunction'
        raw_function = self.parseText(txt).entries[0]
        proc_function = self.conv.process(raw_function)
        self.assertEqual(proc_function.name, 'myFunction')
        self.assertEqual(proc_function.kind, 'global_function')

    def testConvertMinimalInterface(self):
        txt = '@fn Klass#myFunction'
        raw_function = self.parseText(txt).entries[0]
        proc_function = self.conv.process(raw_function)
        self.assertEqual(proc_function.name, 'Klass#myFunction')
        self.assertEqual(proc_function.kind, 'interface_function')

    def testConvertMinimalMember(self):
        txt = '@fn Klass::myFunction'
        raw_function = self.parseText(txt).entries[0]
        proc_function = self.conv.process(raw_function)
        self.assertEqual(proc_function.name, 'Klass::myFunction')
        self.assertEqual(proc_function.kind, 'member_function')

    def testConvertFullGlobal(self):
        txt = ('@fn myFunction\n'
               '@headerfile <seqan/header.h>\n'
               '@headerfile <seqan/header2.h>\n'
               '@brief This is the <i>very important</i> class brief.\n'
               '@deprecated Deprecation msg.\n'
               '@signature template <typename T1>\n'
               '           TReturn foo(T1 x);\n'
               '@param[in] x       The parameter\n'
               '@tparam    T1      The type of the first template parameter.\n'
               '@return    TReturn The return value.\n'
               '\n'
               'This is the first paragraph.\n'
               '@section First <em>heading</em>\n'
               '\n'
               'Second paragraph\n'
               '@see Link Target\n')
        raw_function = self.parseText(txt).entries[0]
        proc_function = self.conv.process(raw_function)
        self.assertEqual(proc_function.name, 'myFunction')
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
        # brief
        txt = '<div>This is the <i>very important</i> class brief.</div>'
        self.assertEqual(proc_function.brief.toHtmlLike(), txt)
        self.assertEqual(len(proc_function.deprecation_msgs), 1)
        self.assertEqual(proc_function.deprecation_msgs[0].toHtmlLike(), '<div>Deprecation msg.</div>')
        txt = ('<div>'
               '<p>This is the first paragraph.</p>'
               '<h0>First <em>heading</em></h0>'
               '<p>Second paragraph</p>'
               '</div>'
               )
        self.assertEqual(proc_function.body.toHtmlLike(), txt)
        txt = '<a href="seqan:Link Target">Link Target</a>'
        self.assertEqual(len(proc_function.sees), 1)
        self.assertEqual(proc_function.sees[0].toHtmlLike(), txt)

    def testConvertFullInterface(self):
        txt = ('@fn Klass#myFunction\n'
               '@headerfile <seqan/header.h>\n'
               '@headerfile <seqan/header2.h>\n'
               '@brief This is the <i>very important</i> class brief.\n'
               '@deprecated Deprecation msg.\n'
               '@signature template <typename T1>\n'
               '           TReturn foo(T1 x);\n'
               '@param[in] x       The parameter\n'
               '@tparam    T1      The type of the first template parameter.\n'
               '@return    TReturn The return value.\n'
               '\n'
               'This is the first paragraph.\n'
               '@section First <em>heading</em>\n'
               '\n'
               'Second paragraph\n'
               '@see Link Target\n')
        raw_function = self.parseText(txt).entries[0]
        proc_function = self.conv.process(raw_function)
        self.assertEqual(proc_function.name, 'Klass#myFunction')
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
        # brief
        txt = '<div>This is the <i>very important</i> class brief.</div>'
        self.assertEqual(proc_function.brief.toHtmlLike(), txt)
        self.assertEqual(len(proc_function.deprecation_msgs), 1)
        self.assertEqual(proc_function.deprecation_msgs[0].toHtmlLike(), '<div>Deprecation msg.</div>')
        txt = ('<div>'
               '<p>This is the first paragraph.</p>'
               '<h0>First <em>heading</em></h0>'
               '<p>Second paragraph</p>'
               '</div>'
               )
        self.assertEqual(proc_function.body.toHtmlLike(), txt)
        txt = '<a href="seqan:Link Target">Link Target</a>'
        self.assertEqual(len(proc_function.sees), 1)
        self.assertEqual(proc_function.sees[0].toHtmlLike(), txt)

    def testConvertFullMember(self):
        txt = ('@fn Klass::myFunction\n'
               '@headerfile <seqan/header.h>\n'
               '@headerfile <seqan/header2.h>\n'
               '@brief This is the <i>very important</i> class brief.\n'
               '@deprecated Deprecation msg.\n'
               '@signature template <typename T1>\n'
               '           TReturn foo(T1 x);\n'
               '@param[in] x       The parameter\n'
               '@tparam    T1      The type of the first template parameter.\n'
               '@return    TReturn The return value.\n'
               '\n'
               'This is the first paragraph.\n'
               '@section First <em>heading</em>\n'
               '\n'
               'Second paragraph\n'
               '@see Link Target\n')
        raw_function = self.parseText(txt).entries[0]
        proc_function = self.conv.process(raw_function)
        self.assertEqual(proc_function.name, 'Klass::myFunction')
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
        # brief
        txt = '<div>This is the <i>very important</i> class brief.</div>'
        self.assertEqual(proc_function.brief.toHtmlLike(), txt)
        self.assertEqual(len(proc_function.deprecation_msgs), 1)
        self.assertEqual(proc_function.deprecation_msgs[0].toHtmlLike(), '<div>Deprecation msg.</div>')
        txt = ('<div>'
               '<p>This is the first paragraph.</p>'
               '<h0>First <em>heading</em></h0>'
               '<p>Second paragraph</p>'
               '</div>'
               )
        self.assertEqual(proc_function.body.toHtmlLike(), txt)
        txt = '<a href="seqan:Link Target">Link Target</a>'
        self.assertEqual(len(proc_function.sees), 1)
        self.assertEqual(proc_function.sees[0].toHtmlLike(), txt)

        
class TestConvertMacro(TestConverterBase):
    def setUp(self):
        self.proc = proc_doc.DocProcessor()
        self.conv = self.proc.converters['macro']

    def testConvertMinimalGlobal(self):
        txt = '@macro MACRO'
        raw_macro = self.parseText(txt).entries[0]
        proc_macro = self.conv.process(raw_macro)
        self.assertEqual(proc_macro.name, 'MACRO')
        self.assertEqual(proc_macro.kind, 'macro')

    def testConvertMinimalGrouped(self):
        txt = '@fn Group#MACRO'
        raw_macro = self.parseText(txt).entries[0]
        proc_macro = self.conv.process(raw_macro)
        self.assertEqual(proc_macro.name, 'Group#MACRO')
        self.assertEqual(proc_macro.kind, 'grouped_macro')

    def testConvertFullGlobal(self):
        txt = ('@macro MACRO\n'
               '@headerfile <seqan/header.h>\n'
               '@headerfile <seqan/header2.h>\n'
               '@brief This is the <i>very important</i> macro brief.\n'
               '@deprecated Deprecation msg.\n'
               '@signature MACRO(param, param2)\n'
               '@param param  The parameter.\n'
               '@param param2 The second parameter.\n'
               '@return    TReturn The return value.\n'
               '\n'
               'This is the first paragraph.\n'
               '@section First <em>heading</em>\n'
               '\n'
               'Second paragraph\n'
               '@see Link Target\n')
        raw_macro = self.parseText(txt).entries[0]
        proc_macro = self.conv.process(raw_macro)
        self.assertEqual(proc_macro.name, 'MACRO')
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
        # brief
        txt = '<div>This is the <i>very important</i> macro brief.</div>'
        self.assertEqual(proc_macro.brief.toHtmlLike(), txt)
        self.assertEqual(len(proc_macro.deprecation_msgs), 1)
        self.assertEqual(proc_macro.deprecation_msgs[0].toHtmlLike(), '<div>Deprecation msg.</div>')
        txt = ('<div>'
               '<p>This is the first paragraph.</p>'
               '<h0>First <em>heading</em></h0>'
               '<p>Second paragraph</p>'
               '</div>'
               )
        self.assertEqual(proc_macro.body.toHtmlLike(), txt)
        txt = '<a href="seqan:Link Target">Link Target</a>'
        self.assertEqual(len(proc_macro.sees), 1)
        self.assertEqual(proc_macro.sees[0].toHtmlLike(), txt)

    def testConvertFullGrouped(self):
        txt = ('@macro Group#MACRO\n'
               '@headerfile <seqan/header.h>\n'
               '@headerfile <seqan/header2.h>\n'
               '@brief This is the <i>very important</i> class brief.\n'
               '@deprecated Deprecation msg.\n'
               '@signature MACRO(param)\n'
               '@param param The parameter\n'
               '@return    TReturn The return value.\n'
               '\n'
               'This is the first paragraph.\n'
               '@section First <em>heading</em>\n'
               '\n'
               'Second paragraph\n'
               '@see Link Target\n')
        raw_macro = self.parseText(txt).entries[0]
        proc_macro = self.conv.process(raw_macro)
        self.assertEqual(proc_macro.name, 'Group#MACRO')
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
        # brief
        txt = '<div>This is the <i>very important</i> class brief.</div>'
        self.assertEqual(proc_macro.brief.toHtmlLike(), txt)
        self.assertEqual(len(proc_macro.deprecation_msgs), 1)
        self.assertEqual(proc_macro.deprecation_msgs[0].toHtmlLike(), '<div>Deprecation msg.</div>')
        txt = ('<div>'
               '<p>This is the first paragraph.</p>'
               '<h0>First <em>heading</em></h0>'
               '<p>Second paragraph</p>'
               '</div>'
               )
        self.assertEqual(proc_macro.body.toHtmlLike(), txt)
        txt = '<a href="seqan:Link Target">Link Target</a>'
        self.assertEqual(len(proc_macro.sees), 1)
        self.assertEqual(proc_macro.sees[0].toHtmlLike(), txt)

        
class TestConvertMetafunction(TestConverterBase):
    def setUp(self):
        self.proc = proc_doc.DocProcessor()
        self.conv = self.proc.converters['global_metafunction']

    def testConvertMinimalGlobal(self):
        txt = '@mfn Metafunktion'
        raw_function = self.parseText(txt).entries[0]
        proc_function = self.conv.process(raw_function)
        self.assertEqual(proc_function.name, 'Metafunktion')
        self.assertEqual(proc_function.kind, 'global_metafunction')

    def testConvertMinimalInterface(self):
        txt = '@mfn Klass#Metafunktion'
        raw_function = self.parseText(txt).entries[0]
        proc_function = self.conv.process(raw_function)
        self.assertEqual(proc_function.name, 'Klass#Metafunktion')
        self.assertEqual(proc_function.kind, 'interface_metafunction')

    def testConvertFullGlobal(self):
        txt = ('@mfn Metafunktion\n'
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
               '<h0>First <em>heading</em></h0>'
               '<p>Second paragraph</p>'
               '</div>'
               )
        self.assertEqual(proc_metafunction.body.toHtmlLike(), txt)
        txt = '<a href="seqan:Link Target">Link Target</a>'
        self.assertEqual(len(proc_metafunction.sees), 1)
        self.assertEqual(proc_metafunction.sees[0].toHtmlLike(), txt)

    def testConvertFullInterface(self):
        txt = ('@fn Klass#Metafunktion\n'
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
               '<h0>First <em>heading</em></h0>'
               '<p>Second paragraph</p>'
               '</div>'
               )
        self.assertEqual(proc_metafunction.body.toHtmlLike(), txt)
        txt = '<a href="seqan:Link Target">Link Target</a>'
        self.assertEqual(len(proc_metafunction.sees), 1)
        self.assertEqual(proc_metafunction.sees[0].toHtmlLike(), txt)


class TestConvertVariable(TestConverterBase):
    def setUp(self):
        self.proc = proc_doc.DocProcessor()
        self.conv = self.proc.converters['variable']

    def testConvertMinimalGlobal(self):
        txt = '@var Type myVar'
        raw_function = self.parseText(txt).entries[0]
        proc_function = self.conv.process(raw_function)
        self.assertEqual(proc_function.name, 'myVar')
        self.assertEqual(proc_function.type, 'Type')
        self.assertEqual(proc_function.kind, 'variable')

    def testConvertMinimalMember(self):
        txt = '@var Type Klass::myVar'
        raw_function = self.parseText(txt).entries[0]
        doc = proc_doc.ProcDoc
        proc_function = self.conv.process(raw_function)
        self.assertEqual(proc_function.name, 'Klass::myVar')
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
        self.assertEqual(proc_variable.type, 'Type')
        self.assertEqual(proc_variable.kind, 'variable')
        # brief
        txt = '<div>This is the <i>very important</i> variable brief.</div>'
        self.assertEqual(proc_variable.brief.toHtmlLike(), txt)
        txt = ('<div>'
               '<p>This is the first paragraph.</p>'
               '<h0>First <em>heading</em></h0>'
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
               '<h0>First <em>heading</em></h0>'
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
               '<h0>First <em>heading</em></h0>'
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
               '<h0>First <em>heading</em></h0>'
               '<p>Second paragraph</p>'
               '</div>'
               )
        self.assertEqual(proc_variable.body.toHtmlLike(), txt)
        txt = '<a href="seqan:Link Target">Link Target</a>'
        self.assertEqual(len(proc_variable.sees), 1)
        self.assertEqual(proc_variable.sees[0].toHtmlLike(), txt)


if __name__ == '__main__':
    unittest.main()
