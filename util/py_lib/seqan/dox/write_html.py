#!/usr/bin/env python
"""Writing for HTML pages."""

import os
import os.path
import shutil
import sys
import xml.sax.saxutils

import jinja2
import proc_doc
import raw_doc


def escapeForXml(s):
    """Return escaped XML of s."""
    return xml.sax.saxutils.escape(s)


def escapeName(name):
    """Escape a name such that it is safe to use for files and anchors."""
    escape = '_'
    xs = []
    for c in name:
        if c.isalpha() or c in ['-']:
            xs.append(c)
        else:
            xs += [escape, str(ord(c))]
    return ''.join(xs)


class PathManager(object):
    """Handles the path and link generation."""
    
    def __init__(self, out_dir):
        self.out_dir = out_dir
        self.this_dir = os.path.dirname(os.path.realpath(__file__))

    def getTplPath(self, path):
        """Returns path to template."""
        return os.path.join(self.this_dir, 'tpl', path)

    def getEntryPath(self, entry):
        """Returns target path for page for entry."""
        path = '%s_%s.html' % (entry.kind, escapeName(entry.name))
        return os.path.join(self.out_dir, path)

    def getTopFramesetPath(self):
        """Returns target path for top frameset."""
        return os.path.join(self.out_dir, 'index.html')

    def getListPath(self, kind):
        """Returns target path for list."""
        return os.path.join(self.out_dir, 'lists', '%s.html' % kind)


class TextNodeToHtml(object):
    def __init__(self, text_node, skip_top_tag=None, start_heading=2):
        self.skip_top_tag = skip_top_tag
        self.text_node = text_node
        self.res = []
        self.start_heading = start_heading
        self.heading_table = {}
        for i in range(0, 10):
            self.heading_table['h%d' % i] = 'h%d' % (i + start_heading)

    def openTag(self, text_node, **kwargs):
        res = ['<', self.heading_table.get(text_node.type, text_node.type)]
        for key, value in text_node.attrs.iteritems():
            res += [' ', key, '=', '"', repr(value)[1:-1], '"']
        for key, value in kwargs.iteritems():
            res += [' ', key, '=', '"', value, '"']
        res.append('>')
        return res

    def closeTag(self, text_node):
        return ['</', self.heading_table.get(text_node.type, text_node.type), '>']

    def convertCode(self, source_code):
        # TODO(holtgrew): Interpret source type.
        try:
            import pygments, pygments.lexers, pygments.formatters
            return pygments.highlight(
                source_code, pygments.lexers.CppLexer(),
                pygments.formatters.HtmlFormatter(style='friendly'))
        except ImportError:
            return '<pre class="code">' + escapeForXml(source_code) + '</pre>'
        return 
        
    def handleTag(self, text_node):
        if text_node.type == '<text>':
            self.res.append(text_node.text)
        elif text_node.type == 'code':
            self.res.append(self.convertCode(text_node.children[0].text))
        else:
            self.res += self.openTag(text_node)
            for c in text_node.children:
                self.handleTag(c)
            self.res += self.closeTag(text_node)
        
    def convert(self):
        if not self.text_node:
            return None
        if not self.skip_top_tag:
            self.res += self.openTag(self.text_node)
        for c in self.text_node.children:
            self.handleTag(c)
        if not self.skip_top_tag:
            self.res += self.closeTag(self.text_node)
        return ''.join(self.res)


def toDox(proc_entry, line_length=110):
    """Process a ProcEntry into the dox-like format."""
    formatter = raw_doc.DoxFormatter()
    result = []
    result.append(proc_entry.raw_entry.getFormatted(formatter))
    for key, lst in proc_entry.subentries.iteritems():
        for elem in lst:
            result.append(elem.raw_entry.getFormatted(formatter))
    return '\n'.join(result)


def transTextNode(text_node, top=True, start_heading=3, **kwargs):
    #return text_node.toHtmlLike(skip_top_tag=not top)
    converter = TextNodeToHtml(text_node, skip_top_tag=not top, start_heading=start_heading)
    return converter.convert() or ''

def createTransLink(doc):
    link_converter = LinkConverter(doc)
    def transLink(entity_name):
        text_node = proc_doc.TextNode(type='a')
        text_node.attrs['href'] = 'seqan:%s' % entity_name
        text_node.children = [proc_doc.TextNode(text=entity_name)]
        link_converter.visit(text_node)
        return transTextNode(text_node)
    return transLink

def createNameToPath(doc):
    path_converter = PathConverter(doc)
    def convertPath(entry_name):
        return path_converter.convert(entry_name)[0]
    return convertPath


class TemplateManager(object):
    def __init__(self, path_manager, doc):
        self.path_manager = path_manager
        self.env = jinja2.Environment(loader=jinja2.FileSystemLoader(os.path.join(self.path_manager.this_dir, 'tpl')))
        def identity(x):
            return x
        self.env.filters['escape_name'] = escapeName
        self.env.filters['transtext'] = transTextNode
        self.env.filters['to_dox'] = toDox
        self.env.filters['translink'] = createTransLink(doc)
        self.env.filters['name_to_path'] = createNameToPath(doc)
        self.tpls = {}
        for path in ['page.html', 'concept.html']:
            self.loadTemplate(path)

    def loadTemplate(self, path):
        """Load template string at path."""
        self.tpls[path] = self.env.get_template(path)

    def render(self, path, **kwargs):
        if not path in self.tpls:
            self.loadTemplate(path)
        return self.tpls[path].render(**kwargs)


class PathConverter(object):
    """Convert entry names to URL fragments (filename + anchor)."""

    def __init__(self, doc):
        self.doc = doc

    def convert(self, name):
        """Return None, None on failure path, title otherwise."""
        if self.doc.top_level_entries.get(name):
            entry = self.doc.top_level_entries.get(name)
            path = '%s_%s.html' % (entry.kind, escapeName(entry.name))
            title = None
            if entry.kind == 'page':
                title = list(entry.title.children)
            return path, title
        elif self.doc.entries.get(name):
            print name
            first, second = proc_doc.splitSecondLevelEntry(name)
            entry = self.doc.top_level_entries.get(first)
            path = '%s_%s.html#%s' % (entry.kind, escapeName(entry.name), escapeName(name))
            return path, None
        else:
            return None, None


class LinkConverter(proc_doc.TextNodeVisitor):
    """Convert raw links to HTML-like links.

    Raw links are links of the form <a href="seqan:$target">$label</a>.
    """
    
    def __init__(self, doc):
        self.doc = doc
        self.path_converter = PathConverter(doc)

    def visit(self, text_node):
        if not text_node or text_node.type == '<text>':
            return
        if text_node.type == 'a':
            self._translateLink(text_node)
        else:
            for i, c in enumerate(text_node.children):
                text_node.children[i] = self._replaceNode(c)

    def _translateLink(self, a_node):
        if not a_node.attrs.get('href', '').startswith('seqan:'):
            return
        target_path, target_title = self.path_converter.convert(a_node.attrs['href'][6:])
        if target_path is not None:
            a_node.attrs['href'] = target_path
            if target_title:
                a_node.children = target_title
        else:
            class_attr = a_node.attrs.get('class', '')
            if class_attr:
                class_attr += ' '
            class_attr += 'error'
            a_node.attrs['class'] = class_attr
            if a_node.attrs.get('href'):
                del a_node.attrs['href']
            
    def _replaceNode(self, text_node):
        if text_node.type == '<text>':
            return text_node
        if text_node.type == 'a':
            self._translateLink(text_node)
        for i, c in enumerate(text_node.children):
            text_node.children[i] = self._replaceNode(c)
        return text_node


class HtmlWriter(object):
    def __init__(self, doc, out_dir='html'):
        self.doc = doc
        self.out_dirs = {}
        # Normalize path.
        out_dir = os.path.abspath(out_dir)
        # Generate path names.
        self.out_dirs['root'] = out_dir
        self.out_dirs['css'] = os.path.join(out_dir, 'css')
        self.out_dirs['img'] = os.path.join(out_dir, 'img')
        self.out_dirs['js'] = os.path.join(out_dir, 'js')
        self.out_dirs['lists'] = os.path.join(out_dir, 'lists')
        self.out_dirs['docs'] = os.path.join(out_dir, 'docs', 'seqan')
        # Create managers.
        self.path_manager = PathManager(out_dir)
        self.tpl_manager = TemplateManager(self.path_manager, doc)

    def generateFor(self):
        self.log('Generating HTML documentation')
        self.log('Output Directory: %s', self.out_dirs['root'])
        self.makedirs()
        self.copyFiles()
        self.generateTopFrameSet()
        self.generateLists(self.doc)
        self.translateLinks(self.doc)
        self.generatePages(self.doc)

    def makedirs(self):
        for path in self.out_dirs.values():
            if not os.path.exists(path):
                self.log('Creating directory %s', path)
                os.makedirs(path)

    def copyFiles(self):
        """Copy static files."""
        for kind in ['css', 'js', 'img']:
            in_dir = os.path.join(self.path_manager.this_dir, 'tpl/%s' % kind)
            files = [f for f in os.listdir(in_dir) if os.path.isfile(os.path.join(in_dir, f)) ]
            for f in files:
                in_path = os.path.join(in_dir, f)
                out_path = os.path.join(self.out_dirs[kind], f)
                self.log('  Copying %s => %s', in_path, out_path)
                shutil.copyfile(in_path, out_path)

    def generateTopFrameSet(self):
        """Generate frameset."""
        html = self.tpl_manager.render('index.html')  # TODO(holtgrew): Add title.
        with open(self.path_manager.getTopFramesetPath(), 'w') as f:
            f.write(html)

    def generateLists(self, doc):
        """Generate top level/second level/page index."""
        with open(self.path_manager.getListPath('first'), 'w') as f:
            f.write(self.tpl_manager.render('list_first.html', doc=doc))
        with open(self.path_manager.getListPath('second'), 'w') as f:
            f.write(self.tpl_manager.render('list_second.html', doc=doc))
        with open(self.path_manager.getListPath('page'), 'w') as f:
            f.write(self.tpl_manager.render('list_page.html', doc=doc))

    def translateLinks(self, doc):
        link_converter = LinkConverter(doc)
        for proc_entry in doc.entries.values():
            self.log('    * %s', proc_entry.name)
            proc_entry.visitTextNodes(link_converter)

    def generatePages(self, doc):
        """Generate pages for proc_doc.Documentation entries."""
        for entry in doc.top_level_entries.values():
            path = self.path_manager.getEntryPath(entry)
            self.log('Creating %s', path)
            self.generatePage(entry, path, doc)

    def generatePage(self, entry, path, doc):
        """Generate page for entry to file at path."""
        if entry.kind == 'page':
            html = self.tpl_manager.render('page.html', page=entry, doc=doc)
        elif entry.kind == 'concept':
            html = self.tpl_manager.render('concept.html', concept=entry, doc=doc)
        elif entry.kind == 'class':
            html = self.tpl_manager.render('class.html', klass=entry, doc=doc)
        elif entry.kind == 'enum':
            html = self.tpl_manager.render('enum.html', enum=entry, doc=doc)
        elif entry.kind == 'adaption':
            html = self.tpl_manager.render('adaption.html', adaption=entry, doc=doc)
        elif entry.kind == 'shortcut':
            html = self.tpl_manager.render('shortcut.html', shortcut=entry, doc=doc)
        elif entry.kind in ['global_function', 'member_function',
                            'interface_function']:
            html = self.tpl_manager.render('function.html', function=entry, doc=doc)
        elif entry.kind in ['global_metafunction', 'interface_metafunction']:
            html = self.tpl_manager.render('metafunction.html', metafunction=entry, doc=doc)
        elif entry.kind == 'group':
            html = self.tpl_manager.render('group.html', group=entry, doc=doc)
        elif entry.kind == 'tag':
            html = self.tpl_manager.render('tag.html', tag=entry, doc=doc)
        elif entry.kind == 'macro':
            html = self.tpl_manager.render('macro.html', macro=entry, doc=doc)
        elif entry.kind == 'global_typedef':
            html = self.tpl_manager.render('typedef.html', typedef=entry, doc=doc)
        elif entry.kind == 'global_variable':
            html = self.tpl_manager.render('variable.html', variable=entry, doc=doc)
        else:
            assert False, entry.kind
        with open(self.path_manager.getEntryPath(entry), 'w') as f:
            f.write(html)

    def log(self, s, *args):
        print >>sys.stderr, s % args

