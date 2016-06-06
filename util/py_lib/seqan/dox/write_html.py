#!/usr/bin/env python2
"""Writing for HTML pages."""

import distutils.dir_util
import json
import os
import os.path
import shutil
import sys
import xml.sax.saxutils
import urllib

import jinja2
import proc_doc
import raw_doc


def escapeForXml(s):
    """Return escaped XML of s."""
    return xml.sax.saxutils.escape(s)


def escapeName(name):
    """Escape a name such that it is safe to use for files and anchors."""
    """TODO(rmaerker): Encode special chars using urllib.quote(c.encode('utf8'))"""
    escape = '_'
    xs = []
    for c in name:
        if c.isalnum() or c in "-":
            xs.append(c)
        else:
            xs += [escape, str(ord(c))]
    return ''.join(xs)


def escapeAnchor(name):
    """Escape a name such that it is safe to use for anchors."""
    return name


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

    def getListPath(self):
        """Returns target path for list."""
        return os.path.join(self.out_dir, 'list.html')

    def translateDemoPath(self, path):
        """Translate demo path."""
        return path


class TextNodeToHtml(object):
    def __init__(self, text_node, skip_top_tag=None, start_heading=2, path_mgr=None):
        self.skip_top_tag = skip_top_tag
        self.text_node = text_node
        self.res = []
        self.start_heading = start_heading
        self.heading_table = {}
        self.path_mgr = path_mgr
        for i in range(0, 10):
            self.heading_table['h%d' % i] = 'h%d' % (i + start_heading - 1)

    def openTag(self, text_node, **kwargs):
        if text_node.raw_html:
            res = ['<', text_node.type]
        else:
            res = ['<', self.heading_table.get(text_node.type, text_node.type)]
        for key, value in text_node.attrs.iteritems():
            res += [' ', key, '=', '"', repr(value)[1:-1], '"']
        for key, value in kwargs.iteritems():
            res += [' ', key, '=', '"', value, '"']
        res.append('>')
        return res

    def closeTag(self, text_node):
        if text_node.raw_html:
            return ['</', text_node.type, '>']
        else:
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
        elif text_node.type == 'dox:code':
            target_path = text_node.attrs.get('path')
            if text_node.attrs.get('type') in ['.cpp', '.h']:
                self.res.append('<div data-src-path="%s">%s' %
                                (target_path, self.convertCode(text_node.children[0].text)))
                if self.path_mgr:
                    target_path = self.path_mrg.translateDemoPath(self.path_mgr)
                if text_node.attrs.get('source') == 'snippet':
                    self.res.append(
                        '<div class="path_label"><span class="label">Snippet from:'
                        '</span> <a href="%s" target="_top">%s</a></div>' %
                        (target_path, text_node.attrs.get('path')))
                elif text_node.attrs.get('source') == 'include':
                    self.res.append(
                        '<div class="path_label"><span class="label">Demo:'
                        '</span> <a href="%s" target="_top">%s</a></div>' %
                        (target_path, text_node.attrs.get('path')))
                self.res.append('</div>')
            elif text_node.attrs.get('type') in ['.console', '.stdout', '.stderr']:
                self.res.append('<pre class="console" data-src-path="%s">%s</pre>' %
                                (target_path, escapeForXml(text_node.children[0].text)))
            else:
                self.res.append('<pre class="code" data-src-path="%s">%s</pre>' %
                                (target_path, escapeForXml(text_node.children[0].text)))
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


def toDox(proc_entry, line_length=110, in_comment=False):
    """Process a ProcEntry into the dox-like format."""
    formatter = raw_doc.DoxFormatter()
    result = []
    result.append(proc_entry.raw_entry.getFormatted(formatter))
    for key, lst in proc_entry.subentries.iteritems():
        for elem in lst:
            result.append(elem.raw_entry.getFormatted(formatter))
    if in_comment:
        result = [' * ' + l for line in result for l in line.splitlines(False)]
        while result and result[-1] == ' * ':
            result.pop(-1)
        result = ['/*!'] + result + [' */']
    return '\n'.join(result)


def transTextNode(text_node, top=True, start_heading=3, path_mgr=None, **kwargs):
    #return text_node.toHtmlLike(skip_top_tag=not top)
    converter = TextNodeToHtml(text_node, skip_top_tag=not top, start_heading=start_heading, path_mgr=path_mgr)
    return converter.convert() or ''


def createTransLink(doc, path_mgr):
    link_converter = LinkConverter(doc)
    def transLink(entity_name, text=None):
        if not text:
            text = entity_name
        text_node = proc_doc.TextNode(type='a')
        text_node.attrs['href'] = 'seqan:%s' % entity_name
        text_node.children = [proc_doc.TextNode(text=text)]
        link_converter.visit(text_node)
        return transTextNode(text_node, path_mgr)
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
        self.env.filters['escape_anchor'] = escapeAnchor
        self.env.filters['url_encode'] = urllib.quote
        self.env.filters['transtext'] = transTextNode
        self.env.filters['to_dox'] = toDox
        self.env.filters['translink'] = createTransLink(doc, self)
        self.env.filters['name_to_path'] = createNameToPath(doc)
        self.env.filters['tojson'] = json.dumps
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
                title = entry.title
            return path, title, entry
        elif self.doc.entries.get(name):
            first, second = proc_doc.splitSecondLevelEntry(name)
            father = self.doc.top_level_entries.get(first)
            entry = self.doc.second_level_entries.get(name)
            path = '%s_%s.html#%s' % (father.kind, escapeName(father.name), escapeAnchor(name))
            return path, name, entry
        else:
            return None, None, None


# TODO(holtgrew): Should be doable in a simpler way than recursing ourselves here.  Visitor pattern for TextNode?
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
        target = a_node.attrs['href'][6:]
        target_path, target_title, target_obj = self.path_converter.convert(target)
        # Shorten path title if not manually specified.
        if (a_node.children and a_node.plainText == target_title and
           self.doc.local_name_counter.get(target_title, 1) <= 1):
            short_title = proc_doc.splitSecondLevelEntry(target_title)[1]
            a_node.children = [proc_doc.TextNode(text=short_title)]
        if target_title:
            target_title = proc_doc.TextNode(text=target_title)
        else:
            target_title = proc_doc.TextNode(text=target)
        # TODO(holtgrew): Catch target_title being None, target_path not found!
        if target_path is not None:
            if target_obj:
                a_node.attrs['data-lang-entity'] = target_obj.kind
            a_node.attrs['href'] = target_path
            if not a_node.children:
                a_node.addChild(target_title)
        else:
            class_attr = a_node.attrs.get('class', '')
            if class_attr:
                class_attr += ' '
            class_attr += 'error'
            a_node.attrs['class'] = class_attr
            if a_node.attrs.get('href'):
                del a_node.attrs['href']
            #a_node.addChild(target_title)
    
    def _replaceNode(self, text_node):
        if text_node.type == '<text>':
            return text_node
        if text_node.type == 'a':
            self._translateLink(text_node)
        for i, c in enumerate(text_node.children):
            text_node.children[i] = self._replaceNode(c)
        return text_node


class ImagePathUpdater(proc_doc.TextNodeVisitor):
    """Update image paths to target image path."""
    
    def __init__(self, doc, prefix):
        self.doc = doc
        self.prefix = prefix

    def visit(self, text_node):
        if not text_node or text_node.type == '<text>':
            return
        if text_node.type == 'img':
            self._updateImagePath(text_node)
        else:
            for i, c in enumerate(text_node.children):
                text_node.children[i] = self._replaceNode(c)

    def _updateImagePath(self, img_node):
        if not img_node.attrs.get('src'):
            return  # No path.
        img_node.attrs['src'] = os.path.join(self.prefix, img_node.attrs['src'])

    def _replaceNode(self, text_node):
        if text_node.type == '<text>':
            return text_node
        if text_node.type == 'img':
            self._updateImagePath(text_node)
        for i, c in enumerate(text_node.children):
            text_node.children[i] = self._replaceNode(c)
        return text_node


class HtmlWriter(object):
    def __init__(self, doc, args, config):
        self.doc = doc
        self.out_dirs = {}
        self.args = args
        self.config = config
        # Normalize path.
        out_dir = args.out_dir
        # Generate path names.
        self.out_dirs['root'] = out_dir
        self.out_dirs['css'] = os.path.join(out_dir, 'css')
        self.out_dirs['img'] = os.path.join(out_dir, 'img')
        self.out_dirs['js'] = os.path.join(out_dir, 'js')
        self.out_dirs['lib'] = os.path.join(out_dir, 'lib')
        self.out_dirs['lists'] = os.path.join(out_dir, 'lists')
        self.out_dirs['docs'] = os.path.join(out_dir, 'docs', 'seqan')
        # Create managers.
        self.path_manager = PathManager(out_dir)
        self.tpl_manager = TemplateManager(self.path_manager, doc)
        self.path_converter = PathConverter(doc)

    def generateFor(self):
        self.log('Generating HTML documentation')
        self.log('Output Directory: %s', self.out_dirs['root'])
        self.makedirs()
        self.copyStaticFiles()
        self.copyDocImages()
        self.generateTopFrameSet()
        self.generateLists(self.doc)
        self.translateLinks(self.doc)
        self.updateImagePaths(self.doc)
        self.generatePages(self.doc, self.config)
        self.generateDemoPages(self.doc)
        self.generateSearchIndex(self.doc)
        self.generateLinkData(self.doc)
        self.generateLanguageEntities()

    def makedirs(self):
        for path in self.out_dirs.values():
            if not os.path.exists(path):
                #self.log('Creating directory %s', path)
                os.makedirs(path)

    def copyStaticFiles(self):
        """Copy static files."""
        for kind in ['css', 'js', 'img', 'lib']:
            in_dir = os.path.join(self.path_manager.this_dir, 'tpl/%s' % kind)
            out_path = self.out_dirs[kind]
            self.log('  Copying %s => %s', in_dir, out_path)
            distutils.dir_util.copy_tree(in_dir, out_path, verbose=True)

    def copyDocImages(self):
        """Copy images from paths given in --image-dir parameter."""
        for image_dir in self.args.image_dirs:
            join = os.path.join  # shortcut
            files = [f for f in os.listdir(image_dir)
                     if os.path.isfile(join(image_dir, f))]
            for f in files:
                in_path = join(image_dir, f)
                out_path = os.path.join(self.out_dirs['img'], f)
                #self.log('  Copying %s => %s', in_path, out_path)
                shutil.copy(in_path, out_path)

    def generateTopFrameSet(self):
        """Generate frameset."""
        html = self.tpl_manager.render('index.html',
                                       development=self.args.development)  # TODO(holtgrew): Add title.
        with open(self.path_manager.getTopFramesetPath(), 'w') as f:
            f.write(html)

    def generateLists(self, doc):
        """Generate top level/second level/page index."""
        with open(self.path_manager.getListPath(), 'w') as f:
            f.write(self.tpl_manager.render('list.html', doc=doc, config=self.config,
                                            development=self.args.development))

    def translateLinks(self, doc):
        link_converter = LinkConverter(doc)
        for proc_entry in doc.entries.values():
            #self.log('    * %s', proc_entry.name)
            proc_entry.visitTextNodes(link_converter)

    def updateImagePaths(self, doc):
        """Appends image output directory to src attributes."""
        updater = ImagePathUpdater(doc, 'img')
        for proc_entry in doc.entries.values():
            #self.log('    * %s', proc_entry.name)
            proc_entry.visitTextNodes(updater)

    def generatePages(self, doc, config):
        """Generate pages for proc_doc.Documentation entries."""
        try:
            import pygments, pygments.lexers, pygments.formatters
            pygments_style = pygments.formatters.HtmlFormatter().get_style_defs('.highlight')
        except ImportError:
            pygments_style = '<!-- pygments not available -->'

        for entry in doc.top_level_entries.values():
            path = self.path_manager.getEntryPath(entry)
            #self.log('Creating %s', path)
            self.generatePage(entry, path, doc, config, pygments_style)

    def generatePage(self, entry, path, doc, config, pygments_style):
        """Generate page for entry to file at path."""

        common_kwargs = {'doc': doc,
                         'config': config,
                         'development': self.args.development,
                         'pygments_style': pygments_style,
                         'entry_kind': entry.kind,
                         'entry_name': entry.name}
        if entry.kind == 'page':
            html = self.tpl_manager.render('page.html', page=entry,  **common_kwargs)
        elif entry.kind == 'concept':
            html = self.tpl_manager.render('concept.html', concept=entry,  **common_kwargs)
        elif entry.kind in ['class', 'specialization']:
            html = self.tpl_manager.render('class.html', klass=entry,  **common_kwargs)
        elif entry.kind == 'enum':
            html = self.tpl_manager.render('enum.html', enum=entry,  **common_kwargs)
        elif entry.kind == 'adaption':
            html = self.tpl_manager.render('adaption.html', adaption=entry,  **common_kwargs)
        elif entry.kind == 'shortcut':
            html = self.tpl_manager.render('shortcut.html', shortcut=entry,  **common_kwargs)
        elif entry.kind in ['global_function', 'member_function', 'interface_function']:
            html = self.tpl_manager.render('function.html', function=entry,  **common_kwargs)
        elif entry.kind in ['global_metafunction', 'interface_metafunction']:
            html = self.tpl_manager.render('metafunction.html', metafunction=entry,  **common_kwargs)
        elif entry.kind == 'group':
            html = self.tpl_manager.render('group.html', group=entry,  **common_kwargs)
        elif entry.kind == 'tag':
            html = self.tpl_manager.render('tag.html', tag=entry,  **common_kwargs)
        elif entry.kind == 'macro':
            html = self.tpl_manager.render('macro.html', macro=entry,  **common_kwargs)
        elif entry.kind == 'global_typedef':
            html = self.tpl_manager.render('typedef.html', typedef=entry,  **common_kwargs)
        elif entry.kind == 'global_variable':
            html = self.tpl_manager.render('variable.html', variable=entry,  **common_kwargs)
        elif entry.kind == 'variable':
            html = self.tpl_manager.render('variable.html', variable=entry,  **common_kwargs)
        else:
            assert False, entry.kind
        with open(self.path_manager.getEntryPath(entry), 'w') as f:
            f.write(html)

    def generateDemoPages(self, doc):
        """Copy over all demo pages."""
        file_names = set(doc.doc_processor.include_mgr.file_cache.keys() +
                         [x[0] for x in doc.doc_processor.include_mgr.snippet_cache.keys()])
        for path in sorted(file_names):
            self.generateDemoPage(path)

    def generateDemoPage(self, path):
        """Generate a demo page."""
        dirname = os.path.join(self.out_dirs['root'], os.path.dirname(path))
        if not os.path.exists(dirname):
            os.makedirs(dirname)
        in_path = self.doc.doc_processor.include_mgr.resolvePath(path)
        to_path = os.path.join(self.out_dirs['root'], path)
        #print >>sys.stderr, in_path, '=>', to_path
        shutil.copyfile(in_path, to_path)

    def generateSearchIndex(self, doc):
        """Generate the search index."""
        js = ['window.searchData = [']
        js_module = ['window.searchDataModule = [']
        for entry in doc.top_level_entries.itervalues():
            akas, subentries, headerfile = '', '', ''
            if hasattr(entry, 'akas'):
                akas = ','.join(entry.akas)
            if hasattr(entry, 'subentries'):
                subentries = []
                for t in entry.subentries.values():
                    for s in t:
                        sID = s.title
                        title = proc_doc.splitSecondLevelEntry(s.title)[1]
                        subentries.append({'type': s.kind, 'name': s.name, 'title': title, 'id': sID})
            if hasattr(entry, 'headerfiles') and len(entry.headerfiles) > 0 :
                headerfile = entry.headerfiles[0]
                headerfile = headerfile[headerfile.find("/")+1:-3]

            if entry in self.doc.doc_processor.topLevelEntry_filenames:
                delimiter = "/include/seqan/"
                srcfile = self.doc.doc_processor.topLevelEntry_filenames[entry]
                srcfile = srcfile[srcfile.find(delimiter)+len(delimiter):]
            else :
                srcfile = ""

            js.append('  {title:%s,name:%s,text:%s,akas:%s,subentries:%s,loc:%s,langEntity:%s},' %
                      (repr(entry.title), repr(entry.name), repr(""), repr(akas), repr(subentries),
                       repr(self.path_converter.convert(entry.name)[0]), repr(entry.kind)))
            js_module.append('  {definedIn:%s,srcfile:%s},' % (repr(headerfile), repr(srcfile)))
        js.append('];')
        js_module.append('];')

        with open(os.path.join(self.out_dirs['js'], 'search.data.js'), 'wb') as f:
            f.write('\n'.join(js))
        with open(os.path.join(self.out_dirs['js'], 'search.data.module.js'), 'wb') as f:
            f.write('\n'.join(js_module))

    def generateLinkData(self, doc):
        """Generate the Data for top level entry links."""
        js = ['window.lookup = {']
        for entry in doc.top_level_entries.itervalues():
            js.append('    \'%(name)s\': \'%(kind)s_%(name)s\',' % { 'name': entry.name,
                                                                     'kind': entry.kind })
        js.append('};')
        with open(os.path.join(self.out_dirs['js'], 'link.data.js'), 'wb') as f:
            f.write('\n'.join(js))

    def generateLanguageEntities(self):
        """Generate language entities JavaScript file."""
        js = self.tpl_manager.render('js/lang_entities.js', config=self.config,
                                     development=self.args.development)
        with open(os.path.join(self.out_dirs['js'], 'lang_entities.js'), 'wb') as f:
            f.write(js)
            
    def log(self, s, *args):
        print >>sys.stderr, s % args


