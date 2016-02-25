#!/usr/bin/env python2

import sys
import os
import os.path
import shutil
import json
import string
import StringIO
import re
import datetime
import time

import pyratemp
import json

import core

# TODO(holtgrew): Take this from dddoc tree ;)
TPL_FILES = [
    # Javascript Files
    'js/jquery-1.3.2.min.js',
    'js/jquery-effect.js',
    'js/hierarchy-tree.js',
    'js/seqan-dddoc.js',
    'js/main.js',
    'js/searchdoc.js',
    'js/excanvas.compiled.js',
    # CSS Files
    'css/main.css',
    'css/panel.css',
    'css/reset.css',
    'css/hierarchy-tree.css',
    # Images
    'favicon.ico',
    'i/link.png',
    'i/arrows.png',
    'i/results_bg.png',
    'i/tree_bg.png',
    'i/seqan_logo.png',
    # Search Panel
    'panel/tree.js',
    'panel/search_index.js',
    'panel/index.html'
    ]

SITE_TITLE = 'SeqAn'
FRONT_PAGE = 'files/index.html'
FILE_DIRS = ['images']

TYPE_CLASS = 1
TYPE_FUNCTION = 2
TYPE_FILE = 3

def escapeFiles(text):
    """Escaping of file names."""
    text = text.replace("_", "__")

    ret = ""
    for i in range(len(text)):
    	if (text[i] >= 'A') and (text[i] <= 'Z'):
    		ret += "_"
    	ret += text[i]

    ret = ret.replace("\t", "_09")
    ret = ret.replace("\n", "_0a")
    ret = ret.replace("!", "_21")
    ret = ret.replace("\"", "_22")
    ret = ret.replace("#", "_23")
    ret = ret.replace("$", "_24")
    ret = ret.replace("%", "_25")
    ret = ret.replace("&", "_26")
    ret = ret.replace("'", "_27")
    ret = ret.replace("(", "_28")
    ret = ret.replace(")", "_29")
    ret = ret.replace("*", "_2a")
    ret = ret.replace("+", "_2b")
    ret = ret.replace("/", "_2f")
    ret = ret.replace(":", "_3a")
    ret = ret.replace(",", "_2c")
    ret = ret.replace("<", "_3c")
    ret = ret.replace(">", "_3e")
    ret = ret.replace("?", "_3f")
    ret = ret.replace("\\", "_5c")
    ret = ret.replace("|", "_7c")
    ret = ret.replace(" ", "+")
    
    if not ret or ret[0] == '_':
        return ret
    else:
        return '.' + ret


def translate(text):
    """Replaces some escape sequences with HTML/XML entities."""
    # Line Breaks.
    text = text.replace("\\br", "<br>")

    # Umlaute
    text = text.replace("\\\"a", "&auml;")
    text = text.replace("\\\"o", "&ouml;")
    text = text.replace("\\\"u", "&uuml;")
    text = text.replace("\\\"A", "&Auml;")
    text = text.replace("\\\"O", "&Ouml;")
    text = text.replace("\\\"U", "&Uuml;")
    text = text.replace("\\3", "&szlig;")
    text = text.replace("\\ss", "&szlig;")

    # Escaping Special Characters
    text = text.replace("\\colon", ":")
    text = text.replace("\\dot", ".")
    text = text.replace("\\at", "@")
    text = text.replace("\\pipe", "|")
    text = text.replace("\\dollar", "$")
    text = text.replace("\\quote", "\"")
    text = text.replace("\\backslash", "\\")

    # Comparisons
    text = text.replace("\\neq", "&#x2260;")
    text = text.replace("\\approx", "&#x2248;")
    text = text.replace("\\equiv", "&#x2261;")
    text = text.replace("\\leq", "&#x2264;")
    text = text.replace("\\geq", "&#x2265;")

    # Binary Operators
    text = text.replace("\\pm", "&plusmn;") #"<font face=\"symbol\">&#xb1;</font>")
    text = text.replace("\\times", "&times;") #"<font face=\"symbol\">&#xb4;</font>")
    text = text.replace("\\cdot", "&middot;") #"<font face=\"symbol\">&#xd7;</font>")
    text = text.replace("\\div", "&divide;") #"<font face=\"symbol\">&#xb8;</font>")
    text = text.replace("\\ast", "*")
    text = text.replace("\\circ", "&#x25cb;") #"<font face=\"symbol\">&#xb0;</font>")
    text = text.replace("\\otimes", "<font face=\"symbol\">&#xc4;</font>")
    text = text.replace("\\oplus", "<font face=\"symbol\">&#xc5;</font>")

    # Logic
    text = text.replace("\\exists", "<font face=\"symbol\">&#x24;</font>")
    text = text.replace("\\forall", "<font face=\"symbol\">&#x22;</font>")
    text = text.replace("\\neg", "<font face=\"symbol\">&#xd8;</font>")
    text = text.replace("\\vee", "<font face=\"symbol\">&#xda;</font>")
    text = text.replace("\\wedge", "<font face=\"symbol\">&#xd9;</font>")

    # Set Theory
    text = text.replace("\\in", "<font face=\"symbol\">&#xce;</font>")
    text = text.replace("\\ni", "<font face=\"symbol\">&#x27;</font>")
    text = text.replace("\\notin", "<font face=\"symbol\">&#xcf;</font>")
    text = text.replace("\\cup", "<font face=\"symbol\">&#xc8;</font>")
    text = text.replace("\\cap", "<font face=\"symbol\">&#xc7;</font>")
    text = text.replace("\\subset", "<font face=\"symbol\">&#xcc;</font>")
    text = text.replace("\\supset", "<font face=\"symbol\">&#xc9;</font>")
    text = text.replace("\\subseteq", "<font face=\"symbol\">&#xcd;</font>")
    text = text.replace("\\supseteq", "<font face=\"symbol\">&#xca;</font>")

    # Arrows
    text = text.replace("\\leftarrow", "&#x2190;")
    text = text.replace("\\rightarrow", "&#x2192;")
    text = text.replace("\\leftrightarrow", "&#x2194;")
    text = text.replace("\\Leftarrow", "<font face=\"symbol\">&#xdc;</font>")
    text = text.replace("\\Rightarrow", "<font face=\"symbol\">&#xde;</font>")
    text = text.replace("\\Leftrightarrow", "<font face=\"symbol\">&#xdb;</font>")

    # Special Chars
    text = text.replace("\\infty", "&#x221E;")
    text = text.replace("\\ldots", "...")
    text = text.replace("\\squared", "&#x00B2;")
    text = text.replace("\\cubic", "&#x00B3;")

    # Greek
    text = text.replace("\\Gamma", "&#x0393;")
    text = text.replace("\\Delta", "&#x0394;")
    text = text.replace("\\Theta", "&#x0398;")
    text = text.replace("\\Lambda", "&#x039b;")
    text = text.replace("\\Xi", "&#x039e;")
    text = text.replace("\\Pi", "&#x03a0;")
    text = text.replace("\\Sigma", "&#x03a3;")
    text = text.replace("\\Phi", "&#x03a6;")
    text = text.replace("\\Psi", "&#x03a8;")
    text = text.replace("\\Omega", "&#x03a9;")

    text = text.replace("\\alpha", "&#x03b1;")
    text = text.replace("\\beta", "&#x03b2;")
    text = text.replace("\\gamma", "&#x03b3;")
    text = text.replace("\\delta", "&#x03b4;")
    text = text.replace("\\epsilon", "&#x03b5;")
    text = text.replace("\\zeta", "&#x03b6;")
    text = text.replace("\\eta", "&#x03b7;")
    text = text.replace("\\theta", "&#x03b8;")
    text = text.replace("\\iota", "&#x03b9;")
    text = text.replace("\\kappa", "&#x03ba;")
    text = text.replace("\\lambda", "&#x03bb;")
    text = text.replace("\\mu", "&#x03bc;")
    text = text.replace("\\nu", "&#x03bd;")
    text = text.replace("\\xi", "&#x03be;")
    text = text.replace("\\omicron", "&#x03bf;")
    text = text.replace("\\pi", "&#x03c0;")
    text = text.replace("\\rho", "&#x03c1;")
    text = text.replace("\\varsigma", "&#x03c2;")
    text = text.replace("\\sigma", "&#x03c3;")
    text = text.replace("\\tau", "&#x03c4;")
    text = text.replace("\\upsilon", "&#x03c5;")
    text = text.replace("\\phi", "&#x03c6;")
    text = text.replace("\\chi", "&#x03c7;")
    text = text.replace("\\psi", "&#x03c8;")
    text = text.replace("\\omega", "&#x03c9;")

    return text

def getPagePath(cat, subcat, prefix=None):
    result = cat.upper() + escapeFiles(subcat) + ".html"
    if prefix:
        result = os.path.join(prefix, result)
    return result


def getIndexPagePath(tree, cat, prefix=None):
    if tree.find('globals.project.indexcategory').text().startswith(cat):
        result = 'index.html'
    else:
        result = "INDEXPAGE" + escapeFiles(cat) + ".html"
    if prefix:
        result = os.path.join(prefix, result)
    return result


class TplDocsCreator(object):
    """Contains the code for the template-driven documentation creation.

    This

    Attrs:

      tpl_path  Path to and including the tpl directory.
      outPath  Destination path for the documentation.  Existing files are
               overwritten there.
    """

    def __init__(self, error_logger, tree, tpl_path, out_path, include_dirs):
        """Initialize the DocsCreator object.

        Args:

          tree     core.DddocTree to use as information.
          tpl_path  Path to and including the tpl directory.
          out_path  Destination path for the documentation.  Existing files are
                   overwritten there.
        """
        self.error_logger = error_logger
        self.tree = tree
        self.tpl_path = tpl_path
        self.out_path = out_path
        self.include_dirs = include_dirs
        self.html_helper = HtmlHelper(self.error_logger, self.tree, out_path, self.include_dirs)

    def good(self):
        return True

    def copyTemplates(self):
        """Copy template files.

        The files in TPL_FILE are copied verbatimly.
        """
        print >>sys.stderr, 'Copying Template Files'
        for path in TPL_FILES:
            # Get target path name.
            targetPath = os.path.join(self.out_path, os.path.dirname(path))
            if not os.path.exists(targetPath):  # Make sure targetPath exists.
                os.makedirs(targetPath)
            # Copy file.
            srcPath = os.path.join(self.tpl_path, path)
            destPath = os.path.join(self.out_path, path)
            print >>sys.stderr, '  %s => %s' % (srcPath, destPath)
            shutil.copyfile(srcPath, destPath)

    def createRootIndexPage(self):
        """Create the file /index.html.

        The file is basically copied from tpl, but some small replacements are
        done using Python string templates.
        """
        print >>sys.stderr, 'Creating Index Page'
        srcPath = os.path.join(self.tpl_path, 'index.html')
        destPath = os.path.join(self.out_path, 'index.html')
        with open(srcPath, 'rb') as f:
            contents = f.read()
        template = string.Template(contents)
        mapping = {'site_title': SITE_TITLE,
                   'front_page': FRONT_PAGE}
        with open(destPath, 'wb') as f:
            f.write(template.substitute(mapping))

    def _registerCategories(self, cat_type, search_index, long_search_index, infos):
        """Helper function for createPanelSearchIndex().

        Register all categories of type cat_type in search_index,
        long_search_index, infos.
        """
        if not self.tree.find([cat_type]):
            return  # Nothing to document here!
        # TODO(holtgrew): This distinction is bogus and an artefact of the JS expected by the rdoc JS.
        if cat_type == 'Function':
            type_num = TYPE_FUNCTION
        else:
            type_num = TYPE_CLASS
        # Index functions.
        for name, node in sorted(self.tree.find([cat_type]).children.iteritems()):
            ## if node.children.has_key('hidefromindex'):
            ##     continue  # Skip those that are not to appear in search.
            summary = ''
            if node.children.has_key('summary') and node.children['summary'].texts:
                if not ':' in node.children['summary'].texts[0]:
                    continue
                summary = self.html_helper.translateMarkup(node.children['summary'].texts[0].split(':', 1)[1])
                summary = summary.replace('&amp;', '&').replace('&quot;', '"')
                summary = re.sub('<.*?>', '', summary)
                summary = summary[:100]  # TODO(holtgrew): Smarter trimming?
            # Get title of entry type.
            type_ = ' ' + ', '.join(self.tree.find(['globals', 'categories', node.path[0]]).texts)
            search_index.append(name.lower())
            long_search_index.append(name.lower())  # TODO(holtgrew): Could be parent, ...
            # Build filename.
            filename = getPagePath(cat_type, name, 'files')
            # Build list of include headers.
            includes = ''
            if node.children.has_key('include'):
                includes = ' ' + ', '.join(node.children['include'].texts)
            if type_num == TYPE_FUNCTION:
                name = name + '()'
            infos.append([name,  # Entity name.
                          includes,  # Info
                          filename,  # klass.path
                          type_,
                          summary,
                          type_num])
            ## print infos[-1]

    def createPanelSearchIndex(self):
        """Create the file /panel/search_index.js.

        This file contains a dictionary with the entries 'searchIndex',
        'longSearchIndex, 'info' that are used by searchdoc.js for providing the
        search index.
        """
        print >>sys.stderr, 'Creating Search Index'
        destPath = os.path.join(self.out_path, 'panel', 'search_index.js')
        search_index = []
        long_search_index = []
        info = []
        # Get category types to index from DddocTree tree and build the index
        # data for this.
        items = self.tree.find('globals.indexes').children.items()
        key_entries = [(x[0], x[1].entry) for x in items if x[1].entry]
        key_linenos = [(x[0], self.tree.entries[x[1][0]].line_no_begin) for x in key_entries]
        cats = [x[0] for x in sorted(key_linenos, key=lambda a:a[1])]
        for cat_type in cats:
            self._registerCategories(cat_type, search_index, long_search_index, info)
        with open(destPath, 'wb') as f:
            s = 'var search_data = {"index":{"searchIndex":%s,"longSearchIndex":%s,"info":%s}};'
            f.write(s % (json.dumps(search_index), json.dumps(long_search_index), json.dumps(info)))

    def createPanelTree(self):
        """Create the file /panel/tree.js.

        This file contains the JavaScript information for building the naviation
        tree in the left frame.
        """
        print >>sys.stderr, 'Creating Panel Tree'
        # Translate cat/element type from internal to user-readable.
        node = self.tree.find('globals.indexes')
        cats = [(k, v.text(), v) for k, v in node.children.iteritems() if not v.children.has_key('hidefromindex')]
        def getLocation(x):
            return (self.tree.entries[x[2].entry[0]].filename,
                    self.tree.entries[x[2].entry[0]].line_no_begin)
        cats.sort(key=getLocation)
        ## print cats
        trees_data = []
        index = core.buildByTypeAndCatIndex(self.tree)
        for cat, cat_title, cat_node in cats:
            ## print 'XXX', cat, cat_title
            ## if cat not in cats:
            ##     continue  # Skip if not in globals.indexes.
            ## print 'cat =', cat
            subcat_nodes = []
            for subcat in sorted(index[cat].keys()):
                if type(index[cat][subcat]) is core.DddocTreeNode:
                    node = index[cat][subcat]
                    if 'hidefromindex' in node.children:
                        continue  # Skip if not to show in tree/index.
                    ## print cat, node.key
                    filename = getPagePath(cat, node.key, 'files')
                    ## print filename
                    subcat_node = [subcat, filename, '', []]  # is uncategorized element
                    ## print 'subcat_node = ', subcat_node
                else:  # is category
                    assert type(index[cat][subcat]) is list
                    element_nodes = []
                    for node in sorted(index[cat][subcat]):
                        if 'hidefromindex' in node.children:
                            continue  # Skip if not to show in tree/index.
                        filename = getPagePath(cat, node.key, 'files')
                        ## print filename
                        element_nodes.append([node.key, filename, '', []])
                        ## print 'element_node = ', [node.key, '', '', []]
                    filename = getIndexPagePath(self.tree, cat, 'files')
                    anchor = self.html_helper.toGlossaryAnchor(subcat)
                    link = filename + '#' + anchor
                    subcat_node = [subcat, link, '', element_nodes]
                    ## print 'subcat_node = ', subcat_node
                subcat_nodes.append(subcat_node)
            filename = getIndexPagePath(self.tree, cat, 'files')
            trees_data.append([cat_title, filename, '', subcat_nodes])
            if cat_node.children.has_key('hidefromindex'):
                continue
        # Write out tree as JavaScript/JSON.
        ## print 'trees_data =', trees_data
        if not os.path.exists(os.path.join(self.out_path, 'panel')):
            os.makedirs(os.path.join(self.out_path, 'panel'))
        destPath = os.path.join(self.out_path, 'panel', 'tree.js')
        with open(destPath, 'wb') as f:
            f.write('var tree = %s;' % json.dumps(trees_data))


class MarkupParser(object):
    def __init__(self, error_logger, tree, html):
        self.error_logger = error_logger
        self.tree = tree
        self.html = html

    def parse(self, txt):
        """Transform Markup into a parse-tree."""
        ## print 'parse(%s)' % txt
        curr = 0
        result = []
        for m in re.finditer('\$|@|\|', txt):
            ## print 'curr=%d, m.start(0)=%d' % (curr, m.start(0))
            if m.start(0) < curr:
                continue  # Skip if already processed.
            if curr != m.start(0):
                result.append(('TEXT', txt[curr:m.start(0)]))
            if txt[m.start(0)] == '$':
                end_pos = txt.find('$', m.start(0) + 1)
                if end_pos == -1: end_pos = len(txt)
                result.append(('PRE', txt[m.start(0) + 1:end_pos]))
                curr = end_pos + 1
            elif txt[m.start(0)] == '@':
                end_pos = txt.find('@', m.start(0) + 1)
                if end_pos == -1: end_pos = len(txt)
                if m.start(0) + 1 < end_pos:
                    result.append(('LINK', txt[m.start(0) + 1:end_pos]))
                curr = end_pos + 1
            elif txt[m.start(0)] == '|':
                end_pos = m.start(0) + 1
                result.append(('CELL_SPLITTER', None))
                curr = end_pos
        if curr < len(txt):
            result.append(('TEXT', txt[curr:]))
        ## print txt, result
        ## print 'DONE'
        ## print '  ', result
        return result

    def toHtml(self, parseTree, node=None):
        """Translate parse tree into HTML."""
        ## print 'toHtml'
        result = []
        for type_, txt in parseTree:
            ## print '  ', type_, txt
            if type_ == 'TEXT':
                result.append(pyratemp.escape(txt))
            elif type_ == 'LINK':
                result.append(self.html.pageLink(txt, node=node))
            elif type_ == 'PRE':
                result += ['<tt>', pyratemp.escape(txt), '</tt>']
            elif type_ == 'CELL_SPLITTER':
                result.append('|')
            else:
                raise Error('Invalid type: %s' % type_)
        ## print 'DONE'
        ## print parseTree, result
        return ''.join(result)


class HtmlHelper(object):
    def __init__(self, error_logger, tree, out_path, include_dirs):
        self.error_logger = error_logger
        self.tree = tree
        self.markup_parser = MarkupParser(self.error_logger, tree, self)
        self.out_path = out_path
        self.include_dirs = include_dirs

    def classHierarchyJS(self, node):
        def recurseDownwards(node):
            cat, subcat = node.path[0], node.path[1]
            link = cat.upper() + escapeFiles(subcat) + ".html"
            result = {'title': node.key, 'type': node.path[0], 'children': [], 'parents': [], 'link': link}
            for child_key in ['spec', 'childconcept']:
                if not child_key in node.children:
                    continue
                for path in node.children[child_key].texts:
                    child = self.tree.find(path)
                    if not child:
                        self.error_logger.invalidReference(path, [('TODO', -1)])
                        return False
                    recRes = recurseDownwards(child)
                    if recRes:
                        result['children'].append(recRes)
            return result

        def recurseUpwards(node):
            cat, subcat = node.path[0], node.path[1]
            link = cat.upper() + escapeFiles(subcat) + ".html"
            result = {'title': node.key, 'type': node.path[0], 'children': [], 'parents': [], 'link': link}
            for parent_key in ['implements', 'general', 'baseconcept']:
                if not parent_key in node.children:
                    continue
                for path in node.children[parent_key].texts:
                    if '\u0001' in path:
                        continue  # Skip automatically generated upwards links.
                    parent = self.tree.find(path)
                    if not parent:
                        self.error_logger.invalidReference(path, [('TODO', -1)])
                        return False
                    recRes = recurseUpwards(parent)
                    if recRes:
                        result['parents'].append(recRes)
            return result

        res = recurseDownwards(node)
        resUp = recurseUpwards(node)
        if resUp:
            res['parents'] = resUp['parents']
        return json.dumps(res)

    def translateId(self, txt):
        i = txt.find('#');
        if (i >= 0):
            txt = txt[i + 1:]
        i = txt.find('|');
        if (i >= 0):
            txt = txt[i + 1:]
        res = self.markup_parser.toHtml(self.markup_parser.parse(txt))
        res = translate(res)
        return res

    def imageFilename(self, text):
        if '|' in text:
            return pyratemp.escape(text[len('type=image:'):].split('|', 1)[0])
        else:
            return pyratemp.escape(text[len('type=image:'):])

    def imageTitle(self, text):
        if '|' in text:
            return self.translateMarkup(text[len('type=image:'):].split('|', 1)[1])
        return ''

    def toGlossaryAnchor(self, text):
        return re.sub('[^a-zA-Z0-9\-_]', '', text)

    def _glossaryLink(self, txt):
        term = txt[len('glos:'):]
        ## print 'SEARCHING FOR TERM', term
        ## print self.tree.glossary_nodes
        for node in self.tree.glossary_nodes:
            ## print node.path
            if not term in node.children:
                continue  # Skip
            term_node = node.children[term]
            if 'description' in term_node.children:
                title = term_node.children['description'].text()
            elif 'text' in term_node.children:
                title = term_node.children['text'].text()
            else:
                title = term
            anchor = self.toGlossaryAnchor(term)
            return '<a href="#%s" class="glossary" title="%s">%s</a>' % (anchor, title, term)
        self.error_logger.invalidReference(txt, [('TODO', -1)])
        return '<a href="#" class="dead">%s</a>' % term

    def pageLink(self, txt=None, arr=None, node=None):
        """The link can be given as text or as a path.

        If it is given as text then also HTTP/FTP links are allowed, otherwise,
        it can only be a link to an entity in the tree.
        """
        # Compute source file name and line.
        location_candidates = []
        if node and node.entry:
            for entry in node.tree.entries[node.entry[0]:node.entry[1]]:
                if entry.line_no_begin + 1 == entry.line_no_end:
                    line = entry.line_no_begin + 1
                else:
                    line = '%s-%s' % (entry.line_no_begin + 1, entry.line_no_end)
                location_candidates.append((entry.filename, line))
        # Now, switch between txt and arr.
        is_dead = False
        if txt:
            # Split out titles from "$reference|$title".
            title = None
            if '|' in txt:
                txt, title = txt.split('|', 1)
            # Handle the different link types.
            if txt.startswith('glos:'):
                return self._glossaryLink(txt)
            elif txt.split(':')[0] in ['http', 'https', 'ftp']:
                if not title: title = txt
                return '<a href="%s" target="_top">%s</a>' % (pyratemp.escape(txt), pyratemp.escape(title))
            elif txt.startswith('nolink:'):
                if not title: title = txt[len('nolink:'):]
                return self.translateMarkup(title, node=node)
            else:
                # Is not a special link, compute two-element path and title.  We
                # will use the link generation code shared with paths as arrays.
                lst = core.splitKeys(txt[txt.startswith('.'):], '.')  # The startswith removes one leading dot if any.
                lst = core.cleanPath(lst)
                if len(lst) == 1:  # Is link to index.
                    cat, subcat = 'indexpage', lst[0]
                    if not title:
                        if self.tree.find(['globals', 'indexes', subcat]):
                            title = self.tree.find(['globals', 'indexes', subcat]).text()
                        else:
                            title = subcat
                    if not self.tree.find(subcat):
                        is_dead = True
                        self.error_logger.invalidReference(txt, location_candidates)
                else:
                    cat, subcat = lst[0], lst[1]
                    if not title: title = lst[-1]
                    if not self.tree.find([cat, subcat]):
                        is_dead = True
                        self.error_logger.invalidReference(txt, location_candidates)
        else:
            # Code for array paths.
            cat, subcat, title = arr[0], arr[1], arr[1]
        # Shared link generation code.
        title = self.translateId(title)
        filename = cat.upper() + escapeFiles(subcat) + ".html"
        dead_attr = {True: ' class="dead"', False: ''}[is_dead]
        return '<a href="%s"%s>%s</a>' % (pyratemp.escape(filename), dead_attr, title)

    def translateMarkup(self, text, section=None, subsection=None, node=None):
        ## print text, section, subsection
        if not section is None:
            text = text.replace('#', str(section))
        elif not section is None and not subsection is None:
            if text.count('#') == 1:
                text = text.replace('#', str(subsection))
            else:
                p1 = text.find('#')
                p2 = text.find('#', p1 + 1)
                ## print text
                text = ''.join([text[:p1], str(section), text[(p1 + 1):p2], str(subsection), text[(p2 + 1):]])
                ## print '  >', text
        res = self.markup_parser.toHtml(self.markup_parser.parse(text), node=node)
        res = translate(res)  # Translate stuff like \at.
        return res

    def translateTableHeader(self, text):
        return self.translateTableRow(text, cell_tag='th')

    def translateTableRow(self, text, cell_tag='td'):
        parse_tree = self.markup_parser.parse(text)
        cells = [[]]
        for type_, txt in parse_tree:
            if type_ == 'CELL_SPLITTER':
                cells.append([])
            else:
                cells[-1].append((type_, txt))
        res = ['<tr>']
        for c in cells:
            res += ['<', cell_tag, '>', self.markup_parser.toHtml(c), '</', cell_tag, '>']
        res.append('</tr>')
        return ''.join(res)

    def _splitIncludeFile(self, txt):
        """Splits text from included file.

        Returns a list of pairs.  Each pair contains the type of the entry in
        the first and the text of the entry in the second component.  The type
        is either 'TEXT' or 'COD'E.
        """
        result = []
        curr_type = None
        curr_chunk = []
        for line in txt.splitlines():
            if line.startswith('///'):
                if curr_type != 'TEXT':
                    if curr_chunk:
                        result.append((curr_type, curr_chunk))
                    curr_type = 'TEXT'
                    curr_chunk = [line[3:]]
                else:
                    curr_chunk.append(line[3:])
            else:
                if curr_type != 'CODE':
                    if curr_chunk:
                        result.append((curr_type, curr_chunk))
                    curr_type = 'CODE'
                    curr_chunk = [line]
                else:
                    curr_chunk.append(line)
        if curr_chunk:  # Last open chunk.
            result.append((curr_type, curr_chunk))
        return result

    # TODO(holtgrew): Rename into formatCode()?
    def _formatCode(self, txt, linenostart=1):
        try:
            import pygments, pygments.lexers, pygments.formatters
            return pygments.highlight(txt, pygments.lexers.CppLexer(), pygments.formatters.HtmlFormatter(linenos='table', style='friendly', linenostart=linenostart))
        except ImportError:
            return '<pre>' + pyratemp.escape(txt) + '</pre>'

    def highlightCss(self):
        try:
            import pygments.formatters
            return pygments.formatters.HtmlFormatter(linenos='table', style='friendly').get_style_defs('.highlight')
        except ImportError:
            return ''

    def includeFile(self, filename, class_=None, node=None):
        # Get path candidate.
        if not os.path.exists(filename):
            for inc_dir in self.include_dirs:
                if os.path.exists(os.path.join(inc_dir, filename)):
                    filename = os.path.join(inc_dir, filename)
                    break
        # Read in file.
        with open(filename, 'rb') as f:
            fcontents = f.read()
        # Write out file to same directory as filename.
        with open(os.path.join(self.out_path, os.path.basename(filename)), 'wb') as f:
            f.write(fcontents)
        # Create HTML to output.
        chunks = self._splitIncludeFile(fcontents)
        txts = []
        next_lineno = 1
        for type, lines in chunks:
            if type == 'CODE':
                txts.append(self._formatCode('\n'.join(lines), next_lineno))
                next_lineno += len(lines)
            else:
                txts.append('<p class="' + class_ + '">' + self.translateMarkup('\n'.join(lines), node=node) + '</p>')
        
        return '\n'.join(txts)

    def _loadSnippet(self, path, snippet_key):
        result = []
        current_key = None
        current_lines = []
        with open(path, 'rb') as f:
            fcontents = f.read()
        for line in fcontents.splitlines():
            line = line.rstrip()  # Strip line ending and trailing whitespace.
            if line.strip().startswith('//![') and line.strip().endswith(']'):
                key = line.strip()[4:-1].strip()
                if key == current_key:
                    if key == snippet_key:
                        result = current_lines
                    current_lines = []
                    current_key = None
                else:
                    current_key = key
            elif current_key:
                current_lines.append(line)
        return result

    def includeSnippet(self, line, class_=None, node=None):
        filename, snippet_id = line, '<none>'
        if '|' in line:
            filename, snippet_id = line.split('|', 1)
        # Get path candidate.
        if not os.path.exists(filename):
            for inc_dir in self.include_dirs:
                if os.path.exists(os.path.join(inc_dir, filename)):
                    filename = os.path.join(inc_dir, filename)
                    break
        snippet_lines = self._loadSnippet(filename, snippet_id)
        # Read in file.
        with open(filename, 'rb') as f:
            fcontents = f.read()
        # Write out file to same directory as filename.
        with open(os.path.join(self.out_path, os.path.basename(filename)), 'wb') as f:
            f.write(fcontents)
        # Create HTML to output.
        chunks = self._splitIncludeFile(fcontents)
        next_lineno = 1
        txt = self._formatCode('\n'.join(snippet_lines), next_lineno)
        if not class_:
            class_ = ''
        return '<p class="' + class_ + '">' + txt + '</p>'


class DocsCreator(object):
    def __init__(self, error_logger, tree, tpl_path, out_path, include_dirs):
        print >>sys.stderr, 'Setting up Docs Creator'
        self.tree = tree
        self.error_logger = error_logger
        self.tpl_path = tpl_path
        page_tpl_filename = os.path.join(self.tpl_path, 'files', 'page.html')
        ##print page_tpl_filename
        self.page_template = pyratemp.Template(filename=page_tpl_filename)
        index_page_tpl_filename = os.path.join(self.tpl_path, 'files', 'index_page.html')
        self.index_page_template = pyratemp.Template(filename=index_page_tpl_filename)

        self.out_path = os.path.join(out_path, 'files')
        if not os.path.exists(self.out_path):
            os.makedirs(self.out_path)
        self.include_dirs = include_dirs

    def good(self):
        return True  # TODO(holtgrew): Actually look for errors!

    def createIndexPages(self):
        index = core.buildByTypeAndCatIndex(self.tree)
        cat_nodes = self.tree.find('globals.indexes').children
        for cat, node in cat_nodes.iteritems():
            if node.children.has_key('hidefromindex'):
                continue
            print >>sys.stderr, 'Indexes for ' + node.text()
            filename = getIndexPagePath(self.tree, cat, self.out_path)
            with open(filename, 'wb') as f:
                title_node = self.tree.find(['Indexpage', cat, 'title'])
                if title_node:
                    title = title_node.text()
                else:
                    title = self.tree.find(['globals', 'indexes', cat]).text()
                res = self.index_page_template(title=title,
                                               cat=cat,
                                               tree=self.tree,
                                               now=datetime.datetime.utcnow(),
                                               by_subcat=index.get(cat, {}),
                                               time=time,  # for now.strftime()
                                               iterable=lambda x: type(x) is list,
                                               core=core,
                                               html=HtmlHelper(self.error_logger, self.tree, os.path.dirname(filename), self.include_dirs),
                                               json=json)
                f.write(res.encode('utf-8'))

    def copyFiles(self):
        """Copy files in FILE_DIRS."""
        print >>sys.stderr, 'Copying Documentation Files'
        for path in FILE_DIRS:
            entries = os.listdir(path)
            if not os.path.exists(os.path.join(self.out_path, path)):  # Make sure output path exists.
                os.makedirs(os.path.join(self.out_path, path))
            for entry in entries:
                if entry.startswith('.'):
                    continue  # Skip hidden entries.
                source_path = os.path.join(path, entry)
                target_path = os.path.join(self.out_path, path, entry)
                # Copy file.
                print >>sys.stderr, '  %s => %s' % (source_path, target_path)
                shutil.copyfile(source_path, target_path)

    def createPages(self):
        """Create the documentation pages."""
        cat_dict = self.tree.find('globals.categories').children
        for cat, cat_node in cat_dict.iteritems():  # cat=Function,...
            print >>sys.stderr, 'Pages for ' + cat
            if self.tree.find(cat) is None:  # Skip empty categories.
                print >>sys.stderr
                continue
            for subcat, node in self.tree.find(cat).children.iteritems():  # subcat=length,...
                filename = getPagePath(cat, subcat, self.out_path)
                ## print filename
                with open(filename, 'wb') as f:
                    if 'title' in node.children:
                        title = node.children['title'].text()
                    else:
                        title = node.key
                    res = self.page_template(title=title,
                                             cat=cat,
                                             subcat=subcat,
                                             tree=self.tree,
                                             now=datetime.datetime.utcnow(),
                                             time=time,  # for now.strftime()
                                             core=core,
                                             html=HtmlHelper(self.error_logger, self.tree, os.path.dirname(filename), self.include_dirs),
                                             json=json)
                    f.write(res.encode('utf-8'))
                print >>sys.stderr, '.',
            print >>sys.stderr


def createDocs(error_logger, tree, tpl_path, out_path, include_dirs):
    """Fire off the documentation creation."""
    # Things that are created from templates, data JavaScript creation for the
    # client-side JS-driven search.
    creator1 = TplDocsCreator(error_logger, tree, tpl_path, out_path, include_dirs)
    creator1.copyTemplates()
    creator1.createRootIndexPage()
    creator1.createPanelSearchIndex()
    creator1.createPanelTree()
    #return 0
    # Actually create the documentation content; things that will be displayed
    # on the right hand frame.
    creator2 = DocsCreator(error_logger, tree, tpl_path, out_path, include_dirs)
    creator2.copyFiles()
    creator2.createIndexPages()
    creator2.createPages()

    if creator1.good() and creator2.good():
        return 0
    else:
        return 1
