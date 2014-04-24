"""Seqan Doc Links for Trac.

Version 0.1.

Copyright (C) 2010 Manuel Holtgrewe

Install by copying this file into the plugins directory of your trac
work directory.  In your trac.ini, you can use something like this
(the following also shows the defaults).

  [seqan_doc_links]
  prefix = seqan
  base_url = http://www.seqan.de/dddoc/html/
  dox_prefix = dox
  dox_base_url = http://docs.seqan.de/dev3/

Use something like this to test the plugin:

  * {{{[seqan:Page.Sequences]}}} [seqan:Page.Sequences]
  * {{{seqan:Class.Finder}}} seqan:Class.Finder
  * {{{seqan:"Concept.Simple Type"}}} seqan:"Concept.Simple Type"
  * {{{seqan:"Spec.Chunk Pool Allocator}}} seqan:"Spec.Chunk Pool Allocator"
  *
  * {{{dox:ContainerConcept#length}}}
"""
import urllib
import sys

from trac.core import *
import trac.wiki
import genshi.builder as gb
import genshi
from trac.web.chrome import ITemplateProvider, add_stylesheet

def getFilename(cat, item):
    """Get the filename that dddoc would create.

    Args:
      cat   String, category of the link.
      item  String, name of the item.

    Returns:
      File name of the categorized item.
    """
    return cat.upper() + escapeFiles(item) + ".html"


def escapeFiles(text):
    """Escape the file name as dddoc would do it.

    Args:
      text  String with the text to escape.

    Returns:
      Escaped text.
    """
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

    if (len(ret) == 0) or (ret[0] == '_'): return ret
    else: return '.'+ret


class SeqanDocsSyntaxProvider(trac.core.Component):
    """Expands seqan:<Category>.<EntryName> links."""
    trac.core.implements(trac.wiki.IWikiSyntaxProvider)
    implements(ITemplateProvider)

    SECTION_NAME = 'seqan_doc_links'
    DEFAULT_PREFIX = 'seqan'
    DEFAULT_BASE_URL = 'http://www.seqan.de/dddoc/html/'

    def __init__(self):
        # Set defaults.
        self.prefix = self.DEFAULT_PREFIX
        self.base_url = self.DEFAULT_BASE_URL
        # Parse configuration from trac.ini config file.
        for option in self.config.options(self.SECTION_NAME):
            if option[0] == 'prefix':
                self.prefix = option[1]
            if option[0] == 'base_url':
                self.base_url = option[1]

    def get_wiki_syntax(self):
        """Method from IWikiSyntaxProvider.

        Returns empty list, we do not implement any."""
        return []

    def get_link_resolvers(self):
        """Method from IWikiSyntaxProvider.

        Returns iterable (list) of (prefix, function) pairs.
        """
        return [(self.prefix, self.format_doc_link)]

    def format_doc_link(self, formatter, ns, target, label):
        """Function to perform formatting for seqan:XYZ links.

        This roughly follows [1].

        [1] http://trac.edgewall.org/wiki/TracDev/IWikiSyntaxProviderExample
        """
        add_stylesheet(formatter.req, 'doc_links/css/doc_links.css')
        # The following is a heuristic for "no alternative label".
        if ns in label and target in label:
          if '.' in target:
            category, item = tuple(target.split('.', 1))
            label = item
            # Strip everything before and including the first hash.
            if '#' in label:
              label = label.split('#', 1)[1]
          else:
            label = target
        # Ignore if the target does not contain a dot.
        if not '.' in target:
          return target
        # Now, use dddoc's logic to generate the appropriate file name for
        file_name = getFilename(*target.split('.', 1))
        span = [gb.tag.span(genshi.HTML('&nbsp;'), class_='icon'), label]
        title = ' "%s" in SeqAn documentation.' % target
        return gb.tag.a(span, class_='doc-link',
                        href=self.base_url + file_name, title=title)

    ### ITemplateProvider methods

    def get_templates_dirs(self):
        return []

    def get_htdocs_dirs(self):
        from pkg_resources import resource_filename
        return [('doc_links', resource_filename(__name__, 'htdocs'))]


class SeqanDosSyntaxProvider(trac.core.Component):
    """Expands dox:<entry-id> links."""
    trac.core.implements(trac.wiki.IWikiSyntaxProvider)
    implements(ITemplateProvider)

    SECTION_NAME = 'seqan_doc_links'
    DEFAULT_PREFIX = 'dox'
    DEFAULT_BASE_URL = 'http://www.seqan.de/dddoc/html/'

    def __init__(self):
        # Set defaults.
        self.prefix = self.DEFAULT_PREFIX
        self.base_url = self.DEFAULT_BASE_URL
        # Parse configuration from trac.ini config file.
        for option in self.config.options(self.SECTION_NAME):
            if option[0] == 'dox_prefix':
                self.prefix = option[1]
            if option[0] == 'dox_base_url':
                self.base_url = option[1]

    def get_wiki_syntax(self):
        """Method from IWikiSyntaxProvider.

        Returns empty list, we do not implement any."""
        return []

    def get_link_resolvers(self):
        """Method from IWikiSyntaxProvider.

        Returns iterable (list) of (prefix, function) pairs.
        """
        return [(self.prefix, self.format_doc_link)]

    def format_doc_link(self, formatter, ns, target, label):
        """Function to perform formatting for dox:XYZ links.

        This roughly follows [1].

        [1] http://trac.edgewall.org/wiki/TracDev/IWikiSyntaxProviderExample
        """
        # Stylesheet already done for doc_links.
        add_stylesheet(formatter.req, 'doc_links/css/doc_links.css')
        # The following is a heuristic for "no alternative label".
        if not label:
            label = target
        if label.startswith(self.prefix + ':'):
            label = label[len(self.prefix) + 1:]
        # Now, use dddoc's logic to generate the appropriate file name for
        query = '?p=%s' % target  # TODO(holtgrew): url encode
        span = [gb.tag.span(genshi.HTML('&nbsp;'), class_='icon'), label]
        title = ' "%s" in SeqAn documentation.' % target
        return gb.tag.a(span, class_='doc-link',
                        href=self.base_url + query, title=title)

    ### ITemplateProvider methods

    def get_templates_dirs(self):
        return []

    def get_htdocs_dirs(self):
        from pkg_resources import resource_filename
        return [('doc_links', resource_filename(__name__, 'htdocs'))]

