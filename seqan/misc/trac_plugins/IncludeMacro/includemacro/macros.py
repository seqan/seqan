# TracIncludeMacro macros
import re
import urllib2
from StringIO import StringIO

from trac.core import *
from trac.wiki.macros import WikiMacroBase
from trac.wiki.formatter import system_message
from trac.wiki.model import WikiPage
from trac.mimeview.api import Mimeview, get_mimetype, Context
from trac.perm import IPermissionRequestor
from genshi.core import escape
from genshi.input import HTMLParser, ParseError
from genshi.filters.html import HTMLSanitizer

__all__ = ['IncludeMacro']

class IncludeMacro(WikiMacroBase):
    """A macro to include other resources in wiki pages.
    More documentation to follow.
    """
    
    implements(IPermissionRequestor)
    
    # Default output formats for sources that need them
    default_formats = {
        'wiki': 'text/x-trac-wiki',
    }
    
    # IWikiMacroProvider methods
    def expand_macro(self, formatter, name, content):
        req = formatter.req   # Shortcut.
        safe_content = False  # Whether or not to disable cleaning HTML.
        args = [x.strip() for x in content.split(',')]
        if len(args) == 1:
            args.append(None)
        elif len(args) == 3:
            return system_message('args == %s' % args)
            if not args[2].startswith('fragment='):
                msg = ('If three arguments are given, the last one must'
                       ' start with fragment=, but tag content was %s')
                return system_message(msg % content)
        elif len(args) != 2:
            return system_message('Invalid arguments "%s"'%content)

        # Parse out fragment name.
        fragment_name = None
        if args[-1] and args[-1].startswith('fragment='):
            fragment_name = args[-1][len('fragment='):]
            args.pop()
        if len(args) == 1:
            args.append(None)
            
        # Pull out the arguments
        source, dest_format = args
        try:
            source_format, source_obj = source.split(':', 1)
        except ValueError: # If no : is present, assume its a wiki page
            source_format, source_obj = 'wiki', source
            
        # Apply a default format if needed
        if dest_format is None:
            try:
                dest_format = self.default_formats[source_format]
            except KeyError:
                pass
        
        if source_format in ('http', 'https', 'ftp'):
            # Since I can't really do recursion checking, and because this 
            # could be a source of abuse allow selectively blocking it.
            # RFE: Allow blacklist/whitelist patterns for URLS. <NPK>
            # RFE: Track page edits and prevent unauthorized users from ever entering a URL include. <NPK>
            if not req.perm.has_permission('INCLUDE_URL'):
                self.log.info('IncludeMacro: Blocking attempt by %s to include URL %s on page %s', req.authname, source, req.path_info)
                return ''
            try:
                urlf = urllib2.urlopen(source)
                out = urlf.read()  
            except urllib2.URLError, e:
                return system_message('Error while retrieving file', str(e))
            except TracError, e:
                return system_message('Error while previewing', str(e))
            ctxt = Context.from_request(req)
        elif source_format == 'wiki':
            # XXX: Check for recursion in page includes. <NPK>
            if not req.perm.has_permission('WIKI_VIEW'):
                return ''
            page = WikiPage(self.env, source_obj)
            if not page.exists:
                return system_message('Wiki page %s does not exist'%source_obj)
            out = page.text
            ctxt = Context.from_request(req, 'wiki', source_obj)
        elif source_format == 'source':
            if not req.perm.has_permission('FILE_VIEW'):
                return ''
            repo = self.env.get_repository(authname=req.authname)
            node = repo.get_node(source_obj)
            out = node.get_content().read()
            if dest_format is None:
                dest_format = node.content_type or get_mimetype(source_obj, out)
            ctxt = Context.from_request(req, 'source', source_obj)
        # RFE: Add ticket: and comment: sources. <NPK>
        # RFE: Add attachment: source. <NPK>
        else:
            return system_message('Unsupported include source %s'%source)

        # If there was a fragment name given then find the fragment.
        fragment = []
        current_fragment_name = None
        if fragment_name:
            for line in out.splitlines():
                res = re.search(r'FRAGMENT\(([^)]*)\)', line)
                if res:
                    current_fragment_name = res.groups()[0]
                else:
                    if current_fragment_name == fragment_name:
                        fragment.append(line)
            out = '\n'.join(fragment)
            
        # If we have a preview format, use it
        if dest_format:
            # We can trust the output and do not need to call the HTML sanitizer
            # below.  The HTML sanitization leads to whitespace being stripped.
            safe_content = True
            out = Mimeview(self.env).render(ctxt, dest_format, out, force_source=True)
        
        # Escape if needed
        if not safe_content and not self.config.getbool('wiki', 'render_unsafe_content', False):
            try:
                out = HTMLParser(StringIO(out)).parse() | HTMLSanitizer()
            except ParseError:
                out = escape(out)
        
        return out
            
    # IPermissionRequestor methods
    def get_permission_actions(self):
        yield 'INCLUDE_URL'
            
        
        
