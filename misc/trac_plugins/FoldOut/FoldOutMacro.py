# encoding: utf-8
"""Trac Fold Out Macro.

This macro allows to hide sections of text until a "show" link is clicked.

Example:

{{{
#!FoldOut
This is the short description.
It goes to the first horizontal rule.
----
This is the long description.

Both, in the summary and the long description, all wiki text can be used.
}}}
"""

import trac.core
import trac.wiki
import trac.wiki.macros

from genshi.builder import tag
import genshi.core
import uuid
import StringIO

ARROW_RIGHT = u'\u25B6'
ARROW_DOWN = u'\u25BC'

class FoldOutMacro(trac.wiki.macros.WikiMacroBase):
  def expand_macro(self, formatter, name, args):
    summary = []  # Lines with the summary, above the horizontal rule.
    body = []  # Lines with the body, below the horizontal rule.
    sawRule = False  # Flag: Has seen a horizontal rule line.
    # Iterate over the lines in args and split at first line with a horizontal
    # rule.
    for line in args.splitlines():
      if line.startswith('----'):
        sawRule = True
        continue
      if sawRule:
        body.append(line)
      else:
        summary.append(line)
    # Build HTML with summary and toggle'able body.
    hidden_id = uuid.uuid4()
    body_html = self.format_wiki(formatter, '\n'.join(body))
    hidden = tag.div(genshi.core.Markup(body_html), id=hidden_id, style='display:none;')
    toggle_class = uuid.uuid4()
    toggle_js = genshi.core.Markup(u'$(\'#%s\').toggle();$(\'.%s\').toggle();return false;') % (hidden_id, toggle_class)
    toggle_link = tag.a(tag.span(ARROW_RIGHT + ' more...', class_=toggle_class) + tag.span(ARROW_DOWN + ' less...', class_=toggle_class, style='display:none;'), onclick=toggle_js, href='#')
    summary_html = self.format_wiki(formatter, '\n'.join(summary))
    return genshi.core.Markup(summary_html) + toggle_link + genshi.core.Markup(hidden)

  def format_wiki(self, formatter, wiki_string):
    """Format the given string wiki_string to HTML."""
    out = StringIO.StringIO()
    trac.wiki.Formatter(self.env, formatter.context).format(wiki_string, out)
    return out.getvalue()
