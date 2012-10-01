# encoding: utf-8
"""Trac text boxes macro.

This macro allows to create boxes for shell and IDE output.

Example:

{{{
{{{
#!ShellBox
# cd /tmp
# touch foo
}}}

{{{
#!IdeBox
Success!
}}}
}}}
"""

import itertools
import operator
import StringIO

from pkg_resources import resource_filename

from trac.core import *
from trac.web.chrome import ITemplateProvider, add_stylesheet
from trac.web.api import IRequestFilter
from trac.wiki.macros import WikiMacroBase

import trac.wiki
import genshi.core
from genshi.builder import tag

CLASSES = {
    'ShellBox' : 'box shell',
    'IdeBox' : 'box ide',
    'WarningBox': 'box warning',
    'InfoBox': 'box info',
    'ImportantBox': 'box important',
    'AssignmentBox': 'box assignment'
}

ICONS = {
    'WarningBox': 'dialog-warning',
    'InfoBox': 'dialog-information',
    'ImportantBox': 'emblem-important',
    'AssignmentBox': 'Accessories-text-editor',
}

class TextBoxMacro(WikiMacroBase):
    implements(ITemplateProvider, IRequestFilter)
    
    def get_macros(self):
        yield 'ShellBox'
        yield 'IdeBox'
        yield 'WarningBox'
        yield 'InfoBox'
        yield 'ImportantBox'
        yield 'AssignmentBox'
        yield 'TextIcon'
        yield 'MenuTrace'

    def expand_macro(self, formatter, name, content, args):
        add_stylesheet(formatter.req, 'text_boxes/css/text_boxes.css')
        # TODO(holtgrew): Actually, we would like to add a style sheet but this does not work. Thus we use a style tag below.
        #add_stylesheet(formatter.req, 'text_boxes/css/text_boxes.css')
        className = CLASSES.get(name, 'ShellBox')
        if name in ['TextIcon']:
            args = trac.wiki.parse_args(content)
            if not args or not args[0]:  # Handle case of empty entry.
                return None
            DEFAULT = 'gray'
            COLOR_CLASSES = ['green', 'red', 'yellow', 'black', 'gray', 'blue']
            color_class_name = args[1].get('color', DEFAULT)
            if not color_class_name in COLOR_CLASSES:
                color_class_name = DEFAULT
            return tag.span(args[0], class_=('text_icon %s' % color_class_name))
        elif name in ['MenuTrace']:
            args = trac.wiki.parse_args(content)
            if not args[0]:
                return None
            result = tag.span(args[0][0], class_='menu_item')
            for text in args[0][1:]:
                result += tag.span(u' \u25B8 ', class_='arrow')
                result += tag.span(text, class_='menu_item')
            return tag.span(result, class_='menu_trace')
        elif name in ['WarningBox', 'InfoBox', 'ImportantBox', 'AssignmentBox']:
            content_html = self.format_wiki(formatter, content)
            return tag.div(genshi.core.Markup(content_html), class_=className)
        else:
            return tag.pre(content, class_='wiki ' + className)

    def format_wiki(self, formatter, wiki_string):
        """Format the given string wiki_string to HTML."""
        out = StringIO.StringIO()
        trac.wiki.Formatter(self.env, formatter.context).format(wiki_string, out)
        return out.getvalue()

    ### IRequestFilter methods
    
    def pre_process_request(self, req, handler):
        return handler

    def post_process_request(self, req, template, data, content_type):
        add_stylesheet(req, 'text_boxes/css/text_boxes.css')
        return template, data, content_type

    ### ITemplateProvider methods

    def get_templates_dirs(self):
        return []

    def get_htdocs_dirs(self):
        from pkg_resources import resource_filename
        return [('text_boxes', resource_filename(__name__, 'htdocs'))]

