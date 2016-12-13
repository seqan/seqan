#!/usr/bin/env python2
"""Processed version of the documentation.

The documentation from the objects of raw_doc is further processed into
objects from the module proc_doc.  These objects can then be processed
into structured documents such as HTML more easily.
"""

# TODO(holtgrew): Location traceability for entries and text.

import os.path
import HTMLParser
import logging
import re
import sys
import xml.etree.ElementTree
import xml.sax.saxutils

import inc_mgr
import sig_parser
import dox_parser
import dox_tokens
import raw_doc
import validation


def escapeForXml(s):
    """Return escaped XML of s."""
    return xml.sax.saxutils.escape(s)


class DocumentationBuildException(dox_parser.ParserError):
    """Thrown when there is a logical error on building the documentation."""


class LinkResolver(object):
    """Member of ProcDoc for resolving links."""
    
    def __init__(self, proc_doc):
        self.proc_doc = proc_doc


def splitSecondLevelEntry(name):
    """Split second-level entry and return (first, second) pair.

    If name is not a second level entry then (None, name) is returned.
    """
    xs = None
    if name.count('::') > 1 and ' ' in name:
        xs = name.split(' ', 1)
    elif '#' in name:
        xs = name.split('#', 1)
    elif '::' in name:
        xs = name.rsplit('::', 1)
    if xs:
        return xs
    return (None, name)


class ProcDoc(object):
    """Collection of the top-level documentation entries.

    @ivar doc_processor      The DocProcessor that created this ProcDoc.
    @ivar local_name_counter Number of occurrences for local names, used for
                             shortening @link display to second level entries.
    """

    def __init__(self, doc_processor):
        self.doc_processor = doc_processor
        self.top_level_entries = {}
        self.second_level_entries = {}
        self.entries = {}
        self.local_name_counter = {}
        self.top_level_entries_filename = {}
        self.second_level_entries_filename = {}

    def addTopLevelEntry(self, x):
        """Add a top-level-entry."""
        self.registerEntry(x)
        self.top_level_entries[x.name] = x

    def addSecondLevelEntry(self, x):
        """Add a second-level entry."""
        self.registerEntry(x)
        self.second_level_entries[x.name] = x
        first, second = splitSecondLevelEntry(x.name)
        if not first in self.top_level_entries:
            msg = 'Unknown type %s' % first
            token = x.raw_entry.name.tokens[0]
            raise dox_parser.ParserError(msg=msg, token=token)
        if first:
#            print '%s => %s as %s' % (x.name, second, x.kind)
            self.top_level_entries[first].registerSubentry(x)
        # update local name counter
        self.local_name_counter.setdefault(second, 0)
        self.local_name_counter[second] += 1
        
    def addVariable(self, x):
        """Add a variable entry."""
        self.registerEntry(x)
        if '::' in x.name:
            self.second_level_entries[x.name] = x
            first, second = splitSecondLevelEntry(x.name)
            if not first in self.top_level_entries:
                token = x.raw_entry.name.tokens[0]
                self.doc_processor.msg_printer.printTokenError(
                    token, 'Unknown top level entry %s' % first, 'error')
            else:
                self.top_level_entries[first].registerSubentry(x)
        else:
            self.top_level_entries[x.name] = x

    def addEnumValue(self, x):
        """Add an enum value entry."""
        self.registerEntry(x)
        if not x.type in self.top_level_entries:
            token = x.raw_entry.name.tokens[0]
            self.doc_processor.msg_printer.printTokenError(
                token, 'Unknown top level entry %s' % x.type, 'error')
        else:
            self.top_level_entries[x.type].registerSubentry(x)
        
    def registerEntry(self, x):
        """Register an entry."""

        name = x.name
        if name == '':
            msg = 'Entry must not have an empty name.'
            raise DocumentationBuildException(token=x.raw_entry.first_token,
                                              msg=msg)
        if x.kind == 'variable' and self.top_level_entries.get(x.type) and self.top_level_entries[x.type].kind == 'enum':
            name = x.type + '::' + name
        if name.endswith(';'):
            name = name[:-1]
        if name in self.entries:
            old = self.entries[name]
            tpl = ('Trying to define %(new_kind)s %(new_name)s in %(new_file)s:'
                   '%(new_line)s but is previously defined as %(old_kind)s '
                   '%(old_name)s in %(old_file)s:%(old_line)d.')
            vals = {
                'new_kind': x.kind,
                'new_name': name,
                'new_file': old.location[0],
                'new_line': old.location[1],
                'old_kind': old.kind,
                'old_name': old.name,
                'old_file' : x.location[0],
                'old_line' : x.location[1]}
            self.doc_processor.msg_printer.printTokenError(x.raw_entry.name.tokens[0],
                                                           tpl % vals, 'error')
        else:
          self.entries[name] = x
          x.doc = self

    def runTextVisitor(self, v):
        """Run visitor v on all Text members of all entries and sub entries.
        """
        for e in self.entries.itervalues():
            e.runTextVisitor(v)


class TextNode(object):
    """A node represents a part of a processed text.

    Processed text is text generated from tokens lexed from the input file.
    For example, text in the paragraph of a entry's body can be representd by
    TextNode objects.

    TextNode objects are similar to DOM nodes, i.e. they can contain children
    and have attributes.  This means that we can have a link node that has a
    href/target attribute with a target URL and one or more child nodes that
    contain the link's label.

    Additionally, we store the source location (begin and end line/column) of
    the node in its source file.

    We represent plain links, i.e. where the label is the same as the target
    using the representation for "<a href="seqan:$target">$target</a>".

    We represent included code snippets as "<dox:code type='.cpp'>$code</dox:code>."

    @ivar type: The type of the node, as a string.  Reserved values are
                '<text>' for plain text nodes.
    @ivar attrs: A dict object mapping attribute names to string values.
    @ivar children: A list of TextNode objects.
    @ivar text: The text value of a node, a string.
    @ivar tokens: For links, this is the list of tokens in the @link command.
    @ivar raw_html: Comes from HTML in dox, disables h2 to h4 translation, for
                    example.
    """

    def __init__(self, type='<text>', raw_html=False, verbatim=False, text='', attrs={}):
        self.type = type
        self.attrs = dict(attrs)
        self.children = []
        self.raw_html = raw_html
        if verbatim:
            self.text = text
        else:
            self.text = text.replace('<', '&lt;').replace('>', '&gt;')

    def __str__(self):
        attrs = (repr(self.type), repr(self.text), repr(self.attrs), len(self.children))
        return 'TextNode(type=%s, text=%s, attrs=%s, len(children)=%d)' % attrs

    def __repr__(self):
        return str(self)

    @property
    def empty(self):
        if self.type == '<text>':
            return not not self.text
        else:
            return not (self.children or self.attrs)

    def setAttr(self, key, value):
        self.attrs[escapeForXml(key)] = escapeForXml(value)

    def addChild(self, n):
        self.children.append(n)
        return self.children[-1]

    @property
    def X(self):
        """Returns first child, used to retrieve member of top-level <div>."""
        if self.type == '<text>':
            return self
        else:
            return self.children[0]

    @property
    def plainText(self):
        """Converts to HTML and strips tags."""
        def remove_html_markup(s):
            tag = False
            quote = False
            out = []

            for c in s:
                if c == '<' and not quote:
                    tag = True
                elif c == '>' and not quote:
                    tag = False
                elif (c == '"' or c == "'") and tag:
                    quote = not quote
                elif not tag:
                    out.append(c)
            return ''.join(out)
        return remove_html_markup(self.toHtmlLike())

    def toHtmlLike(self, skip_top_tag=False, **kwargs):
        """Returns a string with a HTML-like representation for debuggin.

        @param skip_top_tag: Do ont output top-level tag.
        @param kwargs: Additional attributes to add.
        """
        if self.type == '<text>':
            if self.attrs:
                print >>sys.stderr, 'WARNING: Attributes on text node!'
            return self.text
        else:
            dash = {True: '', False: ' /'}.get(bool(self.children))  # Whether empty
            res = []
            if not skip_top_tag:
                res += ['<', self.type]
                for key, value in self.attrs.iteritems():
                    res += [' ', key, '=', '"', repr(value)[1:-1], '"']
                for key, value in kwargs.iteritems():
                    res += [' ', key, '=', '"', value, '"']
                res.append(dash + '>')
            if self.children:
                res += [x.toHtmlLike() for x in self.children]
                if not skip_top_tag:
                    res += ['</', self.type, '>']
            return ''.join(res)


class ProcEntry(object):
    """A processed representation of a documentation entry.

    A documentation entry has a kind (string), a name (string), a brief
    description (TextNode(type='<text>')), and a list of references/sees to
    other elements (list of TextNode(type='<link>')).  Also, it has a body
    which is a TextNode with children.

    @ivar raw: The raw_doc.Raw* object.
    @ivar kind: The kind of the entry, string.
    @ivar name: The name of the entry, string.
    @ivar title_str: A string with the title.
    @ivar brief: A brief description, a text-typed TextNode or None.
    @ivar body: A TextNode object with children for the documentation body.
    @ivar sees: A list of link-typed TextNode objects, can be empty.
    @ivar doc:  The owning, ProcDoc, set on ProcDoc.registerEntry
    @ivar subentries: Sub entries, dir, grouped by type.
    @ivar raw_entry: The RawEntry object that this ProcEntry was generated from.
    """

    def __init__(self, raw, name, title=None, brief=None, body=None, sees=[]):
        self.raw = raw
        self.name = name
        self.title_str = title
        self.brief = brief
        self.body = body
        self.sees = list(sees)
        self.doc = None
        self.subentries = {}
        self.raw_entry = None
        self._location = None

    def sortedSees(self):
        return sorted(self.sees, key=lambda x: x.plainText)

    def registerSubentry(self, proc_entry):
        self.subentries.setdefault(proc_entry.kind, []).append(proc_entry)

    def hasSubEntry(self, kind, proc_doc):
        """Returns has a subentry of the given kind."""
        if self.subentries.get(kind):
            return True
        if hasattr(self, 'all_extended'):
            for cl in self.all_extended:
                extended = proc_doc.top_level_entries[cl]
                if extended.subentries.get(kind):
                    return True
        if hasattr(self, 'all_implemented'):
            for co in self.all_implemented:
                extended = proc_doc.top_level_entries[co]
                if extended.subentries.get(kind):
                    return True
        return False

    def visitTextNodes(self, visitor):
        """Visit all text nodes using the given visitor."""
        visitor.visit(self.brief)
        visitor.visit(self.body)
        for see in self.sees:
            visitor.visit(see)

    @property
    def title(self):
        return self.title_str or self.name

    @property
    def location(self):
        """Returns pair (path, line)."""
        if not self._location:
            path = '<none>'
            line = -1
            if self.raw_entry.name and self.raw_entry.name.tokens:
                line = self.raw_entry.name.tokens[0].lineno
                path = self.raw_entry.name.tokens[0].file_name
            self._location = (path, line)
        return self._location

    @property
    def kind(self):
        return self.__class__.__name__.replace('Proc', '').lower()


class ProcCodeEntry(ProcEntry):
    """A processed code entry.

    @ivar signatures: A TextNode with the signatures of the entry.  They are
                      properly formatted to be displayed as verbatim text.
    @ivar signature_entries: A list of sig_parser.SigEntry objects.
    @ivar headerfiles: A list of str objects with the arguments to #include.
    @ivar deprecation_msgs: List of TextNode objects with deprecation messages.
    @ivar notes: List of TextNode objects with notes.
    @ivar warnings: List of TextNode objects with warnings.
    @ivar akas: List of strings.
    @ivar internals: List of TextNode objects (possibly empty) with marks as
                     internal.
    """

    def __init__(self, raw, name, brief=None, body=None, sees=[]):
        ProcEntry.__init__(self, raw, name, brief, body, sees)
        self.signatures = []
        self.signature_entries = []
        self.headerfiles = []
        self.deprecation_msgs = []
        self.notes = []
        self.warnings = []
        self.akas = []
        self.internals = []

    def addSignature(self, s):
        self.signatures.append(s)

    def addSignatureEntry(self, e):
        self.signature_entries.append(e)

    def addHeaderfile(self, h):
        self.headerfiles.append(h)
        
    def addDeprecationMsg(self, m):
        self.deprecation_msgs.append(m)

    def addNote(self, n):
        self.notes.append(n)
        
    def addWarning(self, w):
        self.warnings.append(w)
        
    def addAkas(self, a):
        self.akas.append(a)
        
    def addInternal(self, i):
        self.internals.append(i)
        
    def subEntries(self, kind):
        return []

    def visitTextNodes(self, visitor):
        """Visit all text nodes using the given visitor."""
        ProcEntry.visitTextNodes(self, visitor)
        for msg in self.deprecation_msgs:
            visitor.visit(msg)
        for msg in self.notes:
            visitor.visit(msg)
        for msg in self.warnings:
            visitor.visit(msg)


class ProcEnum(ProcCodeEntry):
    """A processed enum documentation.

    @ivar values: A list of ProcVariable entries that represent values
                  of this enum.
    """

    def __init__(self, raw, name, brief=None, body=None, sees=[]):
        ProcCodeEntry.__init__(self, raw, name, brief, body, sees)
        self.values = []


class ProcAdaption(ProcCodeEntry):
    """A processed adaption documentation.

    @ivar values: A list of ProcVariable entries that represent values
                  of this adaption.
    """

    def __init__(self, raw, name, brief=None, body=None, sees=[]):
        ProcCodeEntry.__init__(self, raw, name, brief, body, sees)
        self.values = []


class ProcTypedef(ProcCodeEntry):
    """A processed typedef documentation.

    @ivar values: A list of ProcVariable entries that represent values
                  of this typedef.
    """

    def __init__(self, raw, name, brief=None, body=None, sees=[]):
        ProcCodeEntry.__init__(self, raw, name, brief, body, sees)
        self.values = []

    @property
    def kind(self):
        if '#' in self.name:
            return 'grouped_typedef'
        elif '::' in self.name:
            return 'member_typedef'
        else:
            return 'global_typedef'


class ProcConcept(ProcCodeEntry):
    """A processed concept documentation.

    @ivar extends: A list of str values with the names of the extended
                   concepts.
    @ivar all_extended: A set of str values with the names of all extended
                        concepts, also transitively.
    @ivar all_extending: A set of str values with the names of all extending
                         concepts.
    @ivar all_implementing: A set of str values with the names of all
                            implementing classes.
    """

    def __init__(self, raw, name, brief=None, body=None, sees=[]):
        ProcCodeEntry.__init__(self, raw, name, brief, body, sees)
        self.extends = []
        self.all_extended = set()
        self.all_extending = set()
        self.all_implementing = set()

    def addExtends(self, s):
        self.extends.append(s)

    def __str__(self):
        return 'ProcConcept(%s, brief=%s, body=%s, sees=%s)' % (
            repr(self.name), repr(self.brief), repr(self.body), repr(self.sees))

    def __repr__(self):
        return 'ProcConcept(%s)' % repr(self.name)


class ProcClass(ProcCodeEntry):
    """A processed class documentation.

    @ivar extends: A list of str values with the names of the extended
                   classes.
    @ivar implements: A list of str values with the names of the implemented
                      concepts.
    @ivar all_implemented: Set of str values with the names of all implemented
                           concepts.
    @ivar all_extending: Set of str values with the names of all extending
                         classes.
    @ivar all_extended: Set of str values with the names of all extended classes.
    """

    def __init__(self, raw, name, brief=None, body=None, sees=[]):
        ProcCodeEntry.__init__(self, raw, name, brief, body, sees)
        self.extends = []
        self.implements = []
        self.all_implemented = set()
        self.all_extending = set()
        self.all_extended = set()
        self.tparams = []
        self.typedefs = []

    def visitTextNodes(self, visitor):
        """Visit all text nodes using the given visitor."""
        ProcCodeEntry.visitTextNodes(self, visitor)
        for p in self.tparams:
            p.visitTextNodes(visitor)
        for p in self.typedefs:
            p.visitTextNodes(visitor)

    def addExtends(self, s):
        self.extends.append(s)

    def addImplements(self, s):
        self.implements.append(s)

    def addTParam(self, t):
        self.tparams.append(t)

    def addTypedef(self, t):
        self.typedefs.append(t)

    @property
    def isSpecialization(self):
        return not not self.extends

    @property
    def kind(self):
        if self.isSpecialization:
            return 'specialization'
        else:
            return 'class'


class ProcTag(ProcCodeEntry):
    """A processed tag documentation.
    """

    def __init__(self, raw, name, brief=None, body=None, sees=[]):
        ProcCodeEntry.__init__(self, raw, name, brief, body, sees)
        self.tparams = []

    def addTParam(self, t):
        self.tparams.append(t)

    @property
    def local_name(self):
        """Returns name without group prefix."""
        if '#' in self.name:
            return self.name.split('#', 1)[-1]
        else:
            return self.name

    @property
    def kind(self):
        if '#' in self.name:
            return 'grouped_tag'
        else:
            return 'tag'


class ProcParam(object):
    """Representation of a parameter.

    @ivar raw: Raw representation.
    @ivar name: The name of the parameter. str.
    @ivar in_out: One of IN, OUT, IN_OUT, None.
    @ivar desc: Documentation of the parameter. TextNode.
    """

    def __init__(self, raw):
        self.raw = raw
        self.name = None
        self.in_out = None
        self.desc = TextNode()

    def visitTextNodes(self, visitor):
        """Visit all text nodes using the given visitor."""
        visitor.visit(self.desc)


ProcParam.IN = 'IN'
ProcParam.OUT = 'OUT'
ProcParam.IN_OUT = 'IN_OUT'


class ProcTParam(object):
    """Documentation of a processed template parameter.

    @ivar var: The raw representation.
    @ivar type: The type of the parameter. str
    @ivar desc: Documentation of the parameter. TextNode.
    """

    def __init__(self, raw):
        self.raw = raw
        self.type = None
        self.desc = TextNode()

    def visitTextNodes(self, visitor):
        """Visit all text nodes using the given visitor."""
        visitor.visit(self.desc)


class ProcReturn(object):
    """Documentation of a @return entry.

    @ivar raw: The raw version of this ProcReturn (required for location lookup).
    @ivar type: The return type. str.
    @ivar desc: The documentation of the return value. TextNode.
    """

    def __init__(self, raw):
        self.raw = raw
        self.type = None
        self.desc = TextNode()

    def visitTextNodes(self, visitor):
        """Visit all text nodes using the given visitor."""
        visitor.visit(self.desc)


class ProcThrow(object):
    """Documentation of a @throw entry.

    @ivar raw: The raw version of this ProcThrow (required for location lookup).
    @ivar type: The exception type. str.
    @ivar desc: The documentation of the exception. TextNode.
    """

    def __init__(self, raw):
        self.raw = raw
        self.type = None
        self.desc = TextNode()

    def visitTextNodes(self, visitor):
        """Visit all text nodes using the given visitor."""
        visitor.visit(self.desc)
        

class ProcDataRace(object):
    """Documentation of a @datarace entry.

    @ivar raw: The raw version of this ProcDataRace (required for location lookup).
    @ivar desc: The documentation of the @datarace clause. TextNode.
    """

    def __init__(self, raw):
        self.raw = raw
        self.desc = TextNode()

    def visitTextNodes(self, visitor):
        """Visit all text nodes using the given visitor."""
        visitor.visit(self.desc)


class ProcFunction(ProcCodeEntry):
    """A processed function documentation.

    @ivar params: A list of str values with the names of the extended
                  concepts.
    @ivar tparams:
    @ivar returns:
    @ivar throws: List of ProcThrow objects.
    @ivar dataraces: List of ProcDataRace objects.
    """

    def __init__(self, raw, name, brief=None, body=None, sees=[]):
        ProcCodeEntry.__init__(self, raw, name, brief, body, sees)
        self.params = []
        self.tparams = []
        self.returns = []
        self.throws = []
        self.dataraces = []

    @property
    def kind(self):
        if '#' in self.name:
            return 'interface_function'
        elif '::' in self.name:
            return 'member_function'
        else:
            return 'global_function'

    @property
    def local_name(self):
        """Returns name without class/concept prefix."""
        if '#' in self.name:
            return self.name.split('#', 1)[-1]
        elif '::' in self.name:
            return self.name.split('::', 1)[-1]
        else:
            return self.name

    def visitTextNodes(self, visitor):
        """Visit all text nodes using the given visitor."""
        ProcCodeEntry.visitTextNodes(self, visitor)
        for p in self.params:
            p.visitTextNodes(visitor)
        for p in self.tparams:
            p.visitTextNodes(visitor)
        for p in self.returns:
            p.visitTextNodes(visitor)
        for t in self.throws:
            t.visitTextNodes(visitor)
        for d in self.dataraces:
            d.visitTextNodes(visitor)
        
    def addParam(self, p):
        self.params.append(p)

    def addTParam(self, t):
        self.tparams.append(t)

    def addReturn(self, r):
        self.returns.append(r)

    def addThrow(self, t):
        self.throws.append(t)
    
    def addDataRace(self, t):
        self.dataraces.append(t)


class ProcMacro(ProcCodeEntry):
    """A processed macro documentation.

    @ivar params: A list of str values with the names of the extended
                  concepts.
    @ivar returns: Name displayed for return type.
    @ivar throws: List of ProcThrow objects.
    @ivar throws: List of ProcDataRace objects.
    """

    def __init__(self, raw, name, brief=None, body=None, sees=[]):
        ProcCodeEntry.__init__(self, raw, name, brief, body, sees)
        self.params = []
        self.returns = []
        self.throws = []
        self.dataraces = []

    @property
    def local_name(self):
        """Returns name without group prefix."""
        if '#' in self.name:
            return self.name.split('#', 1)[-1]
        else:
            return self.name

    @property
    def kind(self):
        if '#' in self.name:
            return 'grouped_macro'
        else:
            return 'macro'

    def visitTextNodes(self, visitor):
        """Visit all text nodes using the given visitor."""
        ProcCodeEntry.visitTextNodes(self, visitor)
        for p in self.params:
            p.visitTextNodes(visitor)
        for p in self.returns:
            p.visitTextNodes(visitor)
        for t in self.throws:
            t.visitTextNodes(visitor)
        for d in self.dataraces:
            d.visitTextNodes(visitor)

    def addParam(self, p):
        self.params.append(p)

    def addReturn(self, r):
        self.returns.append(r)

    def addThrow(self, t):
        self.throws.append(t)
        
    def addDataRace(self, d):
        self.dataraces.append(d)


class ProcMetafunction(ProcCodeEntry):
    """A processed function documentation.

    @ivar tparams: A list of str values with the names of the extended
                   concepts.
    @ivar returns: A list of ProcReturn values.
    """

    def __init__(self, raw, name, brief=None, body=None, sees=[]):
        ProcCodeEntry.__init__(self, raw, name, brief, body, sees)
        self.tparams = []
        self.returns = []

    @property
    def kind(self):
        if '#' in self.name:
            return 'interface_metafunction'
        else:
            return 'global_metafunction'

    @property
    def local_name(self):
        """Returns name without class/concept prefix."""
        if '#' in self.name:
            return self.name.split('#', 1)[-1]
        elif '::' in self.name:
            return self.name.split('::', 1)[-1]
        else:
            return self.name

    def visitTextNodes(self, visitor):
        """Visit all text nodes using the given visitor."""
        ProcCodeEntry.visitTextNodes(self, visitor)
        for p in self.tparams:
            p.visitTextNodes(visitor)
        for p in self.returns:
            p.visitTextNodes(visitor)

    def addTParam(self, t):
        self.tparams.append(t)

    def addReturn(self, r):
        self.returns.append(r)


class ProcVariable(ProcCodeEntry):
    """A processed variable documentation.

    @ivar type: A string with the name of a type.
    """

    def __init__(self, raw, name, brief=None, body=None, sees=[]):
        ProcCodeEntry.__init__(self, raw, name, brief, body, sees)
        self.type = None

    @property
    def local_name(self):
        """Returns name without class prefix."""
        if '::' in self.name:
            return self.name.split('::', 1)[-1]
        else:
            return self.name

    @property
    def kind(self):
        if '::' in self.name:
            return 'member_variable'
        else:
            return 'variable'


class ProcEnumValue(ProcVariable):
    """A processed enum value documentation."""

    def __init__(self, raw, name, brief=None, body=None, sees=[]):
        ProcVariable.__init__(self, raw, name, brief, body, sees)

    @property
    def local_name(self):
        """Returns name without class prefix."""
        if '::' in self.name:
            return self.name.split('::', 1)[-1]
        else:
            return self.name

    @property
    def kind(self):
        return 'enum_value'


class ProcPage(ProcEntry):
    """A processed page."""

    def __init__(self, raw, name, brief=None, body=None, sees=[]):
        ProcEntry.__init__(self, raw, name, brief, body, sees)

    def __str__(self):
        return 'Page(name=%s)' % repr(self.name)


class ProcGroup(ProcEntry):
    """A processed group."""

    def __init__(self, raw, name, brief=None, body=None, sees=[]):
        ProcEntry.__init__(self, raw, name, brief, body, sees)
        self.tags = []
        self.typedefs = []

    def __str__(self):
        return 'Group(name=%s)' % repr(self.name)

    def addTypedef(self, t):
        self.typedefs.append(t)


class HtmlTagParser(HTMLParser.HTMLParser):
    """Used for parsing HTML and storing the first tag and its attributes."""

    def __init__(self):
        self.reset()
    
    def reset(self):
        HTMLParser.HTMLParser.reset(self)
        self.tag = None
        self.attrs = None
        self.is_open = None
        self.is_close = None

    def handle_starttag(self, tag, attrs):
        self.tag = tag
        self.attrs = dict(attrs)
        self.is_open = True

    def handle_endtag(self, tag):
        self.tag = tag
        self.is_close = True

    def parse(self, txt):
        self.reset()
        self.feed(txt)


class RawTextToTextNodeConverter(object):
    """Convert raw text including HTML tags to text node.

    @ivar doc_proc: The parent DocProcessor.
    """

    def __init__(self, strip_lt_line_space=False, expected_tags=set(), doc_proc=None):
        self.tag_stack = []
        self.node_stack = []
        self.current = None
        self.strip_lt_line_space = strip_lt_line_space
        self.expected_tags = expected_tags
        # Processing text between inline @begintag @endtag is done by first
        # scanning over the tokens and then processing the tokens in between
        # recursively with a new RawTextToTextNodeConverter.
        self.current_cmd = None  # current command, e.g. '@link'
        self.tokens_cmd = []  # currently scanned tokens
        self.commands = ['COMMAND_LINK', 'COMMAND_ENDLINK']
        self.command_pairs = {'COMMAND_LINK': 'COMMAND_ENDLINK'}
        self.html_parser = HtmlTagParser()
        self.doc_proc = doc_proc

    def handleTag(self, token):
        """Handle a HTML tag.

        The HTML tag is translated into a TextNode and appended to self.current.
        Note that this is meant for parsing one tag only.
        """
        self.html_parser.parse(token.val)
        tag_name = self.html_parser.tag
        if tag_name not in self.expected_tags:
            msg = 'Unknown tag "%s".  Expected one of %s.' % (tag_name, self.expected_tags)
            self.doc_proc.msg_printer.printTokenError(token, msg, 'warning')

        if self.html_parser.is_open:  # Opening tag.
            self.tag_stack.append(self.html_parser.tag)
            self.node_stack.append(self.current)
            tag = TextNode(type=self.html_parser.tag, raw_html=True)
            for key, value in self.html_parser.attrs.items():
                tag.setAttr(key, value)
            self.current = self.current.addChild(tag)

        if self.html_parser.is_close:  # No else, also handle standalone tags.
            if self.tag_stack and self.tag_stack[-1] == tag_name:
                self.tag_stack.pop()  # correct closing tag
            elif self.tag_stack and self.tag_stack[-1] != tag_name:
                # incorrect closing, pop and return
                args = (tag_name, self.tag_stack[-1])
                self.doc_proc.msg_printer.printTokenError(token, 'Closing wrong tag %s instead of %s' % args, 'warning')
                self.tag_stack.pop()
                return
            else:  # not self.tag_stack
                self.doc_proc.msg_printer.printTokenError(token, 'Closing tag without opening %s!' % tag_name, 'warning')
            # Pop from node stack.
            if self.node_stack:
                self.current = self.node_stack[-1]
                self.node_stack.pop()
            else:
                self.doc_proc.msg_printer.printTokenError(token, 'Having closed too many tags!', 'warning')

    def handleCommand(self, token):
        """Handle command for the given token."""
        if self.current_cmd:  # There is a command active
            if token.type == self.command_pairs[self.current_cmd]:  # closing current
                self.handleCommandClosing()  # handle closing of command
            else:  # not closing current
                self.tokens_cmd.append(token)
        else:  # no command active, open
            if not token.type in self.command_pairs.keys():
                expected = map(dox_tokens.transToken, self.command_pairs.keys())
                raise dox_parser.ParserError(token, 'Unexpected token.  Must be one of %s' % expected)
            self.current_cmd = token.type

    def handleCommandClosing(self):
        """Handle closing of current command."""
        assert self.current_cmd == 'COMMAND_LINK', 'Only known command.'
        if self.current_cmd == 'COMMAND_LINK':
            # Trim leading/trailing whitespace tokens
            def isWhitespace(t):
                return t.type in dox_tokens.WHITESPACE
            while self.tokens_cmd and isWhitespace(self.tokens_cmd[0]):
                self.tokens_cmd.pop(0)
            while self.tokens_cmd and isWhitespace(self.tokens_cmd[-1]):
                self.tokens_cmd.pop(-1)
            if not self.tokens_cmd:
                print >>sys.stderr, 'WARNING: Empty @link @endlink.'
                return
            # Get link target.
            target_tokens = []
            while self.tokens_cmd and self.tokens_cmd[0].type not in dox_tokens.WHITESPACE:
                target_tokens.append(self.tokens_cmd.pop(0))
            # Trim leading whitespace again.
            while self.tokens_cmd and isWhitespace(self.tokens_cmd[0]):
                self.tokens_cmd.pop(0)
            # Translate any remaining non-whitespace tokens.
            title_tokens = self.tokens_cmd or list(target_tokens)
            link_text = raw_doc.RawText(title_tokens)
            conv = RawTextToTextNodeConverter(expected_tags=self.expected_tags, doc_proc=self.doc_proc)
            link_text_node = conv.run(link_text)
            link_text_node.type = 'a'
            link_text_node.attrs = {'href': 'seqan:' + ''.join([t.val for t in target_tokens])}
            link_text_node.tokens = target_tokens
            self.current.addChild(link_text_node)
        self.tokens_cmd = []
        self.current_cmd = None

    def run(self, raw_text, verbatim=False):
        """Convert the tokens in raw_text into a tree of TextNode objects."""
        #print >>sys.stderr, '================== %s' % raw_text.text
        #print >>sys.stderr, [(t.type, t.val) for t in raw_text.tokens]
        self.current = TextNode(type='div')
        root = self.current
        at_line_start = True
        for i, t in enumerate(self.fixEntityTokens(raw_text.tokens, verbatim)):
            if self.current_cmd:  # collect token in self.tokens_cmd
                self.handleCommand(t)
                continue

            if t.type in dox_tokens.WHITESPACE:
                if i == 0 or (i + 1) == len(raw_text.tokens):
                    continue  # Ignore leading and trailing whitespace.
                if t.type == 'SPACE' and at_line_start:
                    continue  # Ignore space at the beginning of a line.
                if t.type == 'BREAK':
                    self.current.addChild(TextNode(text='\n'))
                else:
                    self.current.addChild(TextNode(text=' '))
            elif not verbatim and t.type == 'HTML_TAG':
                at_line_start = False
                self.handleTag(t)
            elif not verbatim and t.type in self.commands:
                #print >>sys.stderr, 'command %s' % t
                at_line_start = False
                self.handleCommand(t)
            else:
                at_line_start = False
                # TODO(holtgrew): Escape values.
                self.current.addChild(TextNode(text=t.val))
            at_line_start = t.type in ['EMPTY_LINE', 'BREAK']
        if self.current_cmd:
            self.doc_proc.msg_printer.printTokenError(t, 'Open command %s!' % self.current_cmd, 'warning')
        return root

    def fixEntityTokens(self, tokens, verbatim):
        """Fix entities on a list of tokens.

        Ampersands opening entities that are not closed (i.e. EOF or space are
        hit) are replaced by '&amp;'.
        """
        if verbatim:
            return tokens
        res = []
        buf = []  # Buffer when an ampersand was opened
        in_entity = False
        for token in tokens:
            if token.val == '&':  # type is PUNCTUATION
                in_entity = True
                buf.append(token)
            elif any(c.isspace() for c in token.val):
                if in_entity:
                    buf[0].val = '&amp;'
                    self.doc_proc.msg_printer.printTokenError(buf[0], 'Unclosed entity!', 'warning')
                res += buf
                res.append(token)
                in_entity = False
                buf = []
            elif token.val == ';':  # type is PUNCTUATION
                in_entity = False
                res += buf
                res.append(token)
                buf = []
            elif in_entity:
                buf.append(token)
            else:
                res.append(token)
        if in_entity:
            buf[0].val = '&amp;'
            self.doc_proc.msg_printer.printTokenError(buf[0], 'Unclosed entity!', 'warning')
            res += buf
        return res

    def process(self, raw_entry):
        raise Exception('Not implemented!')


class EntryConverter(object):
    """Base class for the conversion of raw entries processed entries.

    @ivar doc_proc: DocProcessor object.
    @ivar entry_class: The class of the ProcEntry type to create.
    """

    def __init__(self, doc_proc):
        self.doc_proc = doc_proc
        self.entry_class = None
    
    def rawTextToTextNode(self, raw_text, strip_lt_line_space=False, verbatim=False):
        """Convert RawText object into a TextNode object.

        The text node will have the type 'div'.

        @param strip_lt_breaks_lines: Whether or not to remove leading
                                      space for lines.
        @param verbatim: Whether or not to convert HTML tags.
        """
        converter = RawTextToTextNodeConverter(
            strip_lt_line_space, expected_tags=self.doc_proc.expected_tags,
            doc_proc=self.doc_proc)
        return converter.run(raw_text, verbatim)

    def bodyToTextNode(self, raw_body):
        """Convert a RawBody to a TextNode."""
        res = TextNode(type='div')
        for p in raw_body.paragraphs:
            try:
                if p.getType() == 'paragraph':
                    if not p.text.text.strip():
                        continue  # Skip whitespace
                    p = self.rawTextToTextNode(p.text)
                    p.type = 'p'
                    res.addChild(p)
                elif p.getType() == 'section':
                    h = self.rawTextToTextNode(p.heading)
                    h.type = 'h%d' % (p.level + 1)
                    res.addChild(h)
                elif p.getType() == 'include':
                    # Including a whole file.
                    ftype = os.path.splitext(p.path.text)[1]
                    code_text = self.doc_proc.include_mgr.loadFile(p.path.text)
                    proc_include = TextNode(type='dox:code', attrs={'type': ftype, 'source': 'include', 'path': p.path.text})
                    proc_include.addChild(TextNode(text=code_text, verbatim=True))
                    res.addChild(proc_include)
                elif p.getType() == 'snippet':
                    # Including a snippet file.
                    ftype = os.path.splitext(p.path.text)[1]
                    code_text = self.doc_proc.include_mgr.loadSnippet(p.path.text, p.name.text)
                    proc_snippet = TextNode(type='dox:code', attrs={'type': ftype, 'source': 'snippet', 'path': p.path.text})
                    proc_snippet.addChild(TextNode(text=code_text, verbatim=True))
                    res.addChild(proc_snippet)
                elif p.getType() == 'code':
                    code_text = p.text.text
                    type = '.txt'
                    m = re.match(r'^{[^}]+}', code_text)
                    if m:
                        type = m.group(0)[1:-1]
                    code_text = code_text[len(type) + 2:].strip()
                    #print [repr(t.val) for t in p.text.tokens]
                    x = TextNode(type='dox:code', attrs={'type': type})
                    x.addChild(TextNode(text=code_text, verbatim=True))
                    res.addChild(x)
                elif p.getType() == 'htmlonly':
                    res.addChild(TextNode(text=p.text.text, verbatim=True))
            except inc_mgr.IncludeException, e:
                e2 = dox_parser.ParserError(msg=str(e), token=p.text.tokens[0])
                self.doc_proc.msg_printer.printParserError(e2)
                n = TextNode(type='div', attrs={'class': 'note warning'})
                n.children.append(TextNode(text=str(e)))
                res.addChild(n)
        return res

    def process(self, raw_entry):
        entry = self.entry_class(raw_entry, name=raw_entry.name.text)
        # Convert the title
        if raw_entry.title.text:
            entry.title_str = raw_entry.title.text
        # Convert first brief member.  We already warned about duplicate ones
        # elsewhere.
        if raw_entry.briefs:
            entry.brief = self.rawTextToTextNode(raw_entry.briefs[0].text)
        # Convert the body
        if raw_entry.body:
            entry.body = self.bodyToTextNode(raw_entry.body)
        # Convert the sees entries.
        for see in raw_entry.sees:
            link = self.rawTextToTextNode(see.text)
            if see.text.text.startswith('http'):
                link.type = 'a'
                link.attrs['href'] = see.text.text
                link.attrs['target'] = '_top'
            else:
                link = self.rawTextToTextNode(see.text)
                link.type = 'a'
                link.attrs['href'] = 'seqan:%s' % see.text.text
            link.tokens = list(see.text.tokens)
            entry.sees.append(link)
        # Store the raw entry in the processed ones.
        entry.raw_entry = raw_entry
        return entry

    
class CodeEntryConverter(EntryConverter):
    """Base for the processing RawCodeEntry objects into processed entries."""

    def __init__(self, doc_proc):
        EntryConverter.__init__(self, doc_proc)
        self.parse_signature = True
    
    def process(self, raw_entry):
        entry = EntryConverter.process(self, raw_entry)
        # Add headerfile paths as list of strings.
        for s in raw_entry.headerfiles:
            entry.addHeaderfile(s.text.text.strip())
        # Add deprecation messages, notes, warnings, and internal markers as list of TextNodes.
        for s in raw_entry.deprecation_msgs:
            entry.addDeprecationMsg(self.rawTextToTextNode(s.text, strip_lt_line_space=True))
        for s in raw_entry.notes:
            entry.addNote(self.rawTextToTextNode(s.text, strip_lt_line_space=True))
        for s in raw_entry.warnings:
            entry.addWarning(self.rawTextToTextNode(s.text, strip_lt_line_space=True))
        for s in raw_entry.internals:
            entry.addInternal(self.rawTextToTextNode(s.text, strip_lt_line_space=True))
        # Add aka messages as strings.
        for s in raw_entry.akas:
            entry.addAkas(s.text.text.strip())
        # Add signatures as a text node with code.
        for s in raw_entry.signatures:
            entry.addSignature(self.rawTextToTextNode(s.text, strip_lt_line_space=True,
                                                      verbatim=True))
        # Use sig_parser to convert the signature texts to SigEntry objects.
        # They are used for the list of functions/metafunctions for a type.
        if self.parse_signature:
            for s in raw_entry.signatures:
                try:
                    sig_entry = sig_parser.SigParser(s.text.text).parse()
                    entry.addSignatureEntry(sig_entry)
                except sig_parser.SigParseException, e:
                    pass
                    #print >>sys.stderr, '\nWARNING: Could not parse signature: %s' % e
                    #print >>sys.stderr, 'Signature is: %s' % s.text.text.strip()
                
        return entry


class EnumConverter(CodeEntryConverter):
    def __init__(self, doc_proc):
        CodeEntryConverter.__init__(self, doc_proc)
        self.entry_class = ProcEnum
    
    def process(self, raw_entry):
        return CodeEntryConverter.process(self, raw_entry)


class AdaptionConverter(CodeEntryConverter):
    def __init__(self, doc_proc):
        CodeEntryConverter.__init__(self, doc_proc)
        self.entry_class = ProcAdaption
        self.parse_signature = False
    
    def process(self, raw_entry):
        return CodeEntryConverter.process(self, raw_entry)


class TypedefConverter(CodeEntryConverter):
    def __init__(self, doc_proc):
        CodeEntryConverter.__init__(self, doc_proc)
        self.entry_class = ProcTypedef
        self.parse_signature = False
    
    def process(self, raw_entry):
        return CodeEntryConverter.process(self, raw_entry)


class ConceptConverter(CodeEntryConverter):
    def __init__(self, doc_proc):
        CodeEntryConverter.__init__(self, doc_proc)
        self.entry_class = ProcConcept
    
    def process(self, raw_entry):
        concept = CodeEntryConverter.process(self, raw_entry)
        for e in raw_entry.extends:
            concept.addExtends(e.text.text.strip())
        return concept


class ClassConverter(CodeEntryConverter):
    def __init__(self, doc_proc):
        CodeEntryConverter.__init__(self, doc_proc)
        self.entry_class = ProcClass
    
    def process(self, raw_entry):
        klass = CodeEntryConverter.process(self, raw_entry)
        for e in raw_entry.extends:
            klass.addExtends(e.text.text.strip())
        for e in raw_entry.implements:
            klass.addImplements(e.text.text.strip())
        for t in raw_entry.tparams:
            proc_tparam = ProcTParam(t)
            proc_tparam.type = t.name.text
            proc_tparam.desc = self.rawTextToTextNode(t.text)
            klass.addTParam(proc_tparam)
        return klass


class TagConverter(CodeEntryConverter):
    def __init__(self, doc_proc):
        CodeEntryConverter.__init__(self, doc_proc)
        self.entry_class = ProcTag
        self.parse_signature = False

    def process(self, raw_entry):
        tag = CodeEntryConverter.process(self, raw_entry)
        for t in raw_entry.tparams:
            proc_tparam = ProcTParam()
            proc_tparam.type = t.name.text
            proc_tparam.desc = self.rawTextToTextNode(t.text)
            tag.addTParam(proc_tparam)
        return tag


class FunctionConverter(CodeEntryConverter):
    def __init__(self, doc_proc):
        CodeEntryConverter.__init__(self, doc_proc)
        self.entry_class = ProcFunction
        self.in_out_map = {
            'in': ProcParam.IN,
            'out': ProcParam.OUT,
            'in,out': ProcParam.IN_OUT,
            }
    
    def process(self, raw_entry):
        function = CodeEntryConverter.process(self, raw_entry)
        for p in raw_entry.params:
            proc_param = ProcParam(p)
            proc_param.name = p.name.text
            if p.inout:
                proc_param.in_out = self.in_out_map.get(p.inout.val[1:-1])
            proc_param.desc = self.rawTextToTextNode(p.text)
            function.addParam(proc_param)
        for t in raw_entry.tparams:
            proc_tparam = ProcTParam(t)
            proc_tparam.type = t.name.text
            proc_tparam.desc = self.rawTextToTextNode(t.text)
            function.addTParam(proc_tparam)
        for r in raw_entry.returns:
            proc_return = ProcReturn(r)
            proc_return.type = r.name.text
            proc_return.desc = self.rawTextToTextNode(r.text)
            function.addReturn(proc_return)
        for t in raw_entry.throws:
            proc_throw = ProcThrow(t)
            proc_throw.type = t.name.text
            proc_throw.desc = self.rawTextToTextNode(t.text)
            function.addThrow(proc_throw)
        for d in raw_entry.dataraces:
            proc_datarace = ProcDataRace(d)
            proc_datarace.desc = self.rawTextToTextNode(d.text) 
            function.addDataRace(proc_datarace)
        return function


class MacroConverter(CodeEntryConverter):
    def __init__(self, doc_proc):
        CodeEntryConverter.__init__(self, doc_proc)
        self.entry_class = ProcMacro
        self.in_out_map = {
            'in': ProcParam.IN,
            'out': ProcParam.OUT,
            'in,out': ProcParam.IN_OUT,
            }
        self.parse_signature = False
        
    def process(self, raw_entry):
        macro = CodeEntryConverter.process(self, raw_entry)
        for p in raw_entry.params:
            proc_param = ProcParam(p)
            proc_param.name = p.name.text
            if p.inout:
                proc_param.in_out = self.in_out_map.get(p.inout.val[1:-1])
            proc_param.desc = self.rawTextToTextNode(p.text)
            macro.addParam(proc_param)
        for r in raw_entry.returns:
            proc_return = ProcReturn(r)
            proc_return.type = r.name.text
            proc_return.desc = self.rawTextToTextNode(r.text)
            macro.addReturn(proc_return)
        for t in raw_entry.throws:
            proc_throw = ProcThrow(t)
            proc_throw.type = t.name.text
            proc_throw.desc = self.rawTextToTextNode(t.text)
            macro.addThrow(proc_throw)
        for d in raw_entry.dataraces:
            proc_datarace = ProcDataRace(d)
            proc_datarace.desc = self.rawTextToTextNode(d.text) 
            macro.addDataRace(proc_datarace)
        return macro


class MetafunctionConverter(CodeEntryConverter):
    def __init__(self, doc_proc):
        CodeEntryConverter.__init__(self, doc_proc)
        self.entry_class = ProcMetafunction

    def process(self, raw_entry):
        metafunction = CodeEntryConverter.process(self, raw_entry)
        for t in raw_entry.tparams:
            proc_tparam = ProcTParam(t)
            proc_tparam.type = t.name.text
            proc_tparam.desc = self.rawTextToTextNode(t.text)
            metafunction.addTParam(proc_tparam)
        for r in raw_entry.returns:
            proc_return = ProcReturn(r)
            proc_return.type = r.name.text
            proc_return.desc = self.rawTextToTextNode(r.text)
            metafunction.addReturn(proc_return)
        return metafunction


class VariableConverter(CodeEntryConverter):
    def __init__(self, doc_proc):
        CodeEntryConverter.__init__(self, doc_proc)
        self.entry_class = ProcVariable

    def process(self, raw_entry):
        variable = CodeEntryConverter.process(self, raw_entry)
        if raw_entry.type:
            variable.type = raw_entry.type.text
        return variable


class EnumValueConverter(CodeEntryConverter):
    def __init__(self, doc_proc):
        CodeEntryConverter.__init__(self, doc_proc)
        self.entry_class = ProcEnumValue
    
    def process(self, raw_entry):
        enum_value = CodeEntryConverter.process(self, raw_entry)
        if raw_entry.type:
            enum_value.type = raw_entry.type.text
        return enum_value


class TagStack(object):
    """Helper class for processing nested HTML tags."""

    def __init__(self):
        self.stack = []
         
    def push(self, token):
        pass
             
    def pop(self, token):
        pass


class PageConverter(EntryConverter):
    """Process a RawPage into a Page object."""

    def __init__(self, doc_proc):
        EntryConverter.__init__(self, doc_proc)
        self.entry_class = ProcPage


class GroupConverter(EntryConverter):
    """Process a RawGroup into a Group object."""

    def __init__(self, doc_proc):
        EntryConverter.__init__(self, doc_proc)
        self.entry_class = ProcGroup


class TextNodeVisitor(object):
    """Interface/abstract base class for visiting text nodes of entries or such."""

    def visit(self, text_node):
        """Visit TextNode, possibly translating its content.

        @param text_node: TextNode object or None.
        """
        pass


class LinkChecker(TextNodeVisitor):
    """Check raw link targets.

    Raw links are links of the form <a href="seqan:$target">$label</a>.
    """
    
    def __init__(self, doc):
        self.doc = doc

    def visit(self, text_node):
        if not text_node or text_node.type == '<text>':
            return
        if text_node.type == 'a':
            self._checkLink(text_node)
        else:
            for i, c in enumerate(text_node.children):
                self.visit(text_node.children[i])

    def _checkLink(self, a_node):
        if not a_node.attrs.get('href', '').startswith('seqan:'):
            return
        target = a_node.attrs['href'][6:]
        # TODO(holtgrew): Catch target_title being None, target_path not found!
        if target not in self.doc.entries:
            # TODO(holtgrew): Cannot resolve from TextNode to Token :(
            msg = 'Cannot find documentation entry "%s".' % target
            self.doc.doc_processor.msg_printer.printTokenError(
                a_node.tokens[0], msg, 'error')


class DocProcessor(object):
    """Convert a RawDoc object into a ProcDoc object.

    @ivar converters: Dict that maps RawEntry kinds to Converter objects.
    @ivar include_dirs: The bases path for including files that can be used in
                        the @include and @snippet commands.
    @ivar include_mgr: inc_mgr.IncludeManager object for file/snippet
                       inclusion.
    @ivar expected_tags: Iterateable of expected tag names.  Will warn in
                         conversion about unexpected tags if hit.
    @ivar msg_printer: The dox_parser.MessagePrinter instance to use for
                       printing messages.
    """
    
    def __init__(self, logger=None, include_dirs=['.'], expected_tags=[],
                 msg_printer=None):
        self.logger = logger
        self.include_dirs = list(include_dirs)
        self.include_mgr = inc_mgr.IncludeManager(self.include_dirs)
        self.expected_tags = set(expected_tags)
        self.converters = {
            'class': ClassConverter(self),
            'concept': ConceptConverter(self),
            'enum': EnumConverter(self),
            'adaption': AdaptionConverter(self),
            'global_typedef': TypedefConverter(self),
            'member_typedef': TypedefConverter(self),
            'grouped_typedef': TypedefConverter(self),
            'global_function': FunctionConverter(self),
            'global_metafunction': MetafunctionConverter(self),
            'defgroup': GroupConverter(self),
            'grouped_macro' : MacroConverter(self),
            'interface_function': FunctionConverter(self),
            'interface_metafunction': MetafunctionConverter(self),
            'macro' : MacroConverter(self),
            'member_function': FunctionConverter(self),
            'member_variable': VariableConverter(self),
            'enum_value': EnumValueConverter(self),
            'page': PageConverter(self),
            'tag': TagConverter(self),
            'grouped_tag': TagConverter(self),
            'variable': VariableConverter(self),
            }
        self.msg_printer = msg_printer or dox_parser.MessagePrinter()
        self.validators = [x(self.msg_printer) for x in validation.VALIDATORS]
        self.entry_filenames = []
        self.topLevelEntry_filenames = {} 
        self.secondLevelEntry_filenames = {} 

    def run(self, doc):
        res = ProcDoc(self)
        self.entry_filenames = doc.filenames
        self.log('Processing Documentation...')
        self.convertTopLevelEntries(doc, res)
        self.convertSecondLevelEntries(doc, res)
        self.convertVariables(doc, res)
        self.checkLinks(doc, res)
        self.buildInheritanceLists(res)
        self.validate(res)
        return res

    def convertTopLevelEntries(self, doc, res):
        """Convert top level entries.

        Variables are not converted yet.  They are converted in a separate
        step since they might encode enum values.
        """
        self.log('  1) Converting Top-Level Entries.')
        #print 'doc.entries', [e.name.text for e in doc.entries]
        for index,raw_entry in enumerate(doc.entries):
            # Get fitting converter or warn if there is none.
            kind = raw_entry.getType()
            if not kind in ['concept', 'class', 'global_function',
                            'global_metafunction', 'page', 'tag',
                            'defgroup', 'macro', 'adaption', 'global_typedef', 'enum']:
                continue  # Not a top-level entry.
            converter = self.converters.get(kind)
            if not converter:
                self.logWarning('Could not find converter for kind "%s".', kind)
                continue  # Skip if no converter could be found.
            # Perform conversion.
            proc_entry = converter.process(raw_entry)
            # Store object in ProcDoc.
            #self.log('    * %s (%s)' % (proc_entry.name, proc_entry))
            res.addTopLevelEntry(proc_entry)
            self.topLevelEntry_filenames[proc_entry] = self.entry_filenames[index]

    def convertSecondLevelEntries(self, doc, res):
        self.log('  2) Converting Second-Level Entries.')
        for index,raw_entry in enumerate(doc.entries):
            # Get fitting converter or warn if there is none.
            kind = raw_entry.getType()
            if not kind in ['member_function', 'interface_function',
                            'interface_metafunction',
                            'grouped_tag', 'grouped_macro', 'member_typedef',
                            'grouped_typedef']:
                continue  # Not a top-level entry.
            converter = self.converters.get(kind)
            if not converter:
                self.logWarning('Could not find converter for kind "%s".', kind)
                continue  # Skip if no converter could be found.
            # Perform conversion.
            proc_entry = converter.process(raw_entry)
            # Store object in ProcDoc.
            #self.log('    * %s' % proc_entry.name)
            res.addSecondLevelEntry(proc_entry)
            self.secondLevelEntry_filenames[proc_entry] = self.entry_filenames[index]

    def convertVariables(self, doc, res):
        self.log('  3) Converting variable and enum value entries.')
        var_types = ['member_variable', 'grouped_variable', 'variable']
        for raw_entry in [e for e in doc.entries if e.getType() in var_types]:
            kind = raw_entry.getType()
            converter = self.converters.get(kind)
            if not converter:
                self.logWarning('Could not find converter for kind "%s".', kind)
                continue  # Skip if no converter could be found.
            # Perform conversion.
            proc_entry = converter.process(raw_entry)
            # Store object in ProcDoc.
            #self.log('    * %s %s' % (proc_entry.type, proc_entry.name))
            res.addVariable(proc_entry)
        for raw_entry in [e for e in doc.entries if e.getType() == 'enum_value']:
            converter = self.converters.get(raw_entry.getType())
            if not converter:
                self.logWarning('Could not find converter for kind "%s".', kind)
                continue  # Skip if no converter could be found.
            # Perform conversion.
            proc_entry = converter.process(raw_entry)
            # Store object in ProcDoc.
            #self.log('    * %s %s' % (proc_entry.type, proc_entry.name))
            res.addEnumValue(proc_entry)

    def checkLinks(self, doc, res):
        """Check <link> items of text nodes and references.

        References are given either explicitely in items like @extends and
        @implements.
        """
        self.log('  3) Checking References.')
        link_checker = LinkChecker(res)
        for proc_entry in res.entries.values():
            proc_entry.visitTextNodes(link_checker)

    def buildInheritanceLists(self, doc):
        """Build lists regarding the inheritance in the classes and concepts in doc.

        We will build the equivalent to what Javadoc builds.

        For concepts, this is the list of (a) all extended concepts, (b) all
        known extending concepts, (c) all known implementing classes.

        For classes, this is the list of (a) all implemented concepts, (b) all
        direct known specializations, (c) all extended classes.

        @param doc: The ProcDoc object with the classes and concept.

        """
        self.log('  4) Building Inheritance Lists.')
        # Process concepts: All extended and all extending.
        concepts = [x for x in doc.top_level_entries.values()
                    if x.kind == 'concept']
        # Get all concepts that c extends into c.all_extended.
        for c in concepts:
            q = list(c.extends)  # Queue for recursion
            while q:
                name = q[0]
                q.pop(0)
                if name in c.all_extended:
                    continue  # Skip to break loops.
                c.all_extended.add(name)
                q += doc.top_level_entries[name].extends
        # Now, build list of all extending concepts into c.all_extending.
        for c in concepts:
            for name in c.all_extended:
                doc.top_level_entries[name].all_extending.add(c.name)
        # Process classes: All extended and all extending classes.
        classes = [x for x in doc.top_level_entries.values()
                   if x.kind in ['class', 'specialization']]
        # Get all classes that c extends into c.all_extended.
        for c in classes:
            q = list(c.extends)  # Queue for recursion
            while q:
                name = q[0]
                q.pop(0)
                if name in c.all_extended:
                    continue  # Skip to break loops.
                c.all_extended.add(name)
                if name not in doc.top_level_entries:
                    self.logWarning('Could not find entry for extending: %s', name)
                    continue
                q += doc.top_level_entries[name].extends
        # Now, build list of all extending clsses into c.all_extending.
        for c in classes:
            for name in c.all_extended:
                if name not in doc.top_level_entries:
                    self.logWarning('Could not find entry for extending: %s', name)
                    continue
                doc.top_level_entries[name].all_extending.add(c.name)
        # Build list of all direct implementing classes for all concepts.
        for cl in classes:
            for name in cl.implements:
                if '\u0001' in name:
                    continue  # Skip transitive inheritance.
                if not doc.top_level_entries.get(name):
                    self.logWarning('Could not find entry for implementing: %s', name)
                    continue
                co = doc.top_level_entries[name]
                if co.kind != 'concept':
                    self.logWarning('Only concepts can be implemented.')
                    continue
                co.all_implementing.add(cl.name)
                co.all_implementing.update(cl.all_extending)
        # Build list of all implemented concepts for all classes.
        for co in concepts:
            for name in co.all_implementing:
                cl = doc.top_level_entries[name]
                cl.all_implemented.add(co.name)
                cl.all_implemented.update(co.all_extended)  # inheritance
        # Update list of all implementing classes for all concepts (transitive)
        for cl in classes:
            for name in cl.all_implemented:
                co = doc.top_level_entries[name]
                co.all_implementing.add(cl.name)

    def validate(self, doc):
        """Execute validation using the validators from self.validators.
        
        @param doc: The ProcDoc object to validate.
        """
        self.log('  5) Running validation.')
        for name, entry in doc.entries.iteritems():
            for v in self.validators:
                #print v, entry
                v.validate(entry)

    def log(self, msg, *args, **kwargs):
        """Print the given message to the configured logger if any.
        """
        if not self.logger:
            return
        self.logger.info(msg, *args, **kwargs)
        
    def logWarning(self, msg, *args, **kwargs):
        """Print the given message to the configured logger if any.
        """
        if not self.logger:
            return
        self.logger.warning('WARNING: ' + msg, *args, **kwargs)
