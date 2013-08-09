#!/usr/bin/env python
"""SeqAn documentation raw object representation.

This is the direct representation as it can be determined from the embedded
Doxygen-style comments without the interpretation of commands within clauses
and cross-linking.
"""

import textwrap
import dox_tokens


class DoxFormatter(object):
    """Formatter for printing correctly indented and wrapped in Doxygen style.
    """
    
    def __init__(self, width=77):
        self.width = width

    def formatCommand(self, name, text, leading=None):
        """RawReturn string with a formatted command.
        
        The general format is "@$name $leading $text" where the text is wrapped
        to the end of leading.
        """
        if leading:
            res = ['@', name, ' ', leading,  ' ']
        else:
            res = ['@', name, ' ']
        l = len(''.join(res))
        indent = ' ' * l
        wrapped_text = textwrap.wrap(text, self.width - l)
        if wrapped_text:
            res.append(wrapped_text[0])
        for x in wrapped_text[1:]:
            res += ['\n', indent, x]
        return ''.join(res) + '\n'

    def formatParagraph(self, text):
        """Format paragraph."""
        return '\n'.join(textwrap.wrap(text, self.width)) + '\n'


class RawText(object):
    """List of token with easy concatenation into a string.
    
    This type is used for collecting lists of tokens.
    
    @ivar tokens: The list of token objects.
    """
    
    def __init__(self, tokens=[]):
        self.tokens = list(tokens)

    def append(self, token):
        """Append the token to the list of tokens.

        @param token: The lexer.Token object to add.
        @return: Nothing
        """
        self.tokens.append(token)

    @property
    def empty(self):
        """RawReturns whether the token set is empty.

        @return: Whether or not the token list is empty.
        """
        return not bool(self.tokens)

    @property
    def text(self):
        """RawReturns the concatenated tokens' text.

        @return: The concatenated tokens' text.
        """
        return ''.join([x.val for x in self.tokens])

    def __eq__(self, other):
        if not hasattr(other, 'tokens'):
            return False
        return self.tokens == other.tokens


class RawDoc(object):
    """The documentation consists of a number of documentation objects.
    
    @ivar entries List of RawEntry objects.
    """

    def __init__(self):
        self.entries = []

    def merge(self, other_doc):
        for e in other_doc.entries:
            self.addEntry(e)

    def addEntry(self, entry):
        self.entries.append(entry)

    def getFormatted(self, width=77):
        """Get formatted and normalized in dox format."""
        formatter = DoxFormatter(width)
        res = []
        first = True
        for entry in self.entries:
            res.append(entry.getFormatted(formatter))
            first = False
        return '\n\n'.join(res)


class RawEntry(object):
    """One top-level entry of the documentation.
   
    @ivar name  The identifier of the entry.
    @ivar title The title of the entry.
    @ivar brief A string object with a brief summary of the entry.
    @ivar body  A RawBody object with the entry's documentation.
    @ivar sees  A list of RawSee objects.
    @ivar command The name of the command starting the entry type.
    """
    
    def __init__(self, briefs=[], command='<entry>'):
        self.name = RawText()
        self.title = RawText()
        self.briefs = list(briefs)
        self.body = RawBody()
        self.sees = []
        self.command = command

    def addBrief(self, b):
        self.briefs.append(b)

    def addSee(self, see):
        while see.text.tokens and see.text.tokens[-1].type in dox_tokens.WHITESPACE:
            see.text.tokens.pop()
        self.sees.append(see)

    @classmethod
    def entryTypes(cls):
        """RawReturns iterable with all entry types."""
        res = ('concept', 'class', 'function', 'metafunction', 'page', 'enum', 'var',
               'tag', 'defgroup', 'macro')
        return res

    def addParagraph(self, p):
        self.body.addParagraph(p)

    def getFormatted(self, formatter):
        """Get formatted and normalized in dox format."""
        res = []
        if self.title.text:
            res.append(formatter.formatCommand(self.command, self.title.text,
                                               self.name.text))
        else:
            res.append(formatter.formatCommand(self.command, self.name.text))
        if self.briefs:
            res.append('\n')
        for x in self.briefs:
            res.append(x.getFormatted(formatter))
        if not self.body.empty:
            res += ['\n', self.body.getFormatted(formatter)]
        if self.sees:
            res.append('\n')
        for x in self.sees:
            res.append(x.getFormatted(formatter))
        res.append('\n')
        return ''.join(res)


class RawCodeEntry(RawEntry):
    """RawDoc for one code entry concept having a signature.
    
    @ivar signatures A list of RawSignature objects.
    """
    
    def __init__(self, briefs=[], command='<code entry>'):
        RawEntry.__init__(self, briefs=briefs, command=command)
        self.signatures = []
        self.headerfiles = []
        self.deprecation_msgs = []

    def addSignature(self, s):
        self.signatures.append(s)

    def addHeaderfile(self, h):
        self.headerfiles.append(h)

    def addDeprecationMsg(self, h):
        self.deprecation_msgs.append(h)

    def getType(self):
        return 'code'

    def __str__(self):
        res = RawEntry.__str__(self)
        return res + '\n' + '\n'.join(['  @signature %s' % x for x in self.signatures])

    def getFormatted(self, formatter):
        res = []
        if self.title.text:
            res.append(formatter.formatCommand(self.command, self.title.text,
                                               self.name.text))
        else:
            res.append(formatter.formatCommand(self.command, self.name.text))
        if self.headerfiles:
            res.append('\n')
        if self.headerfiles:
            for x in self.headerfiles:
                res.append(x.getFormatted(formatter))
        if self.briefs:
            res.append('\n')
        for x in self.briefs:
            res.append(x.getFormatted(formatter))
        if self.deprecation_msgs:
            res.append('\n')
        for x in self.deprecation_msgs:
            res.append(x.getFormatted(formatter))
        if self.signatures:
            res.append('\n')
        for x in self.signatures:
            res.append(x.getFormatted(formatter))
        if not self.body.empty:
            res.append('\n')
        res += self.body.getFormatted(formatter)
        if self.sees:
            res.append('\n')
        for x in self.sees:
            res.append(x.getFormatted(formatter))
        res.append('\n')
        return ''.join(res)


class RawVariable(RawCodeEntry):
    """RawDoc for one variable or enum constant.

    @ivar type: The type of the variable as a RawText or None.
    """

    def __init__(self, briefs=[]):
        RawCodeEntry.__init__(self, briefs=briefs, command='var')
        self.type = None

    def getType(self):
        if '::' in self.name.text:
            return 'member_variable'
        else:
            return 'variable'

    def getFormatted(self, formatter):
        res = []
        if self.type:
            res.append(formatter.formatCommand(self.command, self.name.text, self.type.text))
        else:
            res.append(formatter.formatCommand(self.command, self.name.text))
        if self.headerfiles:
            res.append('\n')
        if self.headerfiles:
            for x in self.headerfiles:
                res.append(x.getFormatted(formatter))
        if self.briefs:
            res.append('\n')
        for x in self.briefs:
            res.append(x.getFormatted(formatter))
        if self.deprecation_msgs:
            res.append('\n')
        for x in self.deprecation_msgs:
            res.append(x.getFormatted(formatter))
        if self.signatures:
            res.append('\n')
        for x in self.signatures:
            res.append(x.getFormatted(formatter))
        if not self.body.empty:
            res.append('\n')
        res += self.body.getFormatted(formatter)
        if self.sees:
            res.append('\n')
        for x in self.sees:
            res.append(x.getFormatted(formatter))
        res.append('\n')
        return ''.join(res)


class RawTag(RawCodeEntry):
    """RawDoc for one tag."""

    def __init__(self, briefs=[]):
        RawCodeEntry.__init__(self, briefs=briefs, command='tag')
        self.tparams = []

    def addTParam(self, p):
        self.tparams.append(p)

    def getType(self):
        if '#' in self.name.text:
            return 'grouped_tag'
        else:
            return 'tag'
        
    def getFormatted(self, formatter):
        res = []
        if self.title.text:
            res.append(formatter.formatCommand(self.command, self.title.text,
                                               self.name.text))
        else:
            res.append(formatter.formatCommand(self.command, self.name.text))
        if self.headerfiles:
            res.append('\n')
        if self.headerfiles:
            for x in self.headerfiles:
                res.append(x.getFormatted(formatter))
        if self.briefs:
            res.append('\n')
        for x in self.briefs:
            res.append(x.getFormatted(formatter))
        if self.deprecation_msgs:
            res.append('\n')
        for x in self.deprecation_msgs:
            res.append(x.getFormatted(formatter))
        if self.signatures:
            res.append('\n')
        for x in self.signatures:
            res.append(x.getFormatted(formatter))
        for x in self.tparams:
            res.append(x.getFormatted(formatter))
        if not self.body.empty:
            res.append('\n')
        res += self.body.getFormatted(formatter)
        if self.sees:
            res.append('\n')
        for x in self.sees:
            res.append(x.getFormatted(formatter))
        res.append('\n')
        return ''.join(res)


class RawConcept(RawCodeEntry):
    """RawDoc for one concept.
    """

    def __init__(self, briefs=[]):
        RawCodeEntry.__init__(self, briefs=briefs, command='concept')
        self.extends = []

    def addExtends(self, c):
        self.extends.append(c)

    def getType(self):
        return 'concept'

    def getFormatted(self, formatter):
        res = []
        if self.title.text:
            res.append(formatter.formatCommand(self.command, self.title.text,
                                               self.name.text))
        else:
            res.append(formatter.formatCommand(self.command, self.name.text))
        if self.headerfiles:
            res.append('\n')
        if self.headerfiles:
            for x in self.headerfiles:
                res.append(x.getFormatted(formatter))
        if self.extends:
            res.append('\n')
        for x in self.extends:
            res.append(x.getFormatted(formatter))
        if self.briefs:
            res.append('\n')
        for x in self.briefs:
            res.append(x.getFormatted(formatter))
        if self.deprecation_msgs:
            res.append('\n')
        for x in self.deprecation_msgs:
            res.append(x.getFormatted(formatter))
        if self.signatures:
            res.append('\n')
        for x in self.signatures:
            res.append(x.getFormatted(formatter))
        if not self.body.empty:
            res.append('\n')
        res += self.body.getFormatted(formatter)
        if self.sees:
            res.append('\n')
        for x in self.sees:
            res.append(x.getFormatted(formatter))
        res.append('\n')
        return ''.join(res)


class RawEnum(RawCodeEntry):
    """RawDoc for one enum."""

    def __init__(self, briefs=[]):
        RawCodeEntry.__init__(self, briefs=briefs, command='enum')

    def getType(self):
        return 'enum'


class RawTypedef(RawCodeEntry):
    """RawDoc for one typedef."""

    def __init__(self, briefs=[]):
        RawCodeEntry.__init__(self, briefs=briefs, command='typedef')

    def getType(self):
        if '#' in self.name.text:
            return 'grouped_typedef'
        elif '::' in self.name.text:
            return 'member_typedef'
        else:
            return 'global_typedef'


class RawAdaption(RawCodeEntry):
    """RawDoc for one adaption."""

    def __init__(self, briefs=[]):
        RawCodeEntry.__init__(self, briefs=briefs, command='adaption')

    def getType(self):
        return 'adaption'


class RawClass(RawCodeEntry):
    """RawDoc for one class.
    
    @ivar tparams List of RawParameter objects.
    """

    def __init__(self, briefs=[]):
        RawCodeEntry.__init__(self, briefs=briefs, command='class')
        self.extends = []
        self.implements = []
        self.tparams = []

    def addTParam(self, p):
        self.tparams.append(p)

    def addExtends(self, p):
        self.extends.append(p)

    def addImplements(self, p):
        self.implements.append(p)

    def getType(self):
        return 'class'

    def getFormatted(self, formatter):
        res = []
        if self.title.text:
            res.append(formatter.formatCommand(self.command, self.title.text,
                                               self.name.text))
        else:
            res.append(formatter.formatCommand(self.command, self.name.text))
        if self.implements:
            res.append('\n')
        for x in self.implements:
            res.append(x.getFormatted(formatter))
        if self.extends:
            res.append('\n')
        for x in self.extends:
            res.append(x.getFormatted(formatter))
        if self.headerfiles:
            res.append('\n')
        if self.headerfiles:
            for x in self.headerfiles:
                res.append(x.getFormatted(formatter))
        if self.briefs:
            res.append('\n')
        for x in self.briefs:
            res.append(x.getFormatted(formatter))
        if self.deprecation_msgs:
            res.append('\n')
        for x in self.deprecation_msgs:
            res.append(x.getFormatted(formatter))
        if self.signatures:
            res.append('\n')
        for x in self.signatures:
            res.append(x.getFormatted(formatter))
        if self.tparams:
            res.append('\n')
        for x in self.tparams:
            res.append(x.getFormatted(formatter))
        if not self.body.empty:
            res.append('\n')
        res += self.body.getFormatted(formatter)
        if self.sees:
            res.append('\n')
        for x in self.sees:
            res.append(x.getFormatted(formatter))
        res.append('\n')
        return ''.join(res)

    def __str__(self):
        res = RawCodeEntry.__str__(self)
        s = res + '\n'
        if self.tparams:
            s += '\n'.join(['  @tparam %s' % s for s in self.tparams]) + '\n'
        return s


class RawFunction(RawCodeEntry):
    """RawDoc for one function.
    
    @ivar tparams List of RawParameter objects.
    @ivar params List of RawParameter objects.
    @ivar returns List of RawReturn objects.
    """
    
    def __init__(self, briefs=[]):
        RawCodeEntry.__init__(self, briefs=briefs, command='fn')
        self.tparams = []
        self.params = []
        self.returns = []

    def addTParam(self, p):
        self.tparams.append(p)

    def addParam(self, p):
        self.params.append(p)

    def addReturn(self, p):
        self.returns.append(p)

    def getType(self):
        if '#' in self.name.text:
            return 'interface_function'
        elif '::' in self.name.text:
            return 'member_function'
        else:
            return 'global_function'

    def getFormatted(self, formatter):
        res = []
        res.append(formatter.formatCommand(self.command, self.name.text))
        if self.headerfiles:
            res.append('\n')
        if self.headerfiles:
            for x in self.headerfiles:
                res.append(x.getFormatted(formatter))
        if self.briefs:
            res.append('\n')
        for x in self.briefs:
            res.append(x.getFormatted(formatter))
        if self.deprecation_msgs:
            res.append('\n')
        for x in self.deprecation_msgs:
            res.append(x.getFormatted(formatter))
        if self.signatures:
            res.append('\n')
        for x in self.signatures:
            res.append(x.getFormatted(formatter))
        if self.tparams:
            res.append('\n')
        for x in self.tparams:
            res.append(x.getFormatted(formatter))
        if self.params:
            res.append('\n')
        for x in self.params:
            res.append(x.getFormatted(formatter))
        if self.returns:
            res.append('\n')
        for x in self.returns:
            res.append(x.getFormatted(formatter))
        if not self.body.empty:
            res.append('\n')
        res += self.body.getFormatted(formatter)
        if self.sees:
            res.append('\n')
        for x in self.sees:
            res.append(x.getFormatted(formatter))
        res.append('\n')
        return ''.join(res)

    def __str__(self):
        res = RawCodeEntry.__str__(self)
        res += '\n' + '\n'.join(['  @return %s ' % x for x in self.returns])
        res += '\n' + '\n'.join(['  @tparam %s ' % x for x in self.tparams])
        res += '\n' + '\n'.join(['  @param %s ' % x for x in self.params])
        res += '\n'
        return res


class RawMacro(RawCodeEntry):
    """RawDoc for one function.
    
    @ivar params List of RawParameter objects.
    @ivar returns List of RawReturn objects.
    """
    
    def __init__(self, briefs=[]):
        RawCodeEntry.__init__(self, briefs=briefs, command='macro')
        self.params = []
        self.returns = []

    def addParam(self, p):
        self.params.append(p)

    def addReturn(self, p):
        self.returns.append(p)

    def getType(self):
        if '#' in self.name.text:
            return 'grouped_macro'
        else:
            return 'macro'

    def getFormatted(self, formatter):
        res = []
        res.append(formatter.formatCommand(self.command, self.name.text))
        if self.headerfiles:
            res.append('\n')
        if self.headerfiles:
            for x in self.headerfiles:
                res.append(x.getFormatted(formatter))
        if self.briefs:
            res.append('\n')
        for x in self.briefs:
            res.append(x.getFormatted(formatter))
        if self.deprecation_msgs:
            res.append('\n')
        for x in self.deprecation_msgs:
            res.append(x.getFormatted(formatter))
        if self.signatures:
            res.append('\n')
        for x in self.signatures:
            res.append(x.getFormatted(formatter))
        if self.params:
            res.append('\n')
        for x in self.params:
            res.append(x.getFormatted(formatter))
        if self.returns:
            res.append('\n')
        for x in self.returns:
            res.append(x.getFormatted(formatter))
        if not self.body.empty:
            res.append('\n')
        res += self.body.getFormatted(formatter)
        if self.sees:
            res.append('\n')
        for x in self.sees:
            res.append(x.getFormatted(formatter))
        res.append('\n')
        return ''.join(res)

    def __str__(self):
        res = RawCodeEntry.__str__(self)
        res += '\n' + '\n'.join(['  @return %s ' % x for x in self.returns])
        res += '\n' + '\n'.join(['  @tparam %s ' % x for x in self.tparams])
        res += '\n' + '\n'.join(['  @param %s ' % x for x in self.params])
        res += '\n'
        return res


class RawMetafunction(RawCodeEntry):
    """RawDoc for one metafunction.
    
    @ivar tparams List of RawParameter objects.
    @ivar returns List of RawReturn objects.
    """
    
    def __init__(self, briefs=[]):
        RawCodeEntry.__init__(self, briefs=briefs)
        self.tparams = []
        self.returns = []
        self.command = 'mfn'

    def addTParam(self, p):
        self.tparams.append(p)

    def addReturn(self, p):
        self.returns.append(p)

    def getType(self):
        if '#' in self.name.text:
            return 'interface_metafunction'
        else:
            return 'global_metafunction'

    def getFormatted(self, formatter):
        res = []
        res.append(formatter.formatCommand(self.command, self.name.text))
        if self.headerfiles:
            res.append('\n')
        if self.headerfiles:
            for x in self.headerfiles:
                res.append(x.getFormatted(formatter))
        if self.briefs:
            res.append('\n')
        for x in self.briefs:
            res.append(x.getFormatted(formatter))
        if self.deprecation_msgs:
            res.append('\n')
        for x in self.deprecation_msgs:
            res.append(x.getFormatted(formatter))
        if self.signatures:
            res.append('\n')
        for x in self.signatures:
            res.append(x.getFormatted(formatter))
        if self.tparams:
            res.append('\n')
        for x in self.tparams:
            res.append(x.getFormatted(formatter))
        if self.returns:
            res.append('\n')
        for x in self.returns:
            res.append(x.getFormatted(formatter))
        if not self.body.empty:
            res.append('\n')
        res += self.body.getFormatted(formatter)
        if self.sees:
            res.append('\n')
        for x in self.sees:
            res.append(x.getFormatted(formatter))
        res.append('\n')
        return ''.join(res)


class RawPage(RawEntry):
    """A page in the documentation."""

    def __init__(self, briefs=[]):
        RawEntry.__init__(self, briefs=briefs)
        self.command = 'page'
    
    def getType(self):
        return 'page'

    def getFormatted(self, formatter):
        res = []
        if self.title.text:
            res.append(formatter.formatCommand(self.command, self.title.text,
                                               self.name.text))
        else:
            res.append(formatter.formatCommand(self.command, self.name.text))
        if self.briefs:
            res.append('\n')
        for x in self.briefs:
            res.append(x.getFormatted(formatter))
        if not self.body.empty:
            res.append('\n')
        res += self.body.getFormatted(formatter)
        if self.sees:
            res.append('\n')
        for x in self.sees:
            res.append(x.getFormatted(formatter))
        res.append('\n')
        return ''.join(res)


class RawMainPage(RawPage):
    """The main page in the documentation."""
    
    def __init__(self, briefs=[]):
        RawPage.__init__(self, briefs=briefs)
        self.command = 'mainpage'
    
    def getType(self):
        return 'page'

    def getFormatted(self, formatter):
        res = []
        if self.title.text:
            res.append(formatter.formatCommand(self.command, self.title.text))
        else:
            res.append(formatter.formatCommand(self.command, 'NO TITLE'))
        if self.briefs:
            res.append('\n')
        for x in self.briefs:
            res.append(x.getFormatted(formatter))
        if not self.body.empty:
            res.append('\n')
        res += self.body.getFormatted(formatter)
        if self.sees:
            res.append('\n')
        for x in self.sees:
            res.append(x.getFormatted(formatter))
        res.append('\n')
        return ''.join(res)


class RawGroup(RawEntry):
    """A group in the documentation."""

    def __init__(self, briefs=[]):
        RawEntry.__init__(self, briefs=briefs)
        self.command = 'defgroup'
    
    def getType(self):
        return 'defgroup'

    def getFormatted(self, formatter):
        res = []
        if self.title.text:
            res.append(formatter.formatCommand(self.command, self.title.text,
                                               self.name.text))
        else:
            res.append(formatter.formatCommand(self.command, self.name.text))
        if self.briefs:
            res.append('\n')
        for x in self.briefs:
            res.append(x.getFormatted(formatter))
        if not self.body.empty:
            res.append('\n')
        res += self.body.getFormatted(formatter)
        if self.sees:
            res.append('\n')
        for x in self.sees:
            res.append(x.getFormatted(formatter))
        res.append('\n')
        return ''.join(res)


class RawBody(object):
    """A documentation body consists of multiple RawParagraph, RawSection, RawInclude objects.
    
    @ivar entries A list of RawParagraph and RawSection objects.
    """
    
    def __init__(self):
        self.paragraphs = []

    def addParagraph(self, p):
        self.paragraphs.append(p)

    def getFormatted(self, formatter):
        res = []
        for p in self.paragraphs:
            res.append(p.getFormatted(formatter))
        return '\n'.join(res)

    @property
    def empty(self):
        return not bool(self.paragraphs)

    def __eq__(self, other):
        return self.paragraphs == other.paragraphs


class RawSection(object):
    """Represents one section or subsection.
    
    @ivar level An int with the indentation level, starts at 0.
    @ivar heading The text of the heading.
    """

    def getType(self):
        return 'section'

    def __init__(self, heading=RawText(), level=0):
        self.heading = heading
        self.level = level

    def __str__(self):
        if self.level == 0:
            return 'Section(%s)' % (repr(self.heading.text))
        else:
            return 'Sub%ssection(%s)' % (''.join(['sub'] * (self.level - 1)), repr(self.heading.text))

    def getCommand(self):
        if self.level == 0:
            return 'section'
        else:
            return 'sub%ssection' % ''.join(['sub'] * (self.level - 1))

    def getFormatted(self, formatter):
        res = [formatter.formatCommand(self.getCommand(), self.heading.text.strip())]
        return ''.join(res)


class RawInclude(object):
    """An @include statement.

    @ivar path: A RawText object with the path to the included file.
    """

    def __init__(self, path):
        self.path = path

    def getType(self):
        return 'include'

    def __str__(self):
        return 'RawInclude(%s)' % (repr(self.path.text),)

    def getFormatted(self, formatter):
        res = ['@include ', self.path.text.strip(), '\n']
        return ''.join(res)


class RawSnippet(object):
    """A @snippet statement.

    @ivar path: A RawText object with the path to the included file.
    @ivar name: The name of the snippet, a RawText.
    """

    def __init__(self, path, name):
        self.path = path
        self.name = name

    def getType(self):
        return 'snippet'

    def __str__(self):
        return 'RawSnippet(%s, %s)' % (repr(self.path.text),
                                       repr(self.name.text))

    def getFormatted(self, formatter):
        res = ['@snippet ', self.path.text.strip(), ' ',  self.name.text.strip(), '\n']
        return ''.join(res)


class RawParagraph(object):
    """A paragraph in the RawBody of an RawEntry object's documentation.
    
    @ivar text A string with the paragraph's text.
    """

    def __init__(self, text=RawText()):
        self.text = text

    def getType(self):
        return 'paragraph'

    def __str__(self):
        return 'RawParagraph(%s)' % (repr(self.text.text))

    def getFormatted(self, formatter):
        return formatter.formatParagraph(self.text.text)


class RawCode(RawParagraph):
    """A special paragraph that is rendered as code.

    @ivar extension The extension identifying the language.
    """
    
    def __init__(self, text=RawText(), extension='.txt'):
        RawParagraph.__init__(self, text)
        self.extension = extension

    def getType(self):
        return 'code'

    def __str__(self):
        return 'RawCode(%s)' % repr(self.text)

    def getFormatted(self, formatter):
        return '@code%s@endcode' % self.text.text


class RawBrief(object):
    """A representation of a @brief entry.
    
    @ivar text The @brief clauses's text.
    """
    
    def __init__(self, text):
        self.text = text

    def getType(self):
        return 'brief'

    def getFormatted(self, formatter):
        return formatter.formatCommand('brief', self.text.text.strip())

    def __str__(self):
        return 'RawBrief(%s)' % repr(self.text)

    def __eq__(self, other):
        return self.text == other.text


class RawExtends(object):
    """A representation of a @extends entry.
    
    @ivar text The @extends clauses's text.
    """
    
    def __init__(self, text):
        self.text = text

    def getType(self):
        return 'extends'

    def getFormatted(self, formatter):
        return formatter.formatCommand('extends', self.text.text.strip())

    def __str__(self):
        return 'RawExtends(%s)' % repr(self.text)

    def __eq__(self, other):
        return self.text == other.text


class RawImplements(object):
    """A representation of a @implements entry.
    
    @ivar text The @implements clauses's text.
    """
    
    def __init__(self, text):
        self.text = text

    def getType(self):
        return 'implements'

    def getFormatted(self, formatter):
        return formatter.formatCommand('implements', self.text.text.strip())

    def __str__(self):
        return 'RawImplements(%s)' % repr(self.text)

    def __eq__(self, other):
        return self.text == other.text


class RawHeaderfile(object):
    """A representation of a @headerfile entry.
    
    @ivar text The @headerfile clauses's text.
    """
    
    def __init__(self, text):
        self.text = text

    def getType(self):
        return 'headerfile'

    def getFormatted(self, formatter):
        return formatter.formatCommand('headerfile', self.text.text.strip())

    def __str__(self):
        return 'RawHeaderfile(%s)' % repr(self.text)

    def __eq__(self, other):
        return self.text == other.text


class RawDeprecated(object):
    """A representation of a @deprecated entry.
    
    @ivar text The @deprecated clauses's text.
    """
    
    def __init__(self, text):
        self.text = text

    def getType(self):
        return 'deprecated'

    def getFormatted(self, formatter):
        return formatter.formatCommand('deprecated', self.text.text.strip())

    def __str__(self):
        return 'RawDeprecated(%s)' % repr(self.text)

    def __eq__(self, other):
        return self.text == other.text


class RawSee(object):
    """A representation of a @see entry.
    
    @ivar value The @see clauses's parameter.
    """
    
    def __init__(self, text):
        self.text = text

    def getType(self):
        return 'see'

    def getFormatted(self, formatter):
        return formatter.formatCommand('see', self.text.text)


class RawParam(object):
    """RawDoc for one parameter.
    
    @ivar name  Name of the parameter.
    @ivar inout String in {'', 'in', 'out', 'in,out'} describing mutability.
    @ivar text  RawParagraph entry with the documentation of the parameter.
    """

    def getType(self):
        return 'param'
    
    def __init__(self, name=RawText(), text=RawText(), inout=None):
        self.name = name
        self.inout = inout
        self.text = text

    def __str__(self):
        inout = ''
        if self.inout:
            inout = self.inout.val
        return '  @param%s %s %s' % (inout, self.name.text, self.text.text)

    def getFormatted(self, formatter):
        inout = ''
        if self.inout:
            inout = self.inout.val
        return formatter.formatCommand('param%s' % inout, self.text.text, self.name.text)


class RawTParam(RawParam):
    """RawDoc for one template parameter.
    
    @ivar name  Name of the parameter.
    @ivar text  RawParagraph entry with the documentation of the parameter.
    """

    def getType(self):
        return 'tparam'
    
    def __init__(self, name=RawText(), text=RawText(), in_out=None):
        RawParam.__init__(self, name, text)

    def __str__(self):
        return 'RawTParam(%s, %s)' % (repr(self.name.text), repr(self.text.text))

    def getFormatted(self, formatter):
        return formatter.formatCommand('tparam', self.text.text, self.name.text)


class RawReturn(RawParam):
    """RawDoc for one return description.
    
    @ivar type The return type.
    @ivar text RawParagraph entry with the documentation of the parameter.
    """
    # TODO(holtgrew): I don't see the return type here :)
    
    def __init__(self, name=RawText(), text=RawText(), in_out=None):
        RawParam.__init__(self, name, text)

    def getType(self):
        return 'return'

    def getFormatted(self, formatter):
        return formatter.formatCommand('return', self.text.text, self.name.text)


class RawSignature(object):
    """A representation of a @signature entry.
    
    @ivar value The @signature's clauses's parameter. RawText.
    """
    
    def __init__(self, text):
        self.text = text

    def getType(self):
        return 'signature'

    def __str__(self):
        return 'RawSignature(%s)' % repr(self.text.text)

    def getFormatted(self, formatter):
        return formatter.formatCommand('signature', self.text.text.strip())

