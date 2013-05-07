#!/usr/bin/env python
"""Parser for SeqAn Doxygen dialect.

The Doxygen documentation has a regular grammar and thus can be parsed quite
easily.  We do not use a parser generator such as PLY since it is a bit
overkill (it can generate parser for some context-free grammars) and not
as straightforward to use as writing a simple parser by hand.
"""

import operator
import itertools
import re
import sys

import raw_doc
import dox_tokens


class ParserError(Exception):
    """Raised when there is a parser error."""
    
    def __init__(self, token=None, msg=''):
        if msg and token:
            args = (token.lineno, token.column, repr(token.val), msg)
            message = 'Parser error at %d:%d ("%s") : %s' % args
        elif not msg and token:
            args = (token.lineno, token.column, repr(token.val))
            message = 'Parser error at %d:%d ("%s").' % args
        else:
            message = 'Parse error: %s' % (msg,)
        Exception.__init__(self, message)
        self.msg = msg
        self.token = token


def stripWhitespaceTokens(token_list, strip_lt_breaks=False):
    """Strip leading and trailing whitespace tokens from token_list."""
    types = ['SPACE']
    if strip_lt_breaks:
        types.append('BREAK')
    while token_list and token_list[0].type in types:
        token_list.pop(0)
    while token_list and token_list[-1].type in types:
        token_list.pop()


def normalizeWhitespaceTokens(token_list, strip_lt_breaks=False):
    """Normalize whitespace by replacing multiple consecutive spaces by one.
    """
    positions = [i for i, t in enumerate(token_list) if t.type == 'SPACE']
    for i in reversed(positions):
        token_list[i].val = ' '
    stripWhitespaceTokens(token_list, strip_lt_breaks)


class GenericSimpleClauseState(object):
    """Handler used in *DocState for handling simple text clauses clauses.
    """
    
    def __init__(self, parser, parent):
        self.parser = parser
        self.parent = parent
        self.tokens = []
        self.entry_class = None
        # Whether or not to strip leading and trailing breaks.
        self.strip_lt_breaks = False
        # Whether to normalize whitespace tokens in getEntry().
        self.normalize_tokens = True

    def getEntry(self):
        """Returns the Entry for the brief clause."""
        if self.normalize_tokens:
            normalizeWhitespaceTokens(self.tokens, self.strip_lt_breaks)
        return self.entry_class(raw_doc.RawText(self.tokens))

    def handle(self, token):
        # One or more empty lines end such a clause as well as another
        # clause-starting command.
        if token.type in ['EMPTYLINE', 'EOF']:
            self.parent.endClause()
        elif token.type in dox_tokens.CLAUSE_STARTING or \
             token.type in dox_tokens.ITEM_STARTING:
            self.parent.endClause(token)
        elif token.type == 'SPACE' and (not self.tokens or self.tokens[-1].type == 'BREAK'):
            return  # Skip space at beginning or after break
        elif token.type == 'BREAK' and (self.tokens and self.tokens[-1].type == 'SPACE'):
            self.tokens[-1] = token  # Replace space before break
        else:
            #print 'APPEND %s' % repr(token.val)
            self.tokens.append(token)


class ParagraphState(GenericSimpleClauseState):
    """Handler used in *DocState for handling tokens of a paragraph."""
    
    def __init__(self, parser, parent):
        GenericSimpleClauseState.__init__(self, parser, parent)
        self.entry_class = raw_doc.RawParagraph


class SignatureState(GenericSimpleClauseState):
    """Handler used in *DocState for handling tokens of a paragraph."""
    
    def __init__(self, parser, parent):
        GenericSimpleClauseState.__init__(self, parser, parent)
        self.entry_class = raw_doc.RawSignature


class CodeState(GenericSimpleClauseState):
    """Handler used in *DocState for handling tokens of a paragraph."""
    
    def __init__(self, parser, parent):
        GenericSimpleClauseState.__init__(self, parser, parent)
        self.entry_class = raw_doc.RawCode
        self.normalize_tokens = False

    def handle(self, token):
        # Code is only ended by @endcode.
        if token.type == 'COMMAND_ENDCODE':
            self.parent.endClause()
        else:
            self.tokens.append(token)


class DeprecatedState(GenericSimpleClauseState):
    """Handler for the @deprecated clause."""

    def __init__(self, parser, parent):
        GenericSimpleClauseState.__init__(self, parser, parent)
        self.entry_class = raw_doc.RawDeprecated


class HeaderfileState(GenericSimpleClauseState):
    """Handler for the @headerfile clause."""

    def __init__(self, parser, parent):
        GenericSimpleClauseState.__init__(self, parser, parent)
        self.entry_class = raw_doc.RawHeaderfile


class ImplementsState(GenericSimpleClauseState):
    """Handler for the @implements clause."""

    def __init__(self, parser, parent):
        GenericSimpleClauseState.__init__(self, parser, parent)
        self.entry_class = raw_doc.RawImplements


class ExtendsState(GenericSimpleClauseState):
    """Handler for the @extends clause."""

    def __init__(self, parser, parent):
        GenericSimpleClauseState.__init__(self, parser, parent)
        self.entry_class = raw_doc.RawExtends
        self.strip_lt_breaks = True


class BriefState(GenericSimpleClauseState):
    """Handler for the @brief clause."""

    def __init__(self, parser, parent):
        GenericSimpleClauseState.__init__(self, parser, parent)
        self.entry_class = raw_doc.RawBrief


class SeeState(GenericSimpleClauseState):
    """Handler for the @see clause."""

    def __init__(self, parser, parent):
        GenericSimpleClauseState.__init__(self, parser, parent)
        self.entry_class = raw_doc.RawSee


class ParamState(GenericSimpleClauseState):
    """Handler used in *DocState for handling @param clauses.
    """
    
    def __init__(self, parser, parent):
        GenericSimpleClauseState.__init__(self, parser, parent)
        self.entry_class = raw_doc.RawParam
        # Parameters need more than one token list.
        self.in_out = None
        self.name = None

    def handle(self, token):
        # Special handling of the in/out token and the name, if any.
        #print '.... TOKEN', token, token.type == 'PARAM_IN_OUT'
        if token.type == 'PARAM_IN_OUT':
            self.in_out = token
        elif self.name is None:
            if token.type != 'SPACE':
                self.name = token
            # Skipping whitespace.
        else:
            GenericSimpleClauseState.handle(self, token)

    def getEntry(self):
        """Returns the Entry for the parameter."""
        normalizeWhitespaceTokens(self.tokens)
        return self.entry_class(raw_doc.RawText([self.name]),
                                raw_doc.RawText(self.tokens), self.in_out)


class TParamState(ParamState):
    """Handler used in *DocState for handling @tparam clauses.
    """
    
    def __init__(self, parser, parent):
        ParamState.__init__(self, parser, parent)
        self.entry_class = raw_doc.RawTParam


class ReturnState(ParamState):
    """Handler used in *DocState for handling @return clauses.
    """
    
    def __init__(self, parser, parent):
        ParamState.__init__(self, parser, parent)
        self.entry_class = raw_doc.RawReturn


class SectionState(object):
    """Handler used in *DocState for handling @section clauses.
    """
    
    def __init__(self, parser, parent):
        self.parser = parser
        self.parent = parent
        self.tokens = []
        self.level = 0

    def getEntry(self):
        """Returns the Entry for the template parameter."""
        return raw_doc.RawSection(raw_doc.RawText(self.tokens), self.level)

    def handle(self, token):
        # One or more empty lines end a @section raw_doc.Raw
        if token.type in ['EMPTYLINE', 'EOF']:
            self.parent.endClause()
        elif token.type in dox_tokens.CLAUSE_STARTING:
            self.parent.endClause(token)
        else:
            self.tokens.append(token)


class SubsectionState(SectionState):
    """Handler used in *DocState for handling @subsection clauses.
    """
    
    def __init__(self, parser, parent):
        SectionState.__init__(self, parser, parent)
        self.level = 1


class IncludeState(object):
    """Handler used in *DocState for handling @include clauses.
    """
    
    def __init__(self, parser, parent):
        self.parser = parser
        self.parent = parent
        self.tokens = []

    def getEntry(self):
        return raw_doc.RawInclude(raw_doc.RawText(self.tokens))

    def handle(self, token):
        if token.type in dox_tokens.LINE_BREAKS or token.type == 'EOF':
            self.parent.endClause()
        elif token.type in dox_tokens.CLAUSE_STARTING:
            self.parent.endClause(token)
        else:
            if token.type != 'SPACE' or self.tokens:
                self.tokens.append(token)


class SnippetState(object):
    """Handler used in *DocState for handling @snippet clauses.
    """
    
    def __init__(self, parser, parent):
        self.parser = parser
        self.parent = parent
        self.path_done = False  # True after first space after path
        self.path_tokens = []
        self.name_tokens = []

    def getEntry(self):
        return raw_doc.RawSnippet(raw_doc.RawText(self.path_tokens),
                                  raw_doc.RawText(self.name_tokens))

    def handle(self, token):
        if token.type in dox_tokens.LINE_BREAKS or token.type == 'EOF':
            self.parent.endClause()
        elif token.type in dox_tokens.CLAUSE_STARTING:
            self.parent.endClause(token)
        else:
            if token.type == 'SPACE':
                if not self.path_tokens:
                    return  # Ignore beginning space.
                elif self.path_tokens and not self.name_tokens:
                    self.path_done = True
                else:  # append to name
                    self.name_tokens.append(token)
            else:
                if self.path_done:
                    self.name_tokens.append(token)
                else:
                    self.path_tokens.append(token)


class TopLevelState(object):
    """Top level state, expecting a command starting a documentation item.
    
    Whitespace is skipped.
    """

    def __init__(self, parser):
        self.parser = parser

    def handle(self, token):
        # We ignore whitespace.
        if token.type in dox_tokens.WHITESPACE:
            return
        # We expect an item-starting token.
        if token.type in dox_tokens.ITEM_STARTING:
            state_map = {'COMMAND_CLASS' :       'class',
                         'COMMAND_CONCEPT' :     'concept',
                         'COMMAND_FUNCTION' :    'function',
                         'COMMAND_MACRO':        'macro',
                         'COMMAND_METAFUNCTION': 'metafunction',
                         'COMMAND_PAGE':         'page',
                         'COMMAND_DEFGROUP':     'group',
                         'COMMAND_VARIABLE':     'var',
                         'COMMAND_TAG':          'tag',
                         'COMMAND_TYPEDEF':      'typedef',
                         'COMMAND_ADAPTION':     'adaption',
                         'COMMAND_ENUM':         'enum',}
            self.parser.enterState(state_map[token.type])
            return
        msg = 'Expecting one of {%s.} but is %s' % (", ".join(dox_tokens.ITEM_STARTING), token.type)
        raise ParserError(token, msg)


class GenericDocState(object):
    """Generic token handler for top-level entries.
    """

    def __init__(self, parser, entry_class, state_name):
        self.parser = parser
        # Substate is one of "first_line" or "body".
        self.substate = 'first_line'
        self.state_name = state_name
        self.first_line_tokens = None
        self.entry_class = entry_class
        self.clause_state = None
        self.allowed_commands = dox_tokens.CLAUSE_STARTING

    def getEntry(self):
        return self.entry

    def entered(self):
        """Called when the state is entered successfully with @class etc.
        """
        #print >>sys.stderr, ">>>>ENTERING CLASS STATE"
        self.first_line_tokens = []
        self.substate = 'first_line'
        self.entry = self.entry_class()

    def handle(self, token):
        #print 'state = class, substate = %s, clause_state %s' % (self.substate, self.clause_state)
        if self.substate == 'first_line':
            # If we have a line break in the first line then we go to the body
            # of the class documentation.
            if token.type in dox_tokens.LINE_BREAKS or token.type == 'EOF':
                normalizeWhitespaceTokens(self.first_line_tokens)
                self.entry.name = raw_doc.RawText(self.first_line_tokens)
                self.substate = 'body'
                return
            # Otherwise, we collect the token's value.
            self.first_line_tokens.append(token)
        elif self.substate == 'body':
            # If we are already in a clause then continue to use the state.
            if not self.clause_state is None:
                self.clause_state.handle(token)
                return
            # In the body, we look for a item-generating token and will enter
            # the corresponding state if we hit it.
            if token.type in dox_tokens.ITEM_STARTING:
                # Leave the top-level class state and handle the token with the
                # top-level handler.  This will create the appropriate state.
                self.parser.leaveState(self.state_name)
                assert self.parser.states[0] == 'top'
                self.parser.handleToken(token)
                return
            # In a class documentation body, a clause-starting command triggers
            # parsing of that clause.
            if token.type in dox_tokens.CLAUSE_STARTING:
                state_map = {'COMMAND_SIGNATURE' : SignatureState(self.parser, self),
                             'COMMAND_CODE' : CodeState(self.parser, self),
                             'COMMAND_BRIEF' : BriefState(self.parser, self),
                             'COMMAND_EXTENDS' : ExtendsState(self.parser, self),
                             'COMMAND_HEADERFILE' : HeaderfileState(self.parser, self),
                             'COMMAND_DEPRECATED' : DeprecatedState(self.parser, self),
                             'COMMAND_IMPLEMENTS' : ImplementsState(self.parser, self),
                             'COMMAND_SEE' : SeeState(self.parser, self),
                             'COMMAND_RETURN' : ReturnState(self.parser, self),
                             'COMMAND_PARAM' : ParamState(self.parser, self),
                             'COMMAND_TPARAM' : TParamState(self.parser, self),
                             'COMMAND_SECTION' : SectionState(self.parser, self),
                             'COMMAND_SNIPPET' : SnippetState(self.parser, self),
                             'COMMAND_SUBSECTION' : SubsectionState(self.parser, self),
                             'COMMAND_INCLUDE' : IncludeState(self.parser, self),}
                if not token.type in self.allowed_commands:
                    raise ParserError(msg='Invalid command %s' % token.val)
                self.clause_state = state_map[token.type]
                #print '>>> SWITCHING TO CLAUSE STATE %s' % self.clause_state
                return
            # Some commands are explicitely marked as non-paragraph, such as
            # the @endcode token.  These are invalid tokens.
            if token.type in dox_tokens.NON_PARAGRAPH:
                raise ParserError(token, 'Invalid command!')
            # Any other token is an inline-token and part of a paragraph that
            # we will start.  The same is true for any word.
            self.clause_state = ParagraphState(self.parser, self)
            self.clause_state.handle(token)
        else:
            raise ParserError(msg='Invalid substate in @class!')

    def endClause(self, token=None):
        """Ends the current clause and appends its result to the documentation."""
        #print '>>> END CLAUSE(%s)' % token
        if self.clause_state.getEntry():
            entry = self.clause_state.getEntry()
            if entry.getType() in ['paragraph', 'section', 'include', 'code', 'snippet']:
                self.entry.addParagraph(entry)
            elif entry.getType() == 'signature':
                self.entry.addSignature(entry)
            elif entry.getType() == 'tparam':
                self.entry.addTParam(entry)
            elif entry.getType() == 'brief':
                self.entry.addBrief(entry)
            elif entry.getType() == 'param':
                self.entry.addParam(entry)
            elif entry.getType() == 'see':
                self.entry.addSee(entry)
            elif entry.getType() == 'return':
                self.entry.addReturn(entry)
            elif entry.getType() == 'extends':
                self.entry.addExtends(entry)
            elif entry.getType() == 'implements':
                self.entry.addImplements(entry)
            elif entry.getType() == 'headerfile':
                self.entry.addHeaderfile(entry)
            elif entry.getType() == 'deprecated':
                self.entry.addDeprecationMsg(entry)
            else:
                assert False, '%s' % entry
        self.clause_state = None
        if token:
            self.handle(token)


class ClassDocState(GenericDocState):
    """State for documentation of a class."""
    
    def __init__(self, parser):
        GenericDocState.__init__(self, parser, raw_doc.RawClass, 'class')
        self.allowed_commands = set(['COMMAND_SIGNATURE', 'COMMAND_CODE',
                                     'COMMAND_SEE', 'COMMAND_BRIEF', 'COMMAND_TPARAM',
                                     'COMMAND_SECTION', 'COMMAND_SUBSECTION',
                                     'COMMAND_INCLUDE', 'COMMAND_SNIPPET', 'COMMAND_EXTENDS',
                                     'COMMAND_IMPLEMENTS', 'COMMAND_HEADERFILE',
                                     'COMMAND_DEPRECATED'])

class FunctionDocState(GenericDocState):
    """State for documentation of a function."""
    
    def __init__(self, parser):
        GenericDocState.__init__(self, parser, raw_doc.RawFunction, 'function')
        self.allowed_commands = set(['COMMAND_SIGNATURE', 'COMMAND_CODE',
                                     'COMMAND_SEE', 'COMMAND_BRIEF', 'COMMAND_TPARAM',
                                     'COMMAND_PARAM',
                                     'COMMAND_SECTION', 'COMMAND_SUBSECTION',
                                     'COMMAND_INCLUDE', 'COMMAND_SNIPPET', 'COMMAND_RETURN',
                                     'COMMAND_HEADERFILE', 'COMMAND_DEPRECATED'])


class MacroDocState(GenericDocState):
    """State for documentation of a macro."""
    
    def __init__(self, parser):
        GenericDocState.__init__(self, parser, raw_doc.RawMacro, 'macro')
        self.allowed_commands = set(['COMMAND_SIGNATURE', 'COMMAND_CODE',
                                     'COMMAND_SEE', 'COMMAND_BRIEF', 'COMMAND_PARAM',
                                     'COMMAND_SECTION', 'COMMAND_SUBSECTION',
                                     'COMMAND_INCLUDE', 'COMMAND_SNIPPET', 'COMMAND_RETURN',
                                     'COMMAND_HEADERFILE', 'COMMAND_DEPRECATED'])


class MetafunctionDocState(GenericDocState):
    """State for documentation of a metafunction."""
    
    def __init__(self, parser):
        GenericDocState.__init__(self, parser, raw_doc.RawMetafunction, 'metafunction')
        self.allowed_commands = set(['COMMAND_SIGNATURE', 'COMMAND_CODE',
                                     'COMMAND_SEE', 'COMMAND_BRIEF', 'COMMAND_TPARAM',
                                     'COMMAND_SECTION', 'COMMAND_SUBSECTION',
                                     'COMMAND_INCLUDE', 'COMMAND_SNIPPET', 'COMMAND_RETURN',
                                     'COMMAND_HEADERFILE', 'COMMAND_DEPRECATED'])


class ConceptDocState(GenericDocState):
    """State for documentation of a concept."""
    
    def __init__(self, parser):
        GenericDocState.__init__(self, parser, raw_doc.RawConcept, 'concept')
        self.allowed_commands = set(['COMMAND_SIGNATURE', 'COMMAND_CODE',
                                     'COMMAND_SEE', 'COMMAND_BRIEF',
                                     'COMMAND_SECTION', 'COMMAND_SUBSECTION',
                                     'COMMAND_INCLUDE', 'COMMAND_SNIPPET', 'COMMAND_EXTENDS',
                                     'COMMAND_HEADERFILE', 'COMMAND_DEPRECATED'])


class PageState(GenericDocState):
    """State for a documentation page."""
        
    def __init__(self, parser):
        GenericDocState.__init__(self, parser, raw_doc.RawPage, 'page')
        self.allowed_commands = set(['COMMAND_CODE',
                                     'COMMAND_SEE', 'COMMAND_BRIEF',
                                     'COMMAND_SECTION', 'COMMAND_SUBSECTION',
                                     'COMMAND_INCLUDE', 'COMMAND_SNIPPET'])

    def entered(self):
        GenericDocState.entered(self)
        self.name_read = False
        self.name_tokens = []
        self.title_tokens = []

    def handle(self, token):
        # Handle first state here and the remaining in the parent class.
        #print >>sys.stderr, token.type, repr(token.val), self.name_read
        if self.substate == 'first_line':
            # If we have a line break in the first line then we go to the body
            # of the class documentation.
            if token.type in dox_tokens.LINE_BREAKS or token.type == 'EOF':
                #print >>sys.stderr, [v.val for v in self.name_tokens], [v.val for v in self.type_tokens]
                normalizeWhitespaceTokens(self.name_tokens)
                normalizeWhitespaceTokens(self.title_tokens)
                self.entry.name = raw_doc.RawText(self.name_tokens)
                self.entry.title = raw_doc.RawText(self.title_tokens)
                self.substate = 'body'
                return
            # Skip space at the beginning of the type.
            if not self.name_tokens and token.type == 'SPACE':
                return
            # Otherwise, we collect the token's value.
            if self.name_read:
                #print >>sys.stderr, 'NAME', token.type, repr(token.val)
                self.title_tokens.append(token)
            else:
                if token.type == 'SPACE':
                    self.name_read = True
                else:
                    #print >>sys.stderr, 'TYPE', token.type, repr(token.val)
                    self.name_tokens.append(token)
        else:
            GenericDocState.handle(self, token)


class GroupState(PageState):
    def __init__(self, parser):
        GenericDocState.__init__(self, parser, raw_doc.RawGroup, 'group')
        self.allowed_commands = set(['COMMAND_CODE',
                                     'COMMAND_SEE', 'COMMAND_BRIEF',
                                     'COMMAND_SECTION', 'COMMAND_SUBSECTION',
                                     'COMMAND_INCLUDE', 'COMMAND_SNIPPET'])


class VariableState(GenericDocState):
    """State for a variable."""
        
    def __init__(self, parser):
        GenericDocState.__init__(self, parser, raw_doc.RawVariable, 'var')
        self.allowed_commands = set(['COMMAND_SIGNATURE', 'COMMAND_CODE',
                                     'COMMAND_SEE', 'COMMAND_BRIEF',
                                     'COMMAND_SECTION', 'COMMAND_SUBSECTION',
                                     'COMMAND_INCLUDE', 'COMMAND_SNIPPET',
                                     'COMMAND_HEADERFILE', 'COMMAND_DEPRECATED'])

    def entered(self):
        GenericDocState.entered(self)
        self.type_read = False
        self.type_tokens = []
        self.name_tokens = []

    def handle(self, token):
        # Handle first state here and the remaining in the parent class.
        #print >>sys.stderr, token.type, repr(token.val), self.type_read
        if self.substate == 'first_line':
            # If we have a line break in the first line then we go to the body
            # of the class documentation.
            if token.type in dox_tokens.LINE_BREAKS or token.type == 'EOF':
                #print >>sys.stderr, [v.val for v in self.name_tokens], [v.val for v in self.type_tokens]
                normalizeWhitespaceTokens(self.name_tokens)
                normalizeWhitespaceTokens(self.type_tokens)
                self.entry.name = raw_doc.RawText(self.name_tokens)
                self.entry.type = raw_doc.RawText(self.type_tokens)
                self.substate = 'body'
                return
            # Skip space at the beginning of the type.
            if not self.type_tokens and token.type == 'SPACE':
                return
            # Otherwise, we collect the token's value.
            if self.type_read:
                #print >>sys.stderr, 'NAME', token.type, repr(token.val)
                self.name_tokens.append(token)
            else:
                if token.type == 'SPACE':
                    self.type_read = True
                else:
                    #print >>sys.stderr, 'TYPE', token.type, repr(token.val)
                    self.type_tokens.append(token)
        else:
            GenericDocState.handle(self, token)



class TagState(GenericDocState):
    def __init__(self, parser):
        GenericDocState.__init__(self, parser, raw_doc.RawTag, 'tag')
        self.allowed_commands = set(['COMMAND_SIGNATURE', 'COMMAND_CODE',
                                     'COMMAND_SEE', 'COMMAND_BRIEF',
                                     'COMMAND_SECTION', 'COMMAND_SUBSECTION',
                                     'COMMAND_INCLUDE', 'COMMAND_SNIPPET',
                                     'COMMAND_HEADERFILE', 'COMMAND_DEPRECATED'])


class EnumState(GenericDocState):
    """State for an enum."""
        
    def __init__(self, parser):
        GenericDocState.__init__(self, parser, raw_doc.RawEnum, 'enum')
        self.allowed_commands = set(['COMMAND_SIGNATURE', 'COMMAND_CODE',
                                     'COMMAND_SEE', 'COMMAND_BRIEF',
                                     'COMMAND_SECTION', 'COMMAND_SUBSECTION',
                                     'COMMAND_INCLUDE', 'COMMAND_SNIPPET',
                                     'COMMAND_HEADERFILE', 'COMMAND_DEPRECATED'])


class AdaptionState(GenericDocState):
    """State for an adaption."""
        
    def __init__(self, parser):
        GenericDocState.__init__(self, parser, raw_doc.RawAdaption, 'adaption')
        self.allowed_commands = set(['COMMAND_SIGNATURE', 'COMMAND_CODE',
                                     'COMMAND_SEE', 'COMMAND_BRIEF',
                                     'COMMAND_SECTION', 'COMMAND_SUBSECTION',
                                     'COMMAND_INCLUDE', 'COMMAND_SNIPPET',
                                     'COMMAND_HEADERFILE', 'COMMAND_DEPRECATED'])


class TypedefState(GenericDocState):
    """State for an typedef."""
        
    def __init__(self, parser):
        GenericDocState.__init__(self, parser, raw_doc.RawTypedef, 'typedef')
        self.allowed_commands = set(['COMMAND_SIGNATURE', 'COMMAND_CODE',
                                     'COMMAND_SEE', 'COMMAND_BRIEF',
                                     'COMMAND_SECTION', 'COMMAND_SUBSECTION',
                                     'COMMAND_INCLUDE', 'COMMAND_SNIPPET',
                                     'COMMAND_HEADERFILE', 'COMMAND_DEPRECATED'])


class Parser(object):
    """The parser class takes tokens from a lexer.Lexer class.
    
    It generates a raw_doc.RawDoc object from this with
    raw_doc.RawEntry objects.
    """

    def __init__(self):
        self.states = ['top']
        self.handlers = {
            'top':          TopLevelState(self),
            'class':        ClassDocState(self),
            'function':     FunctionDocState(self),
            'macro':        MacroDocState(self),
            'metafunction': MetafunctionDocState(self),
            'concept':      ConceptDocState(self),
            'page':         PageState(self),
            'group':        GroupState(self),
            'var':          VariableState(self),
            'tag':          TagState(self),
            'enum':         EnumState(self),
            'adaption':     AdaptionState(self),
            'typedef':      TypedefState(self),
            }
        self.documentation = raw_doc.RawDoc()

    def parse(self, lexer):
        for token in lexer.tokens():
            self.handleToken(token)
        while len(self.states) > 1:
            self.leaveState(self.states[-1])

    def handleToken(self, token):
        #print 'Handling %s in states %s' % (token, self.states)
        self.handlers[self.states[-1]].handle(token)

    def enterState(self, state):
        #print 'entering state %s' % state
        self.states.append(state)
        self.handlers[state].entered()

    def leaveState(self, state):
        #print 'leaving state %s' % state
        if self.states[-1] == state:
            if self.handlers[state].getEntry():
                self.documentation.addEntry(self.handlers[state].getEntry())
            return self.states.pop()
        raise ParserError(msg='Invalid state "%s", expecting "%s"!' % (state, self.states[-1]))
