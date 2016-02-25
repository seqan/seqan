#!/usr/bin/env python2
"""Parser for the signature supported by the SeqAn Doxygen-style documentation.
"""

# TODO(holtgrew): The parser has become quite complex. Maybe using some external library for parsing is in order.

import sys

import lexer


TOKENS = (
    ('KWD_TEMPLATE',   r'template'),
    ('KWD_TYPENAME',   r'typename'),
    ('KWD_CLASS',      r'class'),
    ('KWD_CONCEPT',    r'concept'),
    ('KWD_STRUCT',     r'struct'),
    ('KWD_ENUM',       r'enum'),

    ('IDENTIFIER',     r'[a-zA-Z_~][a-zA-Z_0-9~]*'),
    ('COMMA',          r','),
    ('NAME_SEP',       r'::'),
    ('HASH',           r'#'),
    ('SEMICOLON',      r';'),
    ('SPACE',          r'[\t ]'),
    ('PROUND_OPEN',    r'\('),
    ('PROUND_OPEN',    r'\('),
    ('PROUND_CLOSE',   r'\)'),
    ('PANGULAR_OPEN' , r'<'),
    ('PANGULAR_CLOSE', r'>'),
    )


class SigParseException(Exception):
    """Raised in the case of signature parsing error."""

    def __init__(self, msg, line=0, column=0):
        Exception.__init__(self, msg)
        self.line = line
        self.column = column
    

class Arg(object):
    """
    @ivar type: The type of the template parameter, e.g. 'typename',
                'class', 'int', 'unsigned' etc. str
    @ivar name: The name of the parameter. str
    """

    def __init__(self, type=None, name=None):
        self.type = type
        self.name = name


class SigEntry(object):
    """A signature entry.

    The following kinds are possible: concept, class, function, variable,
    enum, struct.

    @ivar name: Name of the element.
    @ivar kind: The kind of the element.
    @ivar params: Parameters of the function, in case of function.
    @ivar is_tpl: Whether or not the entry is a template.
    @ivar tparams: Template parameters, in case of templates.
    @ivar return_type: The name of the return type, in case of function.
    @ivar return_name: Name after the :: for Metafunctions.
    @ivar var_type: The type of the variable.

    """

    def __init__(self, name=None, kind=None, params=[], tparams=[],
                 is_tpl=False, return_type=None, return_name=None,
                 var_type=None):
        self.name = name
        self.kind = kind
        self.params = list([])
        self.tparams = list([])
        self.is_tpl = is_tpl
        self.return_type = return_type
        self.return_name = return_name
        self.var_type = var_type

    def toString(self):
        """Convert the SigEntry object back into a string."""
        types = ['concept', 'class', 'struct', 'enum']
        if not self.is_tpl and self.kind in types:
            return '%s %s;' % (self.kind, self.name)
        elif not self.is_tpl and self.kind == 'function':
            params = ', '.join(['%s %s' % (p.type, p.name) for p in self.params])
            if self.return_type:
                return '%s %s(%s);' % (self.return_type, self.name, params)
            else:
                return '%s(%s);' % (self.name, params)
        elif self.is_tpl and self.kind == 'function':
            tparams = ', '.join(['%s %s' % (p.type, p.name) for p in self.tparams])
            params = ', '.join(['%s %s' % (p.type, p.name) for p in self.params])
            return 'template <%s>\n%s %s(%s);' % (tparams, self.return_type, self.name, params)
        elif self.is_tpl and self.kind in ['struct', 'class']:
            tparams = ', '.join(['%s %s' % (p.type, p.name) for p in self.tparams])
            params = ', '.join(['%s %s' % (p.type, p.name) for p in self.params])
            return 'template <%s>\n%s %s;' % (tparams, self.kind, self.name)
        elif self.kind == 'metafunction':
            tparams = ', '.join([p.name for p in self.tparams])
            if self.return_type:
                return '%s %s<%s>::%s;' % (self.return_type, self.name, tparams,
                                           self.return_name)
            else:
                return '%s<%s>::%s;' % (self.name, tparams, self.return_name)
        elif self.kind == 'variable':
            return '%s %s;' % (self.var_type, self.name)


class SigParser(object):
    def __init__(self, buffer):
        self.buffer = buffer
        self.lexer = lexer.Lexer(TOKENS)
        self.lexer.input(buffer)
        self.tokens = self.lexer.tokens()

    def parseTemplate(self, token):
        tparams = []
        # Read <
        t = self.tokens.next()
        self.expectNotEof(t)
        if t.type != 'PANGULAR_OPEN':
            raise SigParseException('Expected opening angular parenthesis')
        # Parse template parameters.
        self.parseParams('PANGULAR_CLOSE', tparams,
                         ['IDENTIFIER', 'KWD_TYPENAME', 'KWD_CLASS'])
        # Parse remaining.
        sig_entry = self.parse()
        sig_entry.is_tpl = True
        sig_entry.tparams = tparams
        return sig_entry

    def parseParams(self, end_token, params_dest, type_tokens):
        t = self.tokens.next()
        self.expectNotEof(t)
        while t.type != end_token:
            if t.type not in type_tokens:
                raise SigParseException('Expected identifier got "%s"' % t.val)
            arg = Arg(type=t.val)
            t = self.tokens.next()
            self.expectNotEof(t)
            if t.type != 'IDENTIFIER':
                raise SigParseException('Expected identifier got "%s"' % t.val)
            arg.name = t.val
            t = self.tokens.next()
            self.expectNotEof(t)
            if t.type not in ['COMMA', end_token]:
                raise SigParseException('Expected COMMA or closing parenthesis')
            if t.type != end_token:
                t = self.tokens.next()
            params_dest.append(arg)

    def parseMetafunctionType(self, name):
        sig_entry = SigEntry(kind='metafunction')
        sig_entry.name = name
        # Expect "#$name" or PANGULAR_CLOSE
        t = self.tokens.next()
        self.expectNotEof(t)
        if t.type == 'HASH':
            sig_entry.name += '#'
            t = self.tokens.next()
            self.expectNotEof(t)
            if t.type != 'IDENTIFIER':
                raise SigParseException('Expecting identifier')
            sig_entry.name += t.val
            t = self.tokens.next()
            self.expectNotEof(t)
        while t.type != 'PANGULAR_CLOSE':
            if t.type != 'IDENTIFIER':
                raise SigParseException('Expecting identifier')
            arg = Arg(name=t.val)
            sig_entry.tparams.append(arg)
            # Read "," or ">"
            t = self.tokens.next()
            self.expectNotEof(t)
            if t.type not in ['PANGULAR_CLOSE', 'COMMA']:
                raise SigParseException('Expecting ">" or ","')
            if t.type == 'COMMA':
                t = self.tokens.next()
                self.expectNotEof(t)
        # Expect "::"
        t = self.tokens.next()
        self.expectNotEof(t)
        if t.type != 'NAME_SEP':
            raise SigParseException('Expecting "::"')
        # Read return_name
        t = self.tokens.next()
        self.expectNotEof(t)
        if t.type != 'IDENTIFIER':
            raise SigParseException('Expecting identifier got %s' % repr(t.val))
        sig_entry.return_name = t.val
        return sig_entry

    def parseMetafunctionValue(self, return_type, name):
        sig_entry = self.parseMetafunctionType(name)
        sig_entry.return_type = return_type
        return sig_entry

    def parseFunction(self, token):
        """Parse a function, variable, or metafunction.

        We started out with an identifier.  The things that this function will
        be triggered for:

          TReturn name(T1 x1, T2 x2, ...)
          TReturn Klass::name(T1 x1, T2 x2, ...)
          TReturn Klass#name(T1 x1, T2 x2, ...)
          TReturn Name<TParam>::VALUE
          TReturn Klass#Name<TParam>::VALUE
          Name<TParam>::Type
          Klass#Name<TParam>::Type
          T var
        
        @param token: lexer.Token object with the previous token.

        """
        is_constructor = False
        other_name = token.val
        # get next token, i sname or "<"
        t = self.tokens.next()
        self.expectNotEof(t)
        if t.type in ['HASH', 'NAME_SEP']:
            other_name += t.val
            t = self.tokens.next()
            self.expectNotEof(t)
            if t.type != 'IDENTIFIER':
                raise SigParseException('Expecting identifier.')
            other_name += t.val
            t = self.tokens.next()
            self.expectNotEof(t)
        if t.type == 'PANGULAR_OPEN':
            return self.parseMetafunctionType(other_name)
        if t.type == 'PROUND_OPEN':
            is_constructor = True
        elif t.type != 'IDENTIFIER':
            raise SigParseException('Expecting identifier as function name.')
        name = t.val
        if not is_constructor:
            t = self.tokens.next()
            self.expectNotEof(t)
            if t.type in ['HASH', 'NAME_SEP']:
                name += t.val
                t = self.tokens.next()
                self.expectNotEof(t)
                if t.type != 'IDENTIFIER':
                    raise SigParseException('Expecting identifier.')
                name += t.val
                t = self.tokens.next()
                self.expectNotEof(t)
            # expect "(" or "<"
            if t.type == 'PANGULAR_OPEN':
                return self.parseMetafunctionValue(other_name, name)
        # Expecting <eof>, ";", or (.  The last triggers generation of a
        # function, the other of a variable.
        if t.type in ['EOF', 'SEMICOLON']:
            sig_entry = SigEntry(kind='variable')
            sig_entry.var_type = other_name
            sig_entry.name = name
            return sig_entry
        sig_entry = SigEntry(kind='function')
        if is_constructor:
            sig_entry.return_type = None
            sig_entry.name = other_name
        else:
            sig_entry.return_type = other_name
            sig_entry.name = name
        if t.type in ['HASH', 'NAME_SEP']:
            sig_entry.name += t.val
            t = self.tokens.next()
            self.expectNotEof(t)
            if t.type != 'IDENTIFIER':
                raise SigParseException('EOF not expected')
            sig_entry.name += t.val
            t = self.tokens.next()
            self.expectNotEof(t)
        if t.type != 'PROUND_OPEN':
            raise SigParseException('Expecting opening parenthesis after '
                                    'function name')
        # parse parameters
        self.parseParams('PROUND_CLOSE', sig_entry.params, ['IDENTIFIER'])
        return sig_entry

    def parseCSE(self, kind):
        """Parse class, struct, enum."""
        sig_entry = SigEntry(kind=kind)
        t = self.tokens.next()
        self.expectNotEof(t)
        if t.type != 'IDENTIFIER':
            raise SigParseException('Expecting identifier after "%s"!' % kind)
        sig_entry.name = t.val
        return sig_entry
        
    def parseClass(self, token):
        return self.parseCSE('class')

    def parseConcept(self, token):
        return self.parseCSE('concept')

    def parseStruct(self, token):
        return self.parseCSE('struct')

    def parseEnum(self, token):
        return self.parseCSE('enum')

    def expectNotEof(self, token):
        if token.type == 'EOF':
            raise SigParseException('Unexpecte EOF!')
        
    def parse(self):
        try:
            t = self.tokens.next()
            self.expectNotEof(t)
            m = {'KWD_TEMPLATE':   self.parseTemplate,
                 'KWD_CLASS':      self.parseClass,
                 'KWD_CONCEPT':    self.parseConcept,
                 'KWD_STRUCT':     self.parseStruct,
                 'KWD_ENUM':       self.parseEnum,
                 'IDENTIFIER':     self.parseFunction}
            if not t.type in m:
                raise SigParseException('Unexpected token of type %s' % t.type)
            return m[t.type](t)
        except lexer.LexerError, e:
            raise SigParseException('Lexer error: %s at pos %s when parsing %s' % (e, e.pos, self.buffer))
