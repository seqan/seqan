#!/usr/bin/env python2
"""Simple regex-based lexer.

Taken from Eli Bendersky [1].

[1] http://stackoverflow.com/questions/133886
"""

import re

class Token(object):
    """ A simple Token structure.
        Contains the token type, value and position.

    TODO(holtgrew): Put file_name before lineno/column.
    """
    def __init__(self, type, val, pos, lineno, column, file_name='<file>'):
        self.type = type
        self.val = val
        self.pos = pos
        self.lineno = lineno
        self.column = column
        self.file_name = file_name

    def __str__(self):
        return '%s(%s) at %d:%d (%d)' % (self.type, repr(self.val), self.lineno, self.column, self.pos)


EOF = Token('EOF', '<eof>', 0, 0, 0)


class LexerError(Exception):
    """ Lexer error exception.

        pos:
            Position in the input line where the error occurred.
    """
    def __init__(self, pos):
        self.pos = pos


class Lexer(object):
    """ A simple regex-based lexer/tokenizer.

        See below for an example of usage.
    """
    def __init__(self, rules, skip_whitespace=True, line_offset=0, col_offset=0):
        """ Create a lexer.

            rules:
                A list of rules. Each rule is a `regex, type`
                pair, where `regex` is the regular expression used
                to recognize the token and `type` is the type
                of the token to return when it's recognized.

            skip_whitespace:
                If True, whitespace (\s+) will be skipped and not
                reported by the lexer. Otherwise, you have to 
                specify your rules for whitespace, or it will be
                flagged as an error.

            line_offset, col_offset:
                Offset for computation.  Used when extracting documentation
                from C/C++ commments and stripping the slashes or the stars.
        """
        self.rules = []

        for name, regex in rules:
            self.rules.append((re.compile(regex), name))

        self.skip_whitespace = skip_whitespace
        self.line_offset = line_offset
        self.col_offset = col_offset
        
        self.re_ws_skip = re.compile('\S')
        
    def input(self, buf, file_name='<mem>', line=0, col=0, offset_col=None):
        """ Initialize the lexer with a buffer as input.
        """
        self.buf = buf
        self.column = col
        self.pos = 0
        self.lineno = line
        self.file_name = file_name
        if offset_col is not None:
            self.col_offset = offset_col

    def token(self):
        """ Return the next token (a Token object) found in the 
            input buffer. None is returned if the end of the 
            buffer was reached. 
            In case of a lexing error (the current chunk of the
            buffer matches no rule), a LexerError is raised with
            the position of the error.
        """
        #import pdb
        #pdb.set_trace()
        if self.pos >= len(self.buf):
            return None
        else:
            if self.skip_whitespace:
                m = self.re_ws_skip.search(self.buf[self.pos:])

                if m:
                    self.lineno += self.buf[self.pos:self.pos + m.start()].count('\n')
                    ln_pos = self.buf.rfind('\n', 0, self.pos + m.start())
                    if ln_pos < 0:
                        ln_pos = -1
                    self.column = self.pos - ln_pos - 1
                    self.pos += m.start()
                else:
                    return None

            for token_regex, token_type in self.rules:
                m = token_regex.match(self.buf[self.pos:])

                if m:
                    value = self.buf[self.pos + m.start():self.pos + m.end()]
                    self.lineno += value.count('\n')
                    ln_pos = self.buf.rfind('\n', 0, self.pos + m.end())
                    if ln_pos < 0:
                        ln_pos = -1
                    self.column = self.pos - ln_pos - 1
                    tok = Token(token_type, value, self.pos, self.lineno + self.line_offset,
                                self.column + self.col_offset, file_name=self.file_name)
                    self.pos += m.end()
                    return tok

            # if we're here, no rule matched
            raise LexerError(self.pos)

    def tokens(self):
        """ Returns an iterator to the tokens found in the buffer.
        """
        while 1:
            tok = self.token()
            if tok is None:
                yield EOF
                break
            yield tok
