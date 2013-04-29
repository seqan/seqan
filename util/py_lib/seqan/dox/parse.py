#!/usr/bin/env python

# ===========================================================================
# Tokenizer
# ===========================================================================

import sys

import ply.lex

tokens = (
    'COMMAND_CLASS',
    'COMMAND_CONCEPT',
    'COMMAND_FUNCTION',
    'COMMAND_METAFUNCTION',
    
    'COMMAND_SIGNATURE',
    'COMMAND_CODE',
    'COMMAND_ENDCODE',
    
    'COMMAND_SEE',
    'COMMAND_RETURN',
    'COMMAND_PARAM',
    'COMMAND_TPARAM',
    'COMMAND_SECTION',
    'COMMAND_SUBSECTION',
    'COMMAND_INCLUDE',
    
    'SPACE',
    'BREAK',
    'EMPTYLINE',
    
    'IDENTIFIER',
    'WORD',
    'NUMBER',
    'HTML_TAG',

    'HASH',
    'NS_SEP',
    'OTHER_NONSPACE',
    )

t_ignore = '\r'

literals = '.,;<>+-/*="!:'

states = (
    ('signature', 'exclusive'),
    ('code', 'exclusive'),
)

# ---------------------------------------------------------------------------
# Commands starting documentation items.
# ---------------------------------------------------------------------------

def t_COMMAND_CONCEPT(t):
    r'@concept'
    return t

def t_COMMAND_CLASS(t):
    r'@class'
    return t

def t_COMMAND_FUNCTION(t):
    r'@fn'
    return t

def t_COMMAND_METAFUNCTION(t):
    r'@mfn'
    return t

def t_COMMAND_PAGE(t):
    r'@page'
    return t

# ---------------------------------------------------------------------------
# Commands starting different modes
# ---------------------------------------------------------------------------

def t_COMMAND_SIGNATURE(t):
    r'@signature'
    t.lexer.begin('signature')
    return t

def t_COMMAND_CODE(t):
    r'@code'
    t.lexer.begin('code')
    return t

def t_COMMAND_ENDCODE(t):
    r'@endcode'
    return t

# ---------------------------------------------------------------------------
# Commands for clauses in documentation items.
# ---------------------------------------------------------------------------

def t_COMMAND_SEE(t):
    r'@see'
    return t

def t_COMMAND_RETURN(t):
    r'@return\|@returns'
    return t

def t_COMMAND_PARAM(t):
    r'@param'
    return t

def t_COMMAND_TPARAM(t):
    r'@tparam'
    return t

def t_COMMAND_SECTION(t):
    r'@section'
    return t

def t_COMMAND_SUBSECTION(t):
    r'@subsection'
    return t

def t_COMMAND_INCLUDE(t):
    r'@include'
    return t

# ---------------------------------------------------------------------------
# Commands for inline use.
# ---------------------------------------------------------------------------

def t_COMMAND_LINK(t):
    r'@link'
    return t

def t_COMMAND_ENDLINK(t):
    r'@endlink'
    return t

# ---------------------------------------------------------------------------
# Space
# ---------------------------------------------------------------------------

SPACE = r'[ \t]+'
EMPTYLINE = r'\n\s*\n\s*'
BREAK = r'\n'

@ply.lex.TOKEN(SPACE)
def t_SPACE(t):
    return t

# Multiple line breaks, at one or more empty lines.
@ply.lex.TOKEN(EMPTYLINE)
def t_EMPTYLINE(t):
    t.lexer.lineno += t.value.count('\n')
    return t

# One line break.
@ply.lex.TOKEN(BREAK)
def t_BREAK(t):
    t.lexer.lineno += 1
    return t

# ---------------------------------------------------------------------------
# Lexing words, identifiers etc.
# ---------------------------------------------------------------------------

IDENTIFIER = r'[_a-zA-Z][_a-zA-Z0-9]*'
HTML_ATTRIBUTE = IDENTIFIER + r'="[^"]*"'
HTML_TAG = r'<(?:/)?\s*' + IDENTIFIER + '(?:\s+' + HTML_ATTRIBUTE + ')*\s*(?:/)?>'

@ply.lex.TOKEN(IDENTIFIER)
def t_IDENTIFIER(t):
    return t

def t_NUMBER(t):
    r'[0-9]+'
    return t

def t_WORD(t):
    r'\S+'
    return t

@ply.lex.TOKEN(HTML_TAG)
def t_HTML_TAG(t):
    return t

# ---------------------------------------------------------------------------
# Lexing punctuation.
# ---------------------------------------------------------------------------

def t_HASH(t):
    r'\#'
    return t

def t_NS_SEP(t):
    r'::'
    return t


# ---------------------------------------------------------------------------
# Lexer Error Handling
# ---------------------------------------------------------------------------

# Compute column. 
#     input is the input text string
#     token is a token instance
def findColumn(input, token):
    last_cr = input.rfind('\n',0,token.lexpos)
    if last_cr < 0:
	    last_cr = 0
    column = (token.lexpos - last_cr) + 1
    return column

def t_error(t):
    vars = (t.value[0], t.lexer.lineno, findColumn(t.lexer.lexdata, t))
    print >>sys.stderr, 'Illegal characters "%s" at line %d column %d.' % vars
    t.lexer.skip(1)

# ---------------------------------------------------------------------------
# State "signature"
# ---------------------------------------------------------------------------

t_signature_ignore = '\r'

t_signature_error = t_error

@ply.lex.TOKEN(EMPTYLINE)
def t_signature_EMPTYLINE(t):
    t.lexer.begin('INITIAL')
    t.lexer.lineno += t.value.count('\n')
    return t

@ply.lex.TOKEN(r'[ \t]+')
def t_signature_SPACE(t):
    return t

def t_signature_COMMAND_SEE(t):
    r'@see'
    t.lexer.begin('INITIAL')
    return t_COMMAND_SEE(t)

def t_signature_COMMAND_RETURN(t):
    r'@return\|@returns'
    t.lexer.begin('INITIAL')
    return t_COMMAND_RETURN(t)

def t_signature_COMMAND_PARAM(t):
    r'@param'
    t.lexer.begin('INITIAL')
    return t_COMMAND_PARAM(t)

def t_signature_COMMAND_TPARAM(t):
    r'@tparam'
    t.lexer.begin('INITIAL')
    return t_COMMAND_TPARAM(t)

def t_signature_COMMAND_SECTION(t):
    r'@section'
    t.lexer.begin('INITIAL')
    return t_COMMAND_SECTION(t)

def t_signature_COMMAND_SUBSECTION(t):
    r'@subsection'
    t.lexer.begin('INITIAL')
    return t_COMMAND_SUBSECTION(t)

def t_signature_COMMAND_INCLUDE(t):
    r'@include'
    t.lexer.begin('INITIAL')
    return t_COMMAND_INCLUDE(t)

# Continue counting lines in signature state.
@ply.lex.TOKEN(BREAK)
def t_signature_BREAK(t):
    t.lexer.lineno += 1
    return t

@ply.lex.TOKEN(r'\S+')
def t_signature_WORD(t):
    return t

# ---------------------------------------------------------------------------
# State "code"
# ---------------------------------------------------------------------------

# In the state code, we parse everything except "@endcode" as word or space.
# We leave the state when reaching the @encode token

t_code_ignore = ''

t_code_error = t_error

def t_code_COMMAND_ENDCODE(t):
    r'@endcode'
    t.lexer.begin('INITIAL')
    return t

def t_code_WORD(t):
    r'\S+'
    return t

def t_code_SPACE(t):
    r'\s+'
    t.lexer.lineno += t.value.count('\n')
    return t

# ===========================================================================
# Parser
# ===========================================================================

INPUT = ''

def p_documentation(p):
    """ documentation : class_doc
                      | empty
    """
    print 'documentation_doc'
    p[0] = ('documentation', p[1])

#def p_concept_doc(p):
#    """ concept_doc : empty
#    """
#    print 'concept_doc'
#    p[0] = ('concept_doc')

def p_class_doc(p):
    """ class_doc : class_def class_body
    """
    print 'class_doc'
    p[0] = ('class_doc', p[1], p[2])

def p_class_def(p):
    """ class_def : COMMAND_CLASS SPACE entity_name stop
    """
    p[0] = ('@class', ''.join(p[3][1]))
    print p[0]

def p_class_body(p):
    """ class_body : signature_clause class_body
                   | tparam_clause class_body
                   | paragraph
                   | empty
    """
    if len(p) > 2:
        p[0] = ('@signature', p[2], p[3])
    print p[0]

def p_tparam_clause(p):
    """ tparam_clause : COMMAND_TPARAM SPACE IDENTIFIER SPACE text EMPTYLINE
                      | COMMAND_TPARAM SPACE IDENTIFIER SPACE text
    """
    p[0] = ('@tparam', p[3], ''.join(p[5][1]))
    print p[0]

def p_signature_clause(p):
    """ signature_clause : COMMAND_SIGNATURE SPACE signature EMPTYLINE
                         | COMMAND_SIGNATURE SPACE signature
    """
    p[0] = ('@signature', ''.join(p[3][1]))
    print p[0]

# TODO(holtgrew): Make the following more compact (p_signature_*).

def p_signature_X_signature(p):
    """ signature : SPACE signature
                  | WORD signature
                  | BREAK signature
    """
    text = [p[1]]
    if len(p[2]) >= 2:
        text += p[2][1]
    p[0] = ('signature', text)

def p_signature_X(p):
    """ signature : WORD
                  | BREAK
                  | SPACE
    """
    p[0] = ('signature', [p[1]])

def p_entity_name_identifier_space_entity_name(p):
    """ entity_name : IDENTIFIER SPACE entity_name
    """
    p[0] = ('entity_name', [p[1], p[2]] + p[3][1])

def p_entity_name_identifier(p):
    """ entity_name : IDENTIFIER
    """
    p[0] = ('entity_name', [p[1]])

def p_text_X_text(p):
    """ text : word_or_identifier text
             | SPACE text
             | BREAK text
    """
    text = [p[1]]
    if p[2] is not None:
        text += p[2][1]
    p[0] = ('text', text)

def p_text_X(p):
    """ text : word_or_identifier
             | BREAK
             | SPACE
    """
    p[0] = ('text', [p[1]])

def p_paragraph(p):
    """ paragraph : text EMPTYLINE
                  | text
    """

def p_word_or_identifier(p):
    """ word_or_identifier : WORD
                           | IDENTIFIER
    """
    p[0] = p[1]

def p_stop(p):
    """ stop : BREAK
             | EMPTYLINE
    """

def p_empty(p):
    """ empty :
    """

def showError(input, line, column):
    print input.splitlines()[line-1]
    print ' ' * (column - 2) + '^'

def p_error(p):
    last_cr = INPUT.rfind('\n', 0, p.lexpos)
    if last_cr < 0:
	    last_cr = 0
    col = (p.lexpos - last_cr) + 1
    print >>sys.stderr, 'Syntax error at line %d column %d: "%s".' % (p.lineno, col, repr(p.value))
    showError(INPUT, p.lineno, col)

# ===========================================================================
# Program
# ===========================================================================

import logging
import optparse
import sys

import ply.lex
import ply.yacc

def doMain(filename, lex_only, debug):
    logging.basicConfig(level=logging.DEBUG)
    logger = logging.getLogger()
    with open(filename, 'r') as f:
        fcontents = f.read()
        global INPUT
        INPUT = fcontents
    if debug:
        print ',-- FILE CONTENTS\n%s`-----------------' % fcontents
    lex = ply.lex.lex(debug=debug, debuglog=logger)
    if lex_only:
        lex.input(fcontents)
        while True:
            tok = lex.token()
            if not tok: break
            print tok
        return 0
    parser = ply.yacc.yacc(debug=debug)
    if debug:
        print parser.parsedebug(fcontents, lexer=lex, debug=logger)
    else:
        print parser.parse(fcontents, lexer=lex)
    return 0

def main():
    parser = optparse.OptionParser()
    parser.add_option('--lex-only', dest='lex_only', help='Lex only.',
                      default=False, action='store_true')
    parser.add_option('--debug', dest='debug', help='Debug.',
                      default=False, action='store_true')
    parser.add_option('-i', dest='input_file', help='Input file.')
    options, args = parser.parse_args()
    if not options.input_file:
        parser.error('Missing input file.')
    return doMain(options.input_file, options.lex_only, options.debug)

if __name__ == '__main__':
    sys.exit(main())

