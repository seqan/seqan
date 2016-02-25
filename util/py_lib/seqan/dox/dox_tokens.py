#!/usr/bin/env python2
"""Tokens for Doxygen-like tokens.
"""

import re

# Regular expressions used in LEXER_TOKENS.
IDENTIFIER = r'[\-_a-zA-Z][\-_a-zA-Z0-9]*'
HTML_ATTRIBUTE = IDENTIFIER + r'="[^"]*"'
HTML_TAG = r'<(?:/)?\s*' + IDENTIFIER + '(?:\s+' + HTML_ATTRIBUTE + ')*\s*(?:/)?>'

LITERALS = r'.,;<>+-/*="!:'

def escapeLiterals(s):
    """Escape string of literals s."""
    return r'\|'.join([re.escape(x) for x in s])

# Shortcuts to token sets.
WHITESPACE = set(['SPACE', 'EMPTYLINE', 'BREAK'])
ITEM_STARTING = set(['COMMAND_CLASS', 'COMMAND_CONCEPT', 'COMMAND_FUNCTION', 'COMMAND_MACRO',
                     'COMMAND_METAFUNCTION', 'COMMAND_PAGE', 'COMMAND_MAINPAGE', 'COMMAND_DEFGROUP',
                     'COMMAND_VARIABLE', 'COMMAND_VALUE', 'COMMAND_TAG', 'COMMAND_ENUM',
                     'COMMAND_TYPEDEF', 'COMMAND_ADAPTION'])
CLAUSE_STARTING = set(['COMMAND_SIGNATURE', 'COMMAND_CODE', 'COMMAND_SEE', 'COMMAND_BRIEF',
                       'COMMAND_RETURN', 'COMMAND_THROW', 'COMMAND_PARAM', 'COMMAND_TPARAM',
                       'COMMAND_SECTION', 'COMMAND_SUBSECTION', 'COMMAND_INCLUDE',
                       'COMMAND_EXTENDS', 'COMMAND_IMPLEMENTS', 'COMMAND_SNIPPET',
                       'COMMAND_HEADERFILE', 'COMMAND_DEPRECATED', 'COMMAND_NOTE', 'COMMAND_WARNING',
                       'COMMAND_INTERNAL', 'COMMAND_AKA', 'COMMAND_HTMLONLY', 'COMMAND_DATARACE'])
NON_PARAGRAPH = set(['COMMAND_ENDCODE', 'COMMAND_ENDHTMONLY'])
LINE_BREAKS = set(['BREAK', 'EMPTYLINE'])

# The lexer tokens.
LEXER_TOKENS = (
    # Lexer tokens for commands starting new documentation items.
    ('COMMAND_CLASS',        r'@class'),
    ('COMMAND_CONCEPT',      r'@concept'),
    ('COMMAND_FUNCTION',     r'@fn'),
    ('COMMAND_MACRO',        r'@macro'),
    ('COMMAND_METAFUNCTION', r'@mfn'),
    ('COMMAND_VARIABLE',     r'@var'),
    ('COMMAND_VALUE',        r'@val'),
    ('COMMAND_TAG',          r'@tag'),
    ('COMMAND_ENUM',         r'@enum'),
    ('COMMAND_PAGE',         r'@page'),
    ('COMMAND_MAINPAGE',     r'@mainpage'),
    ('COMMAND_DEFGROUP',     r'@defgroup'),
    ('COMMAND_TYPEDEF',      r'@typedef'),
    ('COMMAND_ADAPTION',     r'@adaption'),
    
    # Lexer tokens for commands starting different modes.
    ('COMMAND_SIGNATURE',    r'@signature'),
    ('COMMAND_CODE',         r'@code'),
    ('COMMAND_ENDCODE',      r'@endcode'),
    
    # Lexer tokens for commands starting clauses.
    ('COMMAND_HEADERFILE',   r'@headerfile'),
    ('COMMAND_EXTENDS',      r'@extends'),
    ('COMMAND_IMPLEMENTS',   r'@implements'),
    ('COMMAND_BRIEF',        r'@brief'),
    ('COMMAND_SEE',          r'@see'),
    ('COMMAND_RETURN',       r'@return'),
    ('COMMAND_THROW',        r'@throw'),
    ('COMMAND_DATARACE',     r'@datarace'),
    ('COMMAND_PARAM',        r'@param'),
    ('COMMAND_TPARAM',       r'@tparam'),
    ('COMMAND_SECTION',      r'@section'),
    ('COMMAND_SUBSECTION',   r'@subsection'),
    ('COMMAND_INCLUDE',      r'@include'),
    ('COMMAND_SNIPPET',      r'@snippet'),
    ('COMMAND_DEPRECATED',   r'@deprecated'),
    ('COMMAND_NOTE',         r'@note'),
    ('COMMAND_WARNING',      r'@warning'),
    ('COMMAND_INTERNAL',     r'@internal'),
    ('COMMAND_AKA',          r'@aka'),
    
    # Lexer tokens for commands for inline use.
    ('COMMAND_LINK',         r'@link'),
    ('COMMAND_ENDLINK',      r'@endlink'),

    # Lexer tokens for HTML-only blocks.
    ('COMMAND_HTMLONLY',     r'@htmlonly'),
    ('COMMAND_ENDHTMLONLY',  r'@endhtmlonly'),
    
    # Space.
    ('SPACE',                r'[ \t]+'),
    ('EMPTYLINE',            r'\n\s*\n\s*'),
    ('BREAK',                r'\n'),

    # In/out tokens for the @param clause.
    ('PARAM_IN_OUT',         r'\[in\]|\[out\]|\[in,out\]'),

    # Lexing words, identifiers etc.
    ('IDENTIFIER',           IDENTIFIER),
    ('NUMBER',               r'[0-9]+'),
    ('WORD',                 r'\w+'),
    ('HTML_TAG',             HTML_TAG),
    ('PUNCTUATION',          r'\S'),
    
    # Lexing punctuation.
    ('HASH',                 r'\#'),
    ('NS_SEP',               r'::'),
    
    # Literals.
    ('LITERAL',              escapeLiterals(LITERALS)),
    )

# Translation of token name to tag name.
LEXER_TOKENS_DICT = dict(LEXER_TOKENS)

def transToken(name):
    """Translate the token name to the token pattern."""
    return LEXER_TOKENS_DICT[name]
