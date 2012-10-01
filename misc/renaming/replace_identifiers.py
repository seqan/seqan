#!/usr/bin/env python

import sys
import re

from helpers import *

PROGRAM_USAGE = """
SeqAn script to replace invalid identifiers (previously collected) in the SeqAn
codebase.

USAGE: replace_identifiers.py BASE_PATH [REPLACEMENTS]

BASE_PATH is the root path of all the folders to be searched.
REPLACEMENTS is a file of ``"key:value"`` pairs which contain the invalid
    identifier and the replacement string.
    If this file is not given, it is attemped to read the replacements from the
    standard input stream.
""".strip()

def replace_all(text, subst):
    """
    Perform the substitutions given by the dictionary ``subst`` on ``text``.
    """
    for old in subst.keys():
        text = old.sub(subst[old], text)

    return text


def validate_file(file, subst):
    """
    Perform the substitutions given by the dictionary ``subst`` on ``file``.
    """

    #print file

    code = ''
    try:
        f = open(file, 'r') 
    finally:
        code = f.read()

    old_len = len(code)
    replaced = replace_all(code, subst)
    #assert old_len == len(replaced)    
    
    open(file, 'w').write(replaced)


def build_subst_table(file):
    """
    Read the substitutions defined in ``file`` and build a substitution table.
    """
    table = {}

    for line in file:
        old, new = line.rstrip('\r\n').split(':')
        table[re.compile(r'\b%s\b' % old.strip())] = new.strip()

    return table


def main():
    # Either read from stdin or expect a file path in the second argument.
    # Since there is no reliable way of checking for an attached stdin on
    # Windows, just assume good faith if the file name isn't given.
    use_stdin = len(sys.argv) == 2
    if not (len(sys.argv) == 3 or use_stdin):
        print >>sys.stderr, 'ERROR: Invalid number of arguments.'
        print >>sys.stderr, PROGRAM_USAGE
        return 1

    if use_stdin:
        print >>sys.stderr, "Attempting to read from stdin ..."

    project_path = sys.argv[1]
    replacements_file = sys.stdin if use_stdin else open(sys.argv[2], 'r')

    substitutions = build_subst_table(replacements_file)

    for file in all_files(project_path):
        validate_file(file, substitutions)

    return 0


if __name__ == '__main__':
    sys.exit(main())
