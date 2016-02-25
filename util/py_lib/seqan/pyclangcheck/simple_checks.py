#!/usr/bin/env python2
"""Simple source code checks, e.g. trailing whitespace."""

from __future__ import with_statement

import bisect
import re
import sys

import violations

RULE_TRAILING_WHITESPACE = 'whitespace.trailing'
RULE_TEXT_TRAILING_WHITESPACE = 'Trailing whitespace is not allowed.'
 
RULE_TODO_ONE_SPACE = 'whitespace.todo'
RULE_TEXT_TODO_ONE_SPACE= 'There should be exactly one space before TODO.'

RULE_TODO_USERNAME = 'readability.todo'
RULE_TEXT_TODO_USERNAME = 'TODO comments should look like this: "// TODO(username): Text".'

RULE_TODO_SPACE_AFTER = 'whitespace.todo'
RULE_TEXT_TODO_SPACE_AFTER = '"TODO (username):" should be followed by a space.'

RE_TODO = r'^//(\s*)TODO(\(.+?\))?:?(\s|$)?'

class WhitespaceChecker(object):
    """Performs simple white space checks."""

    # TODO(holtgrew): Do not allow tabs.

    def check(self, filename):
        vs = []
        with open(filename, 'rb') as f:
            line_no = 0
            for line in f:
                line_no += 1
                if line.rstrip() == line.rstrip('\r\n'):
                    continue
                v = violations.SimpleRuleViolation(
                    RULE_TRAILING_WHITESPACE, filename, line_no,
                    len(line.rstrip()) + 1, RULE_TEXT_TRAILING_WHITESPACE)
                vs.append(v)
        return dict([(v.key(), v) for v in vs])


class SourceFile(object):
    def __init__(self, name):
        self.name = name


class SourceLocation(object):
    def __init__(self, filename, line, column, offset):
        self.file = SourceFile(filename)
        self.line = line
        self.column = column
        self.offset = offset

    def __str__(self):
        return '%s:%d/%d (@%d)' % (self.file.name, self.line, self.column,
                                   self.offset)

    def __repr__(self):
        return str(self)


def enumerateComments(filename):
    # Read file.
    with open (filename, 'rb') as f:
        lines = f.readlines()
        fcontents = ''.join(lines)
    # Build line break index.
    line_starts = [0]
    for line in lines:
        line_starts.append(line_starts[-1] + len(line))
    #print line_starts
    # Search for all comments.
    pattern = re.compile(
        r'//.*?$|/\*.*?\*/|\'(?:\\.|[^\\\'])*\'|"(?:\\.|[^\\"])*"',
        re.DOTALL | re.MULTILINE)
    for match in re.finditer(pattern, fcontents):
        line_start = bisect.bisect(line_starts, match.start(0))
        line_end = bisect.bisect(line_starts, match.end(0) - 1)
        column_start = match.start(0) - line_starts[line_start - 1]
        column_end = match.end(0) - line_starts[line_end - 1]
        yield (SourceLocation(filename, line_start, column_start + 1, match.start(0)),
               SourceLocation(filename, line_end, column_end + 1, match.end(0)),
               match.group(0))


class CommentChecker(object):
    """Performs the checks on comments."""

    def check(self, filename):
        vs = []
        for cstart, cend, comment in enumerateComments(filename):
            if comment.startswith('//'):
                # Check TODO comments.
                match = re.match(RE_TODO, comment)
                if match:
                    if len(match.group(1)) > 1:
                        print comment
                        v = violations.SimpleRuleViolation(
                            RULE_TODO_ONE_SPACE, filename, cstart.line,
                            cstart.column, RULE_TEXT_TODO_ONE_SPACE)
                        vs.append(v)
                    if not match.group(2):
                        v = violations.SimpleRuleViolation(
                            RULE_TODO_USERNAME, filename, cstart.line,
                            cstart.column, RULE_TEXT_TODO_USERNAME)
                        vs.append(v)
                    if match.group(3) != ' ' and match.group(3) != '':
                        v = violations.SimpleRuleViolation(
                            RULE_TODO_SPACE_AFTER, filename, cstart.line,
                            cstart.column, RULE_TEXT_TODO_SPACE_AFTER)
                        vs.append(v)
        return dict([(v.key(), v) for v in vs])

