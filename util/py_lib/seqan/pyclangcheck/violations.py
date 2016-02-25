#!/usr/bin/env python2
"""Code related to violations and suppressions."""

from __future__ import with_statement

import os
import os.path
import sys

import app
import rules


class RuleViolation(object):
    def __init__(self, rule_id, violator, file, line, column):
        self.rule_id = rule_id
        self.violator = violator
        self.file = file
        self.line = line
        self.column = column
    
    def key(self):
        return (self.file, self.line, self.column, self.rule_id, self.violator)
    
    def __str__(self):
        msg = '[%s:%d/%d] %s "%s": %s'
        return msg % ('/'.join(self.file.split('/')[-2:]), self.line, self.column,
                      self.rule_id, self.violator, rules.RULE_TEXTS[self.rule_id])


class SimpleRuleViolation(object):
    def __init__(self, rule_id, file, line, column, msg):
        self.rule_id = rule_id
        self.file = file
        self.line = line
        self.column = column
        self.msg = msg
    
    def key(self):
        return (self.file, self.line, self.column, self.rule_id)
    
    def __str__(self):
        msg = '[%s:%d/%d] %s : %s'
        return msg % ('/'.join(self.file.split('/')[-2:]), self.line, self.column,
                      self.rule_id, self.msg)


class NolintManager(object):
    """Manage the lines ending in '//nolint'."""

    def __init__(self):
        self.locations = {}

    def hasNolint(self, filename, lineno):
        filename = os.path.abspath(filename)
        # Ensure that the nolint lines are registered in self.locations[filename].
        if not self.locations.has_key(filename):
            line_set = set()
            with open(filename, 'rb') as f:
                line_no = 0
                for line in f:
                    line_no += 1
                    if line.strip().endswith('// nolint'):
                        ## print 'nolint', filename, line_no
                        line_set.add(line_no)
            self.locations[filename] = line_set
        # Query self.locations[filename].
        return lineno in self.locations[filename]


class ViolationPrinter(object):
    def __init__(self, ignore_nolint, show_source):
      self.nolints = NolintManager()
      self.file_cache = app.FileCache()
      self.ignore_nolint = ignore_nolint
      self.show_source = show_source

    def show(self, vs):
        previous = None
        for k in sorted(vs.keys()):
            violation = vs[k]
            if self.ignore_nolint or not self.nolints.hasNolint(violation.file, violation.line):
                print violation
                line = self.file_cache.get(violation.file)[violation.line - 1]
                if self.show_source:
                    print line.rstrip()
                    print '%s^' % (' ' * (violation.column - 1))
                    print
            previous = violation
        print 'Total: %d violations' % len(vs)

