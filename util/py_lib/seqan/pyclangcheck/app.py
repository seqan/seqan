#!/usr/bin/env python2
"""pyclangcheck driver code

This code is the driver code for the pyclangcheck tool.

Copyright: (c) 2010, Knut Reinert, FU Berlin
License:   3-clause BSD (see LICENSE)
"""

from __future__ import with_statement

__author__ = 'Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>'

import datetime
import optparse
import os
import os.path
import sys

import clang.cindex as ci

import simple_checks
import violations
import rules

def _hasFileLocation(node):
    """Return True if node has a file lcoation."""
    if hasattr(node, '_has_file_location'):
        return node._has_file_location
    if not hasattr(node, 'location'):
        node._has_file_location = False
        return False
    if not hasattr(node.location, 'file'):
        node._has_file_location = False
        return False
    if not node.location.file:
        node._has_file_location = False
        return False
    if not hasattr(node.location.file, 'name'):
        node._has_file_location = False
        return False
    if not node.location.file.name:
        node._has_file_location = False
        return False
    node._has_file_location = True
    return True


class FileCache(object):
    def __init__(self):
        self.cache = {}

    def get(self, path):
        if self.cache.has_key(path):
            return self.cache[path]
        with open(path, 'rb') as f:
            fcontents = f.readlines()
        self.cache[path] = fcontents
        return self.cache[path]


class CollectViolationsVisitor(object):
    """Visitor for AST nodes that collects rule violations."""
    
    def __init__(self, options, rules):
        self.options = options
        self.rules = rules
        for rule in self.rules:
            rule.visitor = self
        self.stack = []
        self.violations = {}
        self.file_cache = FileCache()
        self.class_stack = []
        self.seen_files = set()
        self.blocked_files = set()
   
    def seenToBlocked(self):
        """Move seen files to blocked files."""
        self.blocked_files |= self.seen_files

    def enterNode(self, node):
        """Called when a node is entered ("pre-order" traversal)."""
        self.stack.append(node)
        ck = ci.CursorKind
        if node.kind in [ck.CLASS_TEMPLATE, ck.CLASS_TEMPLATE_PARTIAL_SPECIALIZATION, ck.CLASS_DECL]:
            ## print 'PUSH CLASS', node.spelling
            self.class_stack.append(node)

        # Mark file as seen for nodes that are directly below the compilation unit.
        if len(self.stack) <= 2 and _hasFileLocation(node):
            self.seen_files.add(node.location.file.name)
        
        if self.options.verbosity >= 2:
            if node.extent.start.file:
                filename = node.extent.start.file.name
                lines = self.file_cache.get(filename)
                start = "%s:%d:%d" % (os.path.basename(filename), node.extent.start.line-1, node.extent.start.column-1)
                end = "%s:%d:%d" % ('#', node.extent.end.line-1, node.extent.end.column-1)
                lines = [x for x in lines[node.extent.start.line-1:node.extent.end.line]]
                if len(lines) == 1:
                    lines[0] = lines[0][node.extent.start.column - 1:node.extent.end.column-1]
                else:
                    lines[0] = lines[0][node.extent.start.column - 1:]
                    lines[-1] = lines[-1][:node.extent.end.column-1]
                if len(lines) > 1:
                    txt = '<multiline>'
                else:
                    txt = ''.join(lines).replace('\n', '\\n')
                print '  ' * len(self.stack), 'Entering', node.kind, node._kind_id, node.spelling, 'txt="%s"' % txt, "%s-%s" % (start, end)
        violations = []
        for rule in self.rules:
            if rule.allowVisit(node):
                #print ' ', ' ' * len(self.stack), 'Checking rule', rule.rule_id
                vs = rule.check(node)
                ## if self.options.verbosity >= 2:
                ##     for v in vs:
                ##         print 'VIOLATION', v
                violations += vs
        for v in violations:
            ## if self.options.verbosity >= 2:
            ##     print v
            self.violations[v.key()] = v

    def exitNode(self, node):
        """Called when a node is left ("post-order" traversa)."""
        self.stack.pop()
        if self.class_stack and self.class_stack[-1] is node:
            ## print 'POP CLASS', node.spelling
            self.class_stack.pop()

    def getCurrentClassName(self):
        """Returns name of current class."""
        if not self.class_stack:
            ## print 'CURRENT CLASS', None
            return None
        ## print 'CURRENT CLASS', self.class_stack[-1].spelling
        return self.class_stack[-1].spelling


class VisitAllowedRule(object):
    """Decides whether a AST node and its children is visited."""
    
    def __init__(self, options, blocked_files):
        self.options = options
        self.include_dirs = [os.path.abspath(x) for x in options.include_dirs]
        self.cache = {}
        self.blocked_files = blocked_files

    def visitAllowed(self, node):
        """Return True if visiting is allowed."""
        # Visit if translation unit.
        if node.kind == ci.CursorKind.TRANSLATION_UNIT:
            return True
        # Don't visit if it has no location (built-in).
        if not _hasFileLocation(node):
            return False
        # Try to hit cache.
        if self.cache.has_key(node.location.file.name):
            return self.cache[node.location.file.name]
        # Check whether the file is blocked.
        if node.location.file.name in self.blocked_files:
            # print 'Blocked', node.location.file.name
            self.cache[node.location.file.name] = False
            return False
        # Check whether node's location is below the include directories.  It is
        # only visited if this is the case.
        filename = os.path.abspath(node.location.file.name)
        result = False
        for x in self.include_dirs:
            if filename.startswith(x):
                # print filename, x
                result = True
                break
        self.cache[node.location.file.name] = result  # Save in cache.
        return result


class AstTraverser(object):
    """Traverses AST tree and applies given visitor object."""
    
    def __init__(self, node_visitor, options):
        self.node_visitor = node_visitor
        self.options = options
        self.visit_allowed_rule = VisitAllowedRule(options, node_visitor.blocked_files)

    def _recurse(self, node):
        """Recursion helper."""
        if not self.visit_allowed_rule.visitAllowed(node):
            return False  # We did not visit this node.
        self.node_visitor.enterNode(node)
        for c in node.get_children():
            self._recurse(c)
        self.node_visitor.exitNode(node)
        return True

    def run(self, filename):
        """Main entry point."""
        index = ci.Index.create()
        args = ['-I%s' % s for s in self.options.include_dirs]
        # print args
        tu = index.parse(filename, args=args)
        if self.options.verbosity >= 1:
            print 'Translation unit: %s.' % tu.spelling
        return self._recurse(tu.cursor)
    
    @classmethod
    def visitFile(klass, filename, node_visitor, options):
        """Don't instantiate AstTraverser yourself, use this function."""
        if options.verbosity >= 1:
            print >>sys.stderr, 'Checking', filename
        traverser = AstTraverser(node_visitor, options)
        res = traverser.run(filename)
        return res != True


def main():
    # ========================================================================
    # Parse command line arguments.
    # ========================================================================
    parser = optparse.OptionParser("USAGE: %prog [options] file.cpp")
    parser.add_option('-s', '--source-file', dest='source_files', default=[],
                      type='string', help='Specify source (.cpp) files.',
                      action='append')
    parser.add_option('-S', '--source-file-file', dest='source_file_files', default=[],
                      type='string', help='File with path to source files.',
                      action='append')
    parser.add_option('-i', '--include-dir', dest='include_dirs', default=[],
                      type='string', help='Specify include directories',
                      action='append')
    parser.add_option('-e', '--exclude-dir', dest='exclude_dirs', default=[],
                      type='string', help='Violations in these directories are not shown.',
                      action='append')
    parser.add_option('-q', '--quiet', dest='verbosity', default=1,
                      action='store_const', const=0, help='Fewer message.')
    parser.add_option('-v', '--verbose', dest='verbosity', default=1,
                      action='store_const', const=2, help='More messages.')
    parser.add_option('--ignore-nolint', dest='ignore_nolint', default=False,
                      action='store_const', const=True, help='Ignore // nolint statements.')
    parser.add_option('--dont-show-source', dest='show_source', default=True,
                      action='store_const', const=False, help='Suppress source line display')
    options, args = parser.parse_args()

    if len(args) != 0:
        parser.error('Incorrect number of arguments!')
        return 1

    # Load source files given in file of paths.
    for filename in options.source_file_files:
        with open(filename, 'rb') as f:
            options.source_files += [x.strip() for x in f.readlines()]

    # ========================================================================
    # Setup traversal.
    # ========================================================================

    # Recursion Rule: Only check symbols within the include directories.
    recurse_rules = []
    recurse_rules.append(rules.InIncludeDirsRule(options.include_dirs, options.exclude_dirs, options.source_files))
    # Define symbol naming rules.
    R = rules.GenericSymbolNameRule
    r = rules
    ck = ci.CursorKind
    check_rules = [
        R(ck.STRUCT_DECL                          , r.RE_STRUCT       , r.RULE_NAMING_STRUCT                ),
        R(ck.UNION_DECL                           , r.RE_TYPE         , r.RULE_NAMING_UNION                 ),
        R(ck.CLASS_DECL                           , r.RE_TYPE         , r.RULE_NAMING_CLASS                 ),
        R(ck.ENUM_DECL                            , r.RE_TYPE         , r.RULE_NAMING_ENUM                  ),
        R(ck.FIELD_DECL                           , r.RE_VARIABLE     , r.RULE_NAMING_FIELD                 ),
        R(ck.ENUM_CONSTANT_DECL                   , r.RE_CONSTANT     , r.RULE_NAMING_ENUM_CONSTANT         ),
        R(ck.FUNCTION_DECL                        , r.RE_FUNCTION     , r.RULE_NAMING_FUNCTION              ),
        R(ck.PARM_DECL                            , r.RE_VARIABLE     , r.RULE_NAMING_PARAMETER             ),
        R(ck.TYPEDEF_DECL                         , r.RE_TYPE         , r.RULE_NAMING_TYPEDEF               ),
        R(ck.CXX_METHOD                           , r.RE_FUNCTION     , r.RULE_NAMING_CXX_METHOD            ),
        R(ck.TEMPLATE_TYPE_PARAMETER              , r.RE_TYPE         , r.RULE_NAMING_TPL_TYPE_PARAMETER    ),
        R(ck.TEMPLATE_NON_TYPE_PARAMETER          , r.RE_CONSTANT     , r.RULE_NAMING_TPL_NON_TYPE_PARAMETER),
        R(ck.TEMPLATE_TEMPLATE_PARAMTER           , r.RE_TYPE         , r.RULE_NAMING_TPL_TPL_PARAMETER     ),
        #R(ck.FUNCTION_TEMPLATE                    , r.RE_FUNCTION     , r.RULE_NAMING_FUNCTION_TPL          ),
        R(ck.CLASS_TEMPLATE                       , r.RE_TYPE_TEMPLATE, r.RULE_NAMING_CLASS_TPL             ),
        R(ck.CLASS_TEMPLATE_PARTIAL_SPECIALIZATION, r.RE_TYPE_TEMPLATE, r.RULE_NAMING_CLASS_TPL_SPEC        ),
        rules.FunctionTemplateRule(),
        rules.VariableNameRule(),
    ]

    # ========================================================================
    # Perform traversal.
    # ========================================================================

    node_visitor = CollectViolationsVisitor(options, check_rules)
    for filename in options.source_files:
        start = datetime.datetime.now()
        res = AstTraverser.visitFile(filename, node_visitor, options)
        node_visitor.seenToBlocked()
        elapsed = datetime.datetime.now() - start
        print >>sys.stderr, '  took', elapsed.seconds, 's'
        if res:
            break

    # ========================================================================
    # Dumber checks (e.g. whitespace at end of file).
    # ========================================================================

    checkers = [simple_checks.WhitespaceChecker(),
                simple_checks.CommentChecker()]
    vs = {}
    for filename in node_visitor.seen_files:
        for checker in checkers:
            vs.update(checker.check(filename))

    # ========================================================================
    # Print violations.
    # ========================================================================

    print 'VIOLATIONS'
    vs.update(node_visitor.violations)
    printer = violations.ViolationPrinter(options.ignore_nolint,
                                          options.show_source)
    printer.show(vs)
    return len(vs) > 0


if __name__ == '__main__':
    sys.exit(main())
