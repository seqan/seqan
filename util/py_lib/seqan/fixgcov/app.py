#!/usr/bin/env python2
"""Small libclang based app to fix gcov output.

Fix gcov output with templates.  This is done by first parsing in the .cpp files
(compilation units) with libclang.  The AST is then parsed and all lines within
composite statements ({ stmt; stmt; ... }) are memoized as 'interesting' lines.
The resulting interesting lines are serialized to a location file with pickle.
Finally, gcov output files are read and updated.  If a line is interesting but
not marked as covered or uncovered (marker '-'), it is marked as uncovered
(marker '#####').

USAGE: fixgcov.py -i $include_dir -g $gcov_file
USAGE: fixgcov.py -i $include_dir -s $source_file

Copyright: (c) 2010, Knut Reinert, FU Berlin
License:   3-clause BSD (see LICENSE)
"""

from __future__ import with_statement

__author__ = 'Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>'

import optparse
import os
import pickle
import sys

import clang.cindex as ci


def _hasFileLocation(node):
    """Return True if node has a file lcoation."""
    if not hasattr(node, 'location'):
        return False
    if not hasattr(node.location, 'file'):
        return False
    if not node.location.file:
        return False
    if not hasattr(node.location.file, 'name'):
        return False
    if not node.location.file.name:
        return False
    return True


class CollectCompoundStatementNodeVisitor(object):
    """Visitor for AST nodes that collects compound statements."""
    
    def __init__(self, options):
        self.options = options
        self.stack = []
        self.ranges = []
    
    def enterNode(self, node):
        """Called when a node is entered ("pre-order" traversal)."""
        self.stack.append(node)
        ## print '  ' * len(self.stack), node.kind,
        num_children = len([x for x in node.get_children()])
        ## if _hasFileLocation(node):
        ##     print node.location.file.name, '%d-%d' % (node.extent.start.line, node.extent.end.line)
        ## else:
        ##     print
        # Only add range for statements that are no compound statements.  Add
        # for empty compounds.
        if not node.kind.is_statement():
            ## print 'skipping, no statement'
            return
        if node.kind == ci.CursorKind.COMPOUND_STMT and num_children > 0:
            ## print 'skipping, non-empty compound statement', num_children
            return
        if node.kind == ci.CursorKind.DECL_STMT:
            return  # Skip declarations.
        # Only add if has file location.
        if _hasFileLocation(node):
            self.ranges.append((node.location.file.name, node.extent.start.line, node.extent.end.line))

    def exitNode(self, node):
        """Called when a node is left ("post-order" traversa)."""
        self.stack.pop()


class VisitAllowedRule(object):
    """Decides whether a AST node and its children is visited."""
    
    def __init__(self, options):
        self.options = options
        self.include_dirs = [os.path.abspath(x) for x in options.include_dirs]
        self.cache = {}

    def visitAllowed(self, node):
        """Return True if visiting is allowed."""
        # TODO(holtgrew): For this application, stopping at compound statements
        # would be enough.  Visit if translation unit.
        if node.kind == ci.CursorKind.TRANSLATION_UNIT:
            return True
        # Don't visit if it has no location (built-in).
        if not _hasFileLocation(node):
            return False
        # Try to hit cache.
        if self.cache.has_key(node.location.file.name):
            return self.cache[node.location.file.name]
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
        self.visit_allowed_rule = VisitAllowedRule(options)

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
        traverser = AstTraverser(node_visitor, options)
        res = traverser.run(filename)
        return res == True


def main():
    """Main entry point."""
    # ========================================================================
    # Parse command line arguments.
    # ========================================================================
    parser = optparse.OptionParser("USAGE: %prog [options] -s file.cpp")
    parser.add_option('-I', '--include-dir', dest='include_dirs', default=[],
                      type='string', help='Specify include directories',
                      action='append')
    parser.add_option('-s', '--src-file', dest='source_files', default=[],
                      type='string', help='Specify compilation units.',
                      action='append')
    parser.add_option('--src-file-files', dest='source_files_files', default=[],
                      type='string', help='Specify file with paths to compilation units.',
                      action='append')
    parser.add_option('-l', '--location-file', dest='location_file',
                      default='locations.dat', type='string',
                      help='Path to file with compound statement locations.')
    parser.add_option('-g', '--gcov-file', dest='gcov_files', default=[],
                      type='string', help='Specify gcov files to process.',
                      action='append')
    parser.add_option('--gcov-files', dest='gcov_files_files', default=[],
                      type='string', help='Specify gcov files to process.',
                      action='append')
    parser.add_option('-q', '--quiet', dest='verbosity', default=1,
                      action='store_const', const=0, help='Fewer message.')
    parser.add_option('-v', '--verbose', dest='verbosity', default=1,
                      action='store_const', const=2, help='More messages.')
    options, args = parser.parse_args()
    if len(args) != 0:
        parser.error('Incorrect number of arguments!')
        return 1

    options.include_dirs += [os.path.abspath(os.path.dirname(s)) for s in options.source_files]

    # ========================================================================
    # Read in files with paths from arguments.
    # ========================================================================

    for path in options.source_files_files:
      with open(path, 'rb') as f:
        options.source_files += [x.strip() for x in f.readlines()]
    for path in options.gcov_files_files:
      with open(path, 'rb') as f:
        options.gcov_files += [x.strip() for x in f.readlines()]

    if not options.source_files and not options.gcov_files:
        parser.error('Neither source nor gcov file given!')
        return 1

    # ========================================================================
    # Collect interesting lines if any source files given.
    # ========================================================================

    if options.source_files:
        # If any source file is given, all given source files are parsed and all
        # lines with compound statements in all included files are written to
        # the location file.
        if options.verbosity >= 1:
            print >>sys.stderr, 'Building Locations'
        if options.verbosity >= 2:
            print >>sys.stderr, '=================='

        # Fire off AST traversal.
        if options.verbosity >= 1:
            print >>sys.stderr, 'AST Traversal'
        node_visitor = CollectCompoundStatementNodeVisitor(options)
        for src in options.source_files:
            if options.verbosity >= 2:
                print >>sys.stderr, '  Compilation Unit', src
            AstTraverser.visitFile(src, node_visitor, options)

        # Convert locations into points.
        locations = {}
        for filename, start, stop in node_visitor.ranges:
            filename = os.path.abspath(filename)
            for i in range(start, stop + 1):
                locations.setdefault(filename, set()).add(i)

        # Write out the source locations.
        if options.verbosity >= 1:
            print >>sys.stderr, 'Writing out locations to', options.location_file
        with open(options.location_file, 'wb') as f:
            pickle.dump(locations, f)

    # ========================================================================
    # Process GCOV files if any are given.
    # ========================================================================

    if options.gcov_files:
        # If no source files and gcov files are given then
        if options.verbosity >= 1:
            print >>sys.stderr, 'Updating gcov Results'
        if options.verbosity >= 2:
            print >>sys.stderr, '====================='

        if not options.source_files:
            if options.verbosity >= 1:
                print >>sys.stderr, 'Loading locations from', options.location_file
            with open(options.location_file, 'rb') as f:
                locations = pickle.load(f)

        for filename in options.gcov_files:
            filename = os.path.abspath(filename)
            if options.verbosity >= 2:
                print >>sys.stderr, 'Processing', filename
            with open(filename, 'rb') as f:
                lines = f.readlines()
            pos0 = lines[0].find(':')
            pos1 = lines[0].find(':', pos0 + 1)
            source = None
            result = []
            skip = False
            for i, line in enumerate(lines):
                coverage = line[:pos0]
                lineno = int(line[pos0 + 1:pos1].strip())
                slineno = line[pos0 + 1:pos1]
                txt = line[pos1 + 1:]
                if txt.startswith('Source:'):
                    source = os.path.abspath(txt[len('Source:'):].strip())
                    if not locations.has_key(source):
                        if options.verbosity >= 2:
                            print >>sys.stderr, '  Skipping.'
                        skip = True
                        break
                if not source or lineno == 0:
                    result.append(line)
                    continue  # Proceed only if in file.
                if lineno in locations[source] and coverage.strip() == '-':
                    coverage = ('%%%ds' % pos0) % '#####'
                result.append(':'.join([coverage, slineno, txt]))
            # Write back file if not skipped.
            if skip:
                continue
            with open (filename, 'wb') as f:
              f.write(''.join(result))
            #print ''.join(result)

# Entry point if called as program.
if __name__ == '__main__':
    sys.exit(main())
