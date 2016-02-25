#!/usr/bin/env python2
"""SeqAn code generation from templates / skeletons.

This module contains code to help the creation of modules, tests, apps etc.
It can be called directly or imported and the main() function can be called.

It will perform the following replacements:

  %(AUTHOR)s  will be replaced by the author's name, either given on command
              line or taken from environment variable SEQAN_AUTHOR.

  %(NAME)s    will be replaced by the name of the generated code.
  %(TITLE)s   will be replaced by the name of the generated, but centered in
              74 characters, to be used in the file header comment.

  %(YEAR)d    will be replaced by the current year.
  %(DATE)s    will be replaced by the current date.
  %(TIME)s    will be replaced by the current time.

  %(HEADER_GUARD)s  will be replaced by the UPPER_CASE_PATH_H_ to the file.

  %(CMAKE_PROJECT_NAME)s  will be replaced by lower_case_path to the target
                          directory.

  %(CMAKE_PROJECT_PATH)s  will be replaced by the path to the target directory.

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
import string

import paths

# Add os.path.relpath if it is not already there, so we can use Python 2.5, too.
# TODO(holtgrew): This could go into a "compatibility" module.
if not 'relpath' in dir(os.path):
    import posixpath
    from posixpath import curdir, sep, pardir, join

    def relpath(path, start=curdir):
        """Return a relative version of a path"""
        if not path:
            raise ValueError("no path specified")
        start_list = posixpath.abspath(start).split(sep)
        path_list = posixpath.abspath(path).split(sep)
        # Work out how much of the filepath is shared by start and path.
        i = len(posixpath.commonprefix([start_list, path_list]))
        rel_list = [pardir] * (len(start_list)-i) + path_list[i:]
        if not rel_list:
            return curdir
        return join(*rel_list)
    os.path.relpath = relpath

# Length of the header comment.
HEADER_CENTER_WIDTH = 74

# Fallback for author string if neither given on command line or environment
# Variable SEQAN_AUTHOR.
DEFAULT_AUTHOR = 'Your Name <your.email@example.net>'

# Program usage string for command line parser.
USAGE = """
Usage: %prog [options] repository NAME
       %prog [options] [module|test|app|demo|header|lheader] NAME LOCATION
       %prog [options] app_tests LOCATION
""".strip()

# Program description, used for command line parser.  Will be wrapped by, though.
DESCRIPTION = """
The SeqAn code generator.

The first version ("repository") is to be be called to create your new entries
below the directory sandbox.  The second version is to be called to create new
library modules, tests, apps, app tests, and demos inside a sandbox.
""".strip()
#"""
#Example:
#
#  %prog repository sandbox/john_doe
#
#The second version is to be called to create new library modules, tests, apps,
#and demos inside a sandbox.  Example:
#
#  %prog module my_module sandbox/john_doe
#
#This command creates a new library module in sandbox/john_doe/include/seqan.
#It consists of the directory my_module, the files my_module.h and
#my_module/my_module_base.h.
#
#  %prog test my_module sandbox/john_doe
#
#This command creates the tests for module "my_module" in sandbox/john_doe.
#
#  %prog app my_app sandbox/john_doe
#
#This command creates a new application named my_app in sandbox/john_doe/apps.
#
#  %prog demo my_demo sandbox/john_doe
#
#This command creates a new demo in sandbox/john_doe/demos.
#""".strip()

def createDirectory(path, dry_run=False):
    print 'mkdir(%s)' % path
    print
    if not dry_run:
        if not os.path.exists(path):
            os.mkdir(path)

def configureFile(target_file, source_file, replacements, dry_run, options):
    print 'Configuring file.'
    print '  Source:', source_file
    print '  Target:', target_file
    print
    if os.path.exists(target_file) and not options.force:
        msg = 'Target file already exists.  Move it away and call the script again.'
        print >>sys.stderr, msg
        return 1

    with open(source_file, 'rb') as f:
        contents = f.read()
    target_contents = contents % replacements
    if dry_run:
        print 'The contents of the target file are:'
        print '-' * 78
        print target_contents
        print '-' * 78
    else:
        with open(target_file, 'wb') as f:
            f.write(target_contents)
    return 0

def _pathToIdentifier(relative_path):
    result = relative_path.replace('/', '_')
    result = result.replace('\\', '_')
    result = result.replace('-', '_')
    result = result.replace('.', '_')
    result = result.replace(' ', '_')
    return result

def buildReplacements(type_, name, location, target_file, options):
    result = {}
    result['AUTHOR'] = options.author
    result['YEAR'] = datetime.date.today().year
    result['TIME'] = datetime.datetime.now().strftime('%H:%M')
    result['DATE'] = datetime.date.today().strftime('%Y-%m-%d')
    result['NAME'] = name
    result['TITLE'] = name.center(HEADER_CENTER_WIDTH).rstrip()
    path = os.path.relpath(target_file, paths.repositoryRoot())
    guard = _pathToIdentifier(path).upper()
    result['HEADER_GUARD'] = guard + '_'
    path = os.path.relpath(os.path.dirname(target_file),
                           paths.repositoryRoot())
    cmake_project_name = _pathToIdentifier(path)
    result['CMAKE_PROJECT_NAME'] = cmake_project_name
    result['CMAKE_PROJECT_PATH'] = path.replace('\\', '\\\\')
    if type_ == 'repository':
        result['REPOSITORY_PSEUDO_TARGET_NAME'] = name.replace('/', '_').replace('\\', '_').replace(' ', '_')
    if type_ == 'app_tests':
        result['APP_NAME'] = os.path.split(os.path.split(location)[0])[1]
        result['APP_NAME_U'] = result['APP_NAME'].upper()
        result['LOCATION'] = os.path.join(os.path.split(os.path.normpath(location))[0])
    return result

def _checkTargetPaths(target_path, options):
    """Check that the path does not exist but its parent does."""
    # Check that the given path does not exist yet.
    if os.path.exists(target_path) and not options.force:
        msg = 'The path %s already exists. Move it and call this script again.'
        print >>sys.stderr, msg % target_path
        return False
    # Check that the parent path already exists.
    if not os.path.exists(os.path.dirname(target_path)):
        msg = 'The parent of the target path does not exist yet: %s'
        print >>sys.stderr, msg % os.path.dirname(target_path)
        print >>sys.stderr, 'Please create it and call this script again.'
        return False
    return True

def createModule(name, location, options):
    include_path = paths.pathToInclude(location)
    seqan_path = os.path.join(include_path, 'seqan')
    module_path = os.path.join(seqan_path, name)
    header_path = os.path.join(seqan_path, '%s.h' % name)
    print 'Creating module in %s' % module_path
    if options.create_dirs and not _checkTargetPaths(module_path, options):
        return 1
    if options.create_dirs and not _checkTargetPaths(header_path, options):
        return 1
    print '  Module path is: %s' % module_path
    print '  Module header path is: %s' % header_path
    print ''
    if options.create_dirs:
        # Create directory.
        createDirectory(module_path, options.dry_run)
    if options.create_programs:
        # Copy over module header.
        source_file = paths.pathToTemplate('module_template', 'module.h')
        target_file = header_path
        replacements = buildReplacements('module', name, seqan_path, target_file, options)
        res = configureFile(target_file, source_file, replacements, options.dry_run, options)
        if res: return res
        # Copy over header inside module.
        source_file = paths.pathToTemplate('module_template', 'header.h')
        target_file = os.path.join(module_path, '%s_base.h' % name)
        replacements = buildReplacements('module', name, seqan_path, target_file, options)
        res = configureFile(target_file, source_file, replacements, options.dry_run, options)
        if res: return res
    return 0

def createTest(name, location, options):
    target_path = paths.pathToTest(location, name)
    print 'Creating test in %s' % target_path
    if options.create_dirs and not _checkTargetPaths(target_path, options):
        return 1
    print '  Target path is: %s' % target_path
    print ''
    if options.create_dirs:
        # Create directory.
        createDirectory(target_path, options.dry_run)
    if options.create_programs:
        # Copy over .cpp file for test and perform replacements.
        source_file = paths.pathToTemplate('test_template', 'test.cpp')
        target_file = os.path.join(target_path, 'test_%s.cpp' % name)
        replacements = buildReplacements('test', name, location, target_file, options)
        res = configureFile(target_file, source_file, replacements, options.dry_run, options)
        if res: return res
        # Copy over .h file for test and perform replacements.
        source_file = paths.pathToTemplate('test_template', 'test.h')
        target_file = os.path.join(target_path, 'test_%s.h' % name)
        replacements = buildReplacements('test', name, location, target_file, options)
        res = configureFile(target_file, source_file, replacements, options.dry_run, options)
        if res: return res
    if options.create_cmakelists:
        # Copy over CMakeLists.txt file for test and perform replacements.
        source_file = paths.pathToTemplate('test_template', 'CMakeLists.txt')
        target_file = os.path.join(target_path, 'CMakeLists.txt')
        replacements = buildReplacements('test', name, location, target_file, options)
        res = configureFile(target_file, source_file, replacements, options.dry_run, options)
        if res: return res
    return 0

def createApp(name, location, options):
    target_path = paths.pathToApp(location, name)
    print 'Creating app in %s' % target_path
    if options.create_dirs and not _checkTargetPaths(target_path, options):
        return 1
    print '  Target path is: %s' % target_path
    print ''
    if options.create_programs:
        # Create directory.
        createDirectory(target_path, options.dry_run)
        # Copy over .cpp file for app and perform replacements.
        source_file = paths.pathToTemplate('app_template', 'app.cpp')
        target_file = os.path.join(target_path, '%s.cpp' % name)
        replacements = buildReplacements('app', name, location, target_file, options)
        res = configureFile(target_file, source_file, replacements, options.dry_run, options)
        if res: return res
    if options.create_cmakelists:
        # Copy over CMakeLists.txt file for app and perform replacements.
        source_file = paths.pathToTemplate('app_template', 'CMakeLists.txt')
        target_file = os.path.join(target_path, 'CMakeLists.txt')
        replacements = buildReplacements('app', name, location, target_file, options)
        res = configureFile(target_file, source_file, replacements, options.dry_run, options)
        if res: return res
    return 0

def createDemo(name, location, options):
    target_path = paths.pathToDemo(location, name)
    print 'Creating demo in %s' % target_path
    if options.create_dirs and not _checkTargetPaths(target_path, options):
        return 1
    print '  Target path is: %s' % target_path
    print ''
    if options.create_programs:
        # Copy over .cpp file for app and perform replacements.
        source_file = paths.pathToTemplate('demo_template', 'demo.cpp')
        target_file = os.path.join(target_path)
        replacements = buildReplacements('demo', name, location, target_file, options)
        res = configureFile(target_file, source_file, replacements, options.dry_run, options)
        if res: return res
    return 0

def createHeader(name, location, options):
    target_path = paths.pathToHeader(location, name)
    print 'Creating (non-library) header in %s' % target_path
    if not _checkTargetPaths(target_path, options):
        return 1
    print '  Target path is: %s' % target_path
    print ''
    # Copy over .h file for app and perform replacements.
    source_file = paths.pathToTemplate('header_template', 'header.h')
    target_file = os.path.join(target_path)
    replacements = buildReplacements('header', name, location, target_file, options)
    res = configureFile(target_file, source_file, replacements, options.dry_run, options)
    if res: return res
    print 'NOTE: Do not forget to add the header to the CMakeLists.txt file!'
    return 0

def createLibraryHeader(name, location, options):
    target_path = paths.pathToHeader(location, name)
    print 'Creating library header in %s' % target_path
    if not _checkTargetPaths(target_path, options):
        return 1
    print '  Target path is: %s' % target_path
    print ''
    # Copy over .h file for app and perform replacements.
    source_file = paths.pathToTemplate('header_template', 'library_header.h')
    target_file = os.path.join(target_path)
    replacements = buildReplacements('library_header', name, location, target_file, options)
    res = configureFile(target_file, source_file, replacements, options.dry_run, options)
    if res: return res
    return 0

def createRepository(location, options):
    print 'Creating module %s' % location
    target_path = paths.pathToRepository(location)
    if options.create_dirs and not _checkTargetPaths(target_path, options):
        return 1
    print '  Target path is: %s' % target_path
    print ''
    if options.create_dirs:
        # Create directories.
        createDirectory(target_path, options.dry_run)
        createDirectory(os.path.join(target_path, 'apps'), options.dry_run)
        createDirectory(os.path.join(target_path, 'demos'), options.dry_run)
        createDirectory(os.path.join(target_path, 'include'), options.dry_run)
        createDirectory(os.path.join(target_path, 'include', 'seqan'), options.dry_run)
        createDirectory(os.path.join(target_path, 'tests'), options.dry_run)
    if options.create_cmakelists:
        # Copy over file ${REPOSITORY}/CMakeLists.txt.
        target_file = os.path.join(target_path, 'CMakeLists.txt')
        source_file = paths.pathToTemplate('repository_template', 'CMakeLists.txt')
        replacements = buildReplacements('repository', location, target_path, target_file, options)
        configureFile(target_file, source_file, replacements, options.dry_run, options)
        # Copy over file ${REPOSITORY}/apps/CMakeLists.txt.
        target_file = os.path.join(target_path, 'apps', 'CMakeLists.txt')
        source_file = paths.pathToTemplate('repository_template', 'apps_CMakeLists.txt')
        replacements = buildReplacements('repository', location, target_path, target_file, options)
        configureFile(target_file, source_file, replacements, options.dry_run, options)
        # Copy over file ${REPOSITORY}/tests/CMakeLists.txt.
        target_file = os.path.join(target_path, 'tests', 'CMakeLists.txt')
        source_file = paths.pathToTemplate('repository_template', 'tests_CMakeLists.txt')
        replacements = buildReplacements('repository', location, target_path, target_file, options)
        configureFile(target_file, source_file, replacements, options.dry_run, options)
        # Copy over file ${REPOSITORY}/demos/CMakeLists.txt.
        target_file = os.path.join(target_path, 'demos', 'CMakeLists.txt')
        source_file = paths.pathToTemplate('repository_template', 'demos_CMakeLists.txt')
        replacements = buildReplacements('repository', location, target_path, target_file, options)
        configureFile(target_file, source_file, replacements, options.dry_run, options)
    return 0

def createAppTests(location, options):
    print 'Creating app tests at %s' % location
    tests_location = os.path.join(location, 'tests')
    target_path = paths.pathToRepository(tests_location)
    if options.create_dirs and not _checkTargetPaths(target_path, options):
        return 1
    print '  Target path is: %s' % target_path
    print ''

    # Create directories.
    if options.create_dirs:
        createDirectory(target_path, options.dry_run)

    # Copy over file ${APP}/tests/generate_outputs.sh
    target_file = os.path.join(target_path, 'generate_outputs.sh')
    source_file = paths.pathToTemplate('app_tests_template', 'generate_outputs.sh')
    replacements = buildReplacements('app_tests', location, target_path, target_file, options)
    configureFile(target_file, source_file, replacements, options.dry_run, options)
    # Copy over file ${APP}/tests/run_tests.py
    target_file = os.path.join(target_path, 'run_tests.py')
    source_file = paths.pathToTemplate('app_tests_template', 'run_tests.py')
    replacements = buildReplacements('app_tests', location, target_path, target_file, options)
    configureFile(target_file, source_file, replacements, options.dry_run, options)

    print '=' * 80
    print 'Do not forget to add the tests in %s:' % os.path.join(location, 'CMakeLists.txt')
    print ''
    print '# Add app tests if Python interpreter could be found.'
    print 'if(PYTHONINTERP_FOUND)'
    print '  add_test(NAME app_test_%s COMMAND ${PYTHON_EXECUTABLE}' % os.path.split(location)[-1]
    print '    ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_tests.py ${CMAKE_SOURCE_DIR}'
    print '    ${CMAKE_BINARY_DIR})'
    print 'endif(PYTHONINTERP_FOUND)'
    print '=' * 80
    
    return 0

def main():
    # Parse arguments.
    parser = optparse.OptionParser(usage=USAGE, description=DESCRIPTION)
    parser.add_option('-s', '--skel-root', dest='skel_root',
                      help=('Set path to the directory where the skeletons '
                            'live in.  Taken from environment variable '
                            'SEQAN_SKELS if available.'),
                      default=os.environ.get('SEQAN_SKELS',
                                             paths.pathToSkeletons()))
    parser.add_option('-a', '--author', dest='author',
                      help=('Set author to use.  Should have the format USER '
                            '<EMAIL>.  Taken from environment variable '
                            'SEQAN_AUTHOR if it exists.'),
                      default=os.environ.get('SEQAN_AUTHOR', DEFAULT_AUTHOR))
    parser.add_option('-d', '--dry-run', dest='dry_run', action='store_true',
                      help='Do not change anything, just simulate.',
                      default=False)
    parser.add_option('-c', '--cmakelists-only', dest='cmakelists_only',
                      action='store_true',
                      help='Only create CMakeLists.txt files',
                      default=False)
    parser.add_option('--force', dest='force', action='store_true',
                      help='Overwrite existing files and directories.',
                      default=False)
    options, args = parser.parse_args()
    options.create_cmakelists = True
    options.create_infos = True
    options.create_dirs = True
    options.create_programs = True
    if options.cmakelists_only:
        options.create_dirs = False
        options.create_programs = False

    if not args:
        parser.print_help(file=sys.stderr)
        return 1
    if len(args) < 2:
        print >>sys.stderr, 'Invalid argument count!'
        return 1
    if args[0] not in ['module', 'test', 'app', 'demo', 'repository',
                       'header', 'lheader', 'app_tests']:
        print >>sys.stderr, 'Invalid template "%s".' % args[0]
        return 1
    if args[0] in['repository', 'app_tests']:
        if len(args) != 2:
            print >>sys.stderr, 'Invalid argument count!'
            return 1

    if args[0] == 'repository':
        return createRepository(args[1], options)
    elif args[0] == 'app_tests':
        return createAppTests(args[1], options)
    elif len(args) != 3:
        print >>sys.stderr, 'Invalid argument count!'
        return 1
    create_methods = {
        'module' : createModule,
        'test': createTest,
        'app': createApp,
        'demo': createDemo,
        'header': createHeader,
        'lheader': createLibraryHeader,
        }
    return create_methods[args[0]](args[1], args[2], options)

if __name__ == '__main__':
   sys.exit(main())

