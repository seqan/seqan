#!/usr/bin/env python
"""
Automatic building of SeqAn apps and releases.
"""

import subprocess
import optparse
import os.path
import sys
import shutil
import tempfile

# The SVN command to use.
SVN_BINARY='svn'
# The CMake command to use.
CMAKE_BINARY='cmake'

# The default value for the SVN tags.
DEFAULT_TAGS_URL='http://svn.seqan.de/seqan/tags'
# The default minimal revision of that tags must have.
DEFAULT_START_REVISION=13708
# The path to the package repository.
DEFAULT_PACKAGE_DB='.'

class MinisculeSvnWrapper(object):
    """Minimal SVN wrapper."""

    def ls(self, url):
        """Execute 'svn ls ${url}'."""
        print >>sys.stderr, 'Executing "%s %s %s"' % (SVN_BINARY, 'ls -v', url)
        popen = subprocess.Popen([SVN_BINARY, 'ls', '-v', url],
                                 stdout=subprocess.PIPE)
        out_data, err_data = popen.communicate()
        if popen.returncode != 0:
            print >>sys.stderr, 'ERROR during SVN call.'
            return 1
        lines = out_data.splitlines()
        revs_tags = [(int(line.split()[0]), line.split()[-1]) for line in lines]
        res = []
        for rev, tag in revs_tags:
            if tag == './':
                continue  # Skip dot-slash.
            tag2 = tag[:-1]
            res.append((rev, tag2))
        return res

    def co(self, url, dest_dir):
        """Execute 'svn co ${url} ${dest_dir}'."""
        print >>sys.stderr, 'Executing "%s %s %s %s"' % (SVN_BINARY, 'co', url, dest_dir)
        popen = subprocess.Popen([SVN_BINARY, 'co', url, dest_dir])
        popen.wait()
        return popen.returncode


class Package(object):
    """Represent a package with a given name, version, OS and architeture."""

    def __init__(self, name, version, os, word_size, pkg_format):
        self.name = name
        self.version = version
        self.os = os
        self.word_size = word_size
        SYS_NAMES = {'Windows': {'32': 'win32', '64': 'win64'},
                     'Linux': {'32': 'Linux-i686', '64': 'Linux-x86_64'},
                     'Mac': {'32': 'Mac-i686', '64': 'Mac-x86_64'}}
        self.system_name = SYS_NAMES[os][word_size]
        self.pkg_format = pkg_format

    def fileName(self):
        return '%s-%s-%s.%s' % (self.name, self.version, self.system_name,
                                self.pkg_format)


class BuildStep(object):
    """Management of one build step."""

    def __init__(self, path, name, version, os, word_size, pkg_formats,
                 svn_url, make_args, options, tmp_dir=None):
        self.base_path = path
        self.name = name
        self.version = version  # TODO(holtgrew): Unused, overwritten below.
        self.major_version = int(version.split('.')[0])
        self.minor_version = int(version.split('.')[1])
        self.patch_version = 0
        self.version = '%d.%d.%d' % (self.major_version, self.minor_version, self.patch_version)
        if len(version.split('.')) > 2:
            self.patch_version = int(version.split('.')[2])
        self.os = os
        self.word_size = word_size
        self.pkg_formats = pkg_formats
        if name == 'seqan':
            self.packages = [Package(name + suffix, self.version, os, word_size, f)
                             for f in pkg_formats for suffix in ['-apps', '-library']]
        else:
            self.packages = [Package(name, self.version, os, word_size, f)
                             for f in pkg_formats]
        self.svn_url = svn_url
        self.make_args = make_args
        self.options = options
        # If set then this is used instead of a random name in TMPDIR.
        self.tmp_dir = tmp_dir

    def buildNeeded(self):
        """Returns whether one of the package files is missing."""
        for p in self.packages:
            package_path = os.path.join(self.base_path, p.fileName())
            if not os.path.exists(package_path):
                return True
            elif self.options.verbosity >= 1:
                print >>sys.stderr, 'File %s exists.' % package_path
        return False

    def copyArchives(self, build_dir):
        """Copy build packages to base_path directory."""
        for p in self.packages:
            from_ = os.path.join(build_dir, p.fileName())
            if os.path.exists(from_):
                to = os.path.join(self.base_path, os.path.basename(from_))
                print >>sys.stderr, "Copying %s => %s" % (from_, to)
                shutil.copyfile(from_, to)
            else:
                print >>sys.stderr, "%s does not exist (not fatal)" % from_

    def buildSeqAnRelease(self, checkout_dir, build_dir):
        """Build SeqAn release: Apps and library build."""
        # Build seqan-apps.
        #
        # Create build directory.
        print >>sys.stderr, "Creating build directory %s" % (build_dir,)
        os.mkdir(build_dir)
        # Execute CMake.
        cmake_args = [CMAKE_BINARY,
                      "-DSEQAN_BUILD_PACKAGE=SEQAN_RELEASE_APPS"]
        if self.os in ['Linux', 'Mac']:
            cmake_args += ['-m%s' % self.word_size]
        print >>sys.stderr, 'Executing CMake: "%s"' % (' '.join(cmake_args),)
        # Execute Make.
        make_args = [MAKE_BINARY, "package"] + self.make_args
        print >>sys.stderr, 'Executing Make: "%s"' % (' '.join(make_args),)
        # Copy over the archives.
        self.copyArchives(build_dir)
        # Remove build directory.
        print >>sys.stderr, 'Removing build directory %s' % build_dir
        shutil.rmtree(build_dir)
        # Build seqan-library.
        #
        # Create build directory.
        print >>sys.stderr, "Creating build directory %s" % (build_dir,)
        os.mkdir(build_dir)
        # Execute CMake.
        cmake_args = [CMAKE_BINARY,
                      "-DSEQAN_BUILD_SYSTEM=SEQAN_RELEASE_LIBRARY"]
        print >>sys.stderr, 'Executing CMake: "%s"' % (' '.join(cmake_args),)
        # Execute Make.
        make_args = [MAKE_BINARY, "package"] + self.make_args
        print >>sys.stderr, 'Executing Make: "%s"' % (' '.join(make_args),)
        self.copyArchives(build_dir)
        # Remove build directory.
        print >>sys.stderr, 'Removing build directory %s' % build_dir
        shutil.rmtree(build_dir)

    def buildApp(self, checkout_dir, build_dir):
        """Build an application."""
        # Create build directory.
        print >>sys.stderr, "Creating build directory %s" % (build_dir,)
        os.mkdir(build_dir)
        # Execute CMake.
        cmake_args = [CMAKE_BINARY, checkout_dir,# '-G', 'Visual Studio 10',
                      "-DCMAKE_BUILD_TYPE=Release",
                      "-DSEQAN_BUILD_SYSTEM=APP:%s" % self.name,
                      "-DSEQAN_APP_VERSION=%d.%d.%d" %
                        (self.major_version, self.minor_version, self.patch_version)]
        # Use appropriate CMake flags for OS and processor.
        if self.word_size == '32':
            cmake_args.append('-DSEQAN_SYSTEM_PROCESSOR=i686')
            if self.os != 'Windows':
                cmake_args.append('-DCMAKE_CXX_FLAGS=-m32')
            else:
                cmake_args += ['-G', 'Visual Studio 10']
        else:  # self.word_size == '64'
            if self.os == 'Windows':
                cmake_args += ['-G', 'Visual Studio 10 Win64']
        print >>sys.stderr, 'Executing CMake: "%s"' % (' '.join(cmake_args),)
        #for key in sorted(os.environ.keys()):
        #    print key, ': ', os.environ[key]
        popen = subprocess.Popen(cmake_args, cwd=build_dir, env=os.environ.copy())
        out_data, err_data = popen.communicate()
        if popen.returncode != 0:
            print >>sys.stderr, 'ERROR during make call.'
            print out_data
            print err_data
            return 1
        # Build and package project.
        make_args = [CMAKE_BINARY, '--build', build_dir, '--target', 'package']
        if self.options.verbosity > 1:
            make_args.insert(1, 'VERBOSE=1')
        print >>sys.stderr, 'Building with CMake: "%s"' % (' '.join(make_args),)
        popen = subprocess.Popen(make_args, cwd=build_dir)
        out_data, err_data = popen.communicate()
        if popen.returncode != 0:
            print >>sys.stderr, 'ERROR during make call.'
            print out_data
            print err_data
            return 1
        # Copy out archives.
        self.copyArchives(build_dir)
        # Remove build directory.
        print >>sys.stderr, 'Removing build directory %s' % build_dir
        shutil.rmtree(build_dir)

    def tmpDir(self):
      print 'self.tmp_dir = %s' % self.tmp_dir
      if self.tmp_dir:
        if not os.path.exists(self.tmp_dir):
          os.makedirs(self.tmp_dir)
        return self.tmp_dir
      else:
        return tempfile.mkdtemp()

    def execute(self):
        """Execute build step."""
        # Create temporary directory.
        tmp_dir = self.tmpDir()
        print >>sys.stderr, 'Temporary directory is %s' % (tmp_dir,)
        # Create SVN checkout in temporary directory.
        checkout_dir = os.path.join(tmp_dir, os.path.basename(self.svn_url))
        print >>sys.stderr, 'Creating checkout in %s' % checkout_dir
        svn = MinisculeSvnWrapper()
        svn.co(self.svn_url, checkout_dir)
        # Create build directory.
        suffix = '-build-%s-%s' % (self.os, self.word_size)
        build_dir = os.path.join(tmp_dir, os.path.basename(self.svn_url) + suffix)
        if os.path.exists(build_dir):
          print >>sys.stderr, 'Removing build directory %s' % (build_dir,)
          shutil.rmtree(build_dir)
        # Perform the build.  We have to separate between app and whole SeqAn releases.
        if self.name == 'seqan':
            self.buildSeqAnRelease(checkout_dir, build_dir)
        else:
            self.buildApp(checkout_dir, build_dir)
        # Remove temporary directory again.
        if not self.tmp_dir:  # Only remove it not set with --tmp-dir.
          print >>sys.stderr, 'Removing temporary directory %s' % (tmp_dir,)
          shutil.rmtree(tmp_dir)


def work(options):
    """Run the individual steps."""
    # Get the revisions and tag names.
    svn = MinisculeSvnWrapper()
    revs_tags = [(rev, tag) for (rev, tag) in svn.ls(options.tags_url)
                 if rev >= options.start_revision and
                    '-' in tag]
    # Enumerate all package names that we could enumerate.
    print 'revs_tags = %s' % revs_tags
    print 'word_sizes = %s' % options.word_sizes
    for rev, tag in revs_tags:
        name, version = tag.rsplit('-', 1)
        for word_size in options.word_sizes.split(','):
            # Create build step for this package name.
            pkg_formats = options.package_formats.split(',')
            svn_url = options.tags_url + '/' + tag
            build_step = BuildStep(options.package_db, name, version, options.os,
                                   word_size, pkg_formats, svn_url,
                                   options.make_args.split(), options, options.tmp_dir)
            # Check whether we need to build this.
            if not build_step.buildNeeded():
                continue  # Skip
            # Execute build step.
            build_step.execute()
    return 0


def main():
    """Program entry point."""
    # Parse Arguments.
    parser = optparse.OptionParser()
    parser.add_option('-t', '--tags-url', dest='tags_url',
                      default=DEFAULT_TAGS_URL,
                      help='This URL is searched for tags.', metavar='URL')
    parser.add_option('--package-db', dest='package_db', type='string',
                      default=DEFAULT_PACKAGE_DB,
                      help='Path the directory with the packages.')
    parser.add_option('-s', '--start-revision', dest='start_revision',
                      default=DEFAULT_START_REVISION,
                      type='int', help='Ignore all tags with smaller revision.')
    parser.add_option('-v', dest='verbosity', action='count', default=1,
                      help='Increase verbosity.')
    parser.add_option('--package-formats', dest='package_formats',
                      default='tar.bz2,zip',
                      help='Expect the following packages to be created.')
    parser.add_option('--os', dest='os', help='Expect the packages to be created for this OS.',
                      default='Linux')
    parser.add_option('--word-sizes', dest='word_sizes', default='32,64',
                      help='Build binaries with the given word sizes')
    parser.add_option('--make-args', dest='make_args', type="string", default='',
                      help='Arguments for make.')
    parser.add_option('--tmp-dir', dest='tmp_dir', type='string', default=None,
                      help='Temporary directory to use. Use this to reuse the same checkout.')
    parser.epilog = ('The program will use the environment variable TMPDIR as '
                     'the directory for temporary files.')

    options, args = parser.parse_args()
    if args:
        parser.error('No arguments expected!')
        return 1

    # Fire up work.
    print >>sys.stderr, 'Running SeqAn Auto Builder'
    return work(options)
