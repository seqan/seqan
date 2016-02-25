#!/usr/bin/env python2
"""
Automatic building of SeqAn apps and releases.
"""

from __future__ import print_function

import subprocess
import optparse
import os.path
import re
import sys
import shutil
import tempfile

# The git command to use.
GIT_BINARY='git'
# The CMake command to use.
CMAKE_BINARY='cmake'

# The default repository URL.
REPOSITORY_URL='https://github.com/seqan/seqan.git'
# The path to the package repository.
DEFAULT_PACKAGE_DB='.'

# Regular expression to use for tag names.
TAG_RE=r'.*-v\d+\.\d+\.\d(-\w+)?'

class MinisculeGitWrapper(object):
    """Minimal git wrapper."""

    def lsRemote(self, url):
        """Execute 'git ls-remote ${url} --tags'."""
        # Execute ls-remote command.
        print('Executing "%s %s %s"' % (GIT_BINARY, 'ls-remote --tags', url), file=sys.stderr)
        popen = subprocess.Popen([GIT_BINARY, 'ls-remote', '--tags', url],
                                 stdout=subprocess.PIPE)
        out_data, err_data = popen.communicate()
        print('  => %d' % popen.returncode, file=sys.stderr)
        if popen.returncode != 0:
            print('ERROR during git call.', file=sys.stderr)
            return 1
        # Parse out revisions and tags names.
        lines = out_data.splitlines()
        revs_tags = [(line.split()[0], line.split()[-1]) for line in lines]
        res = []
        for rev, tag in revs_tags:
            if '^{}' in tag:
                continue  # Skip with ^{} in tag name
            tag2 = tag[10:]
            res.append((rev, tag2))
        return res

    def checkout(self, path, treeish):
        """Execute "git checkout" in the checkout at path."""
        # Executing git checkout.
        args = [GIT_BINARY, 'checkout', treeish]
        print('Executing "%s" in "%s"' % (' '.join(args), path), file=sys.stderr)
        popen = subprocess.Popen(args, cwd=path)
        out_data, err_data = popen.communicate()
        if popen.returncode != 0:
            print('ERROR during git call.', file=sys.stderr)
        # Executing force resetting to current revision.
        args = [GIT_BINARY, 'rm', '--cached', '.']
        print('Executing "%s" in "%s"' % (' '.join(args), path), file=sys.stderr)
        popen = subprocess.Popen(args, cwd=path)
        out_data, err_data = popen.communicate()
        if popen.returncode != 0:
            print('ERROR during git call.', file=sys.stderr)
        args = [GIT_BINARY, 'reset', '--hard']
        print('Executing "%s" in "%s"' % (' '.join(args), path), file=sys.stderr)
        popen = subprocess.Popen(args, cwd=path)
        out_data, err_data = popen.communicate()
        if popen.returncode != 0:
            print('ERROR during git call.', file=sys.stderr)
        return popen.returncode

    def archive(self, path, treeish, output, prefix):
        """Execute git archive."""
        args = [GIT_BINARY, 'archive', '--prefix=%s/' % prefix, '--output=%s' % output, treeish]
        print('Executing "%s" in "%s"' % (' '.join(args), path), file=sys.stderr)
        popen = subprocess.Popen(args, cwd=path)
        out_data, err_data = popen.communicate()
        if popen.returncode != 0:
            print('ERROR during git call.', file=sys.stderr)
        return popen.returncode

    def clone(self, url, tag, dest_dir):
        """Execute 'git clone ${url} ${dest_dir}' and then get specific tag."""
        # Clone repository
        args = [GIT_BINARY, 'clone', url, dest_dir]
        print('Executing "%s"' % ' '.join(args), file=sys.stderr)
        popen = subprocess.Popen(args)
        popen.wait()
        print('  => %d' % popen.returncode, file=sys.stderr)
        if popen.returncode != 0:
            return popen.returncode
        return self.checkout(dest_dir, tag)


class Package(object):
    """Represent a package with a given name, version, OS and architeture."""

    def __init__(self, name, version, os, word_size, pkg_format):
        self.name = name
        self.version = version
        self.os = os
        self.word_size = word_size
        SYS_NAMES = {'Windows': {'32': 'win32-i686', '64': 'win64-x86_64'},
                     'Linux': {'32': 'Linux-i686', '64': 'Linux-x86_64'},
                     'Mac': {'32': 'Darwin-i686', '64': 'Darwin-x86_64'}}
        self.system_name = SYS_NAMES[os][word_size]
        self.pkg_format = pkg_format

    def fileName(self):
        if self.name == 'seqan-library':
            return '%s-%s.%s' % (self.name, self.version, self.pkg_format)
        else:
            return '%s-%s-%s.%s' % (self.name, self.version, self.system_name,
                                    self.pkg_format)


class BuildStep(object):
    """Management of one build step."""

    def __init__(self, path, treeish, name, version, os, word_size, pkg_formats,
                 repository_url, make_args, options, tmp_dir=None):
        self.base_path = path
        self.treeish = treeish
        self.name = name
        self.version = version.split('-', 1)[0]
        print ( 'Version: %s Self.Version: %s' % (version, self.version), file=sys.stdout)
        self.major_version = int(self.version.split('.')[0])
        print ( 'Major_Version: %s' % self.major_version, file=sys.stdout)
        self.minor_version = int(self.version.split('.')[1])
        print ( 'Minor_Version: %s' % self.minor_version, file=sys.stdout)
        self.patch_version = 0
        if len(self.version.split('.')) > 2:
            self.patch_version = int(self.version.split('.')[2])
        self.version = '%d.%d.%d' % (self.major_version, self.minor_version, self.patch_version)
        print ( 'Self_Version: %s' % self.version, file=sys.stdout)
        self.os = os
        self.word_size = word_size
        self.pkg_formats = pkg_formats
        if name == 'seqan':
            self.packages = [Package(name + suffix, self.version, os, word_size, f)
                             for f in pkg_formats for suffix in ['-apps', '-library']]
        else:
            self.packages = [Package(name, self.version, os, word_size, f)
                             for f in pkg_formats]
        self.repository_url = repository_url
        self.make_args = make_args
        self.options = options
        # If set then this is used instead of a random name in TMPDIR.
        self.tmp_dir = tmp_dir

    def buildNeeded(self):
        """Returns whether one of the package files is missing."""
        for p in self.packages:
            package_path = os.path.join(self.base_path, p.name, p.fileName())
            if 'x86' in package_path and 'x86_64' not in package_path:  # fix processor name
                package_path = package_path.replace('x86', 'x86_64')
            if 'win32' in package_path or 'win64' in package_path:  # fix OS name
                package_path = package_path.replace('win32', 'Windows').replace('win64', 'Windows')
            if 'Darwin' in package_path:  # fix OS name
                package_path = package_path.replace('Darwin', 'Mac')
            if not os.path.exists(package_path):
                if self.options.verbosity >= 1:
                    print('File %s does not exist yet.' % package_path, file=sys.stderr)
                return True
            elif self.options.verbosity >= 1:
                print('File %s exists.' % package_path, file=sys.stderr)
        return False

    def copyArchives(self, build_dir):
        """Copy built packages to base_path directory."""
        for p in self.packages:
            from_ = os.path.join(build_dir, p.fileName())
            if os.path.exists(from_):
                to = os.path.join(self.base_path, p.name, os.path.basename(from_))
                if not os.path.exists(os.path.dirname(to)):  # Create directory if necessary.
                    os.makedirs(os.path.dirname(to))
                print("Copying %s => %s" % (from_, to), file=sys.stderr)
                if 'x86' in to and 'x86_64' not in to:  # fix processor name
                    to = to.replace('x86', 'x86_64')
                if 'win32' in to or 'win64' in to:  # fix OS name
                    to = to.replace('win32', 'Windows').replace('win64', 'Windows')
                if 'Darwin' in to:  # fix OS name
                    to = to.replace('Darwin', 'Mac')
                shutil.copyfile(from_, to)
            else:
                print('%s does not exist (not fatal)' % from_, file=sys.stderr)

    def buildSeqAnRelease(self, checkout_dir, build_dir):
        """Build SeqAn release: Apps and library build."""
        # Build seqan-apps.
        #
        # Create build directory.
        if not os.path.exists(build_dir):
            print('Creating build directory %s' % (build_dir,), file=sys.stderr)
            os.mkdir(build_dir)
        # Execute CMake.
        cmake_args = [CMAKE_BINARY, checkout_dir,
                      '-DSEQAN_BUILD_SYSTEM=SEQAN_RELEASE_APPS']
        # Use appropriate CMake flags for OS and processor.
        # Use appropriate CMake flags for OS and processor.
        if self.word_size == '32':
            cmake_args.append('-DSEQAN_SYSTEM_PROCESSOR=i686')
            if self.os != 'Windows':
                cmake_args.append('-DCMAKE_CXX_FLAGS=-m32')
            else:
                cmake_args += ['-G', 'Visual Studio 10']
        else:  # self.word_size == '64'
            cmake_args.append('-DSEQAN_SYSTEM_PROCESSOR=x86_64')
            if self.os == 'Windows':
                cmake_args += ['-G', 'Visual Studio 10 Win64']
        print('Executing CMake: "%s"' % (' '.join(cmake_args),), file=sys.stderr)
        popen = subprocess.Popen(cmake_args, cwd=build_dir, env=os.environ.copy())
        out_data, err_data = popen.communicate()
        if popen.returncode != 0:
            print('ERROR during make call.', file=sys.stderr)
            print(out_data, file=sys.stderr)
            print(err_data, file=sys.stderr)
            return 1
        # Execute Make.
        cmake_args = [CMAKE_BINARY, '--build', build_dir, '--target', 'package', '--config', 'Release', '--'] + self.make_args
        print('Building with CMake: "%s"' % (' '.join(cmake_args),), file=sys.stderr)
        popen = subprocess.Popen(cmake_args, cwd=build_dir, env=os.environ.copy())
        out_data, err_data = popen.communicate()
        if popen.returncode != 0:
            print('ERROR during make call.', file=sys.stderr)
            print(out_data, file=sys.stderr)
            print(err_data, file=sys.stderr)
            return 1
        # Copy over the archives.
        self.copyArchives(build_dir)
        # Remove build directory.
        if not self.options.keep_build_dir:
            print('Removing build directory %s' % build_dir, file=sys.stderr)
            shutil.rmtree(build_dir)
        # Build seqan-library.
        #
        # Create build directory.
        if not os.path.exists(build_dir):
            print("Creating build directory %s" % (build_dir,), file=sys.stderr)
            os.mkdir(build_dir)
        # Execute CMake.
        cmake_args = [CMAKE_BINARY, checkout_dir,
                      "-DSEQAN_BUILD_SYSTEM=SEQAN_RELEASE_LIBRARY"]
        print('Executing CMake: "%s"' % (' '.join(cmake_args),), file=sys.stderr)
        popen = subprocess.Popen(cmake_args, cwd=build_dir, env=os.environ.copy())
        out_data, err_data = popen.communicate()
        if popen.returncode != 0:
            print('ERROR during make call.', file=sys.stderr)
            print(out_data, file=sys.stderr)
            print(err_data, file=sys.stderr)
            return 1
        # Build Docs
        cmake_args = [CMAKE_BINARY, '--build', build_dir, '--target', 'docs', '--'] + self.make_args
        print('Building with CMake: "%s"' % (' '.join(cmake_args),), file=sys.stderr)
        popen = subprocess.Popen(cmake_args, cwd=build_dir, env=os.environ.copy())
        out_data, err_data = popen.communicate()
        if popen.returncode != 0:
            print('ERROR during make dox call.', file=sys.stderr)
            print(out_data, file=sys.stderr)
            print(err_data, file=sys.stderr)
        # Execute Make.
        cmake_args = [CMAKE_BINARY, '--build', build_dir, '--target', 'package', '--'] + self.make_args
        print('Building with CMake: "%s"' % (' '.join(cmake_args),), file=sys.stderr)
        popen = subprocess.Popen(cmake_args, cwd=build_dir, env=os.environ.copy())
        out_data, err_data = popen.communicate()
        if popen.returncode != 0:
            print('ERROR during make call.', file=sys.stderr)
            print(out_data, file=sys.stderr)
            print(err_data, file=sys.stderr)
            return 1
        self.copyArchives(build_dir)
        # Remove build directory.
        if not self.options.keep_build_dir:
            print('Removing build directory %s' % build_dir, file=sys.stderr)
            shutil.rmtree(build_dir)

    def buildApp(self, checkout_dir, build_dir):
        """Build an application."""
        # Create build directory.
        print("Creating build directory %s" % (build_dir,), file=sys.stderr)
        if not os.path.exists(build_dir):
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
            cmake_args.append('-DSEQAN_SYSTEM_PROCESSOR=x86_64')
            if self.os == 'Windows':
                cmake_args += ['-G', 'Visual Studio 10 Win64']
        print('Executing CMake: "%s"' % (' '.join(cmake_args),), file=sys.stderr)
        #for key in sorted(os.environ.keys()):
        #    print(key, ': ', os.environ[key], file=sys.stderr)
        popen = subprocess.Popen(cmake_args, cwd=build_dir, env=os.environ.copy())
        out_data, err_data = popen.communicate()
        if popen.returncode != 0:
            print('ERROR during make call.', file=sys.stderr)
            print(out_data, file=sys.stderr)
            print(err_data, file=sys.stderr)
            return 1
        # Build and package project.
        make_args = [CMAKE_BINARY, '--build', build_dir, '--target', 'package', '--config', 'Release']
        if self.options.verbosity > 1:
            make_args.insert(1, 'VERBOSE=1')
        print('Building with CMake: "%s"' % (' '.join(make_args),), file=sys.stderr)
        popen = subprocess.Popen(make_args, cwd=build_dir)
        out_data, err_data = popen.communicate()
        if popen.returncode != 0:
            print('ERROR during make call.', file=sys.stderr)
            print(out_data, file=sys.stderr)
            print(err_data, file=sys.stderr)
            return 1
        # Copy out archives.
        self.copyArchives(build_dir)
        # Remove build directory.
        if not self.options.keep_co_dir:
            print('Removing build directory %s' % build_dir, file=sys.stderr)
            shutil.rmtree(build_dir)

    def tmpDir(self):
      print('self.tmp_dir = %s' % self.tmp_dir, file=sys.stderr)
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
        print('Temporary directory is %s' % (tmp_dir,), file=sys.stderr)
        # Create Git checkout in temporary directory.
        checkout_dir = os.path.join(tmp_dir, os.path.basename(self.repository_url))
        print('Creating checkout in %s' % checkout_dir, file=sys.stderr)
        git = MinisculeGitWrapper()
        git.clone(self.repository_url, self.treeish, checkout_dir)
        # Create build directory.
        suffix = '-build-%s-%s' % (self.os, self.word_size)
        build_dir = os.path.join(tmp_dir, os.path.basename(self.repository_url) + suffix)
        if os.path.exists(build_dir) and not self.options.keep_build_dir:
            print('Removing build directory %s' % (build_dir,), file=sys.stderr)
            shutil.rmtree(build_dir)
        # insert app tags
        subprocess.call(['../tag-apps.sh', checkout_dir])
        # Perform the build.  We have to separate between app and whole SeqAn releases.
        if self.name == 'seqan':
            self.buildSeqAnRelease(checkout_dir, build_dir)
        else:
            self.buildApp(checkout_dir, build_dir)
        if not self.options.keep_co_dir:
            print('Removing checkout directory %s' % (checkout_dir,), file=sys.stderr)
            shutil.rmtree(checkout_dir)
        # Remove temporary directory again.
        if self.tmp_dir and not self.options.keep_tmp_dir:
            # Only remove if not explicitely given and not forced to keep.
            print('Removing temporary directory %s' % (tmp_dir,), file=sys.stderr)
            shutil.rmtree(tmp_dir)


def workTags(options):
    """Run the individual steps for tags."""
    # Get the revisions and tag names.
    git = MinisculeGitWrapper()
    #revs_tags = [(rev, tag) for (rev, tag) in git.lsRemote(options.repository_url)
                 #if re.match(TAG_RE, tag)]
    tags = [tag for (rev, tag) in git.lsRemote(options.repository_url)
           if re.match(TAG_RE, tag)]
    tags.extend(subprocess.check_output(['../tag-apps.sh', os.getcwd(), 'printonly']).split('\n'))
    # Enumerate all package names that we could enumerate.
    print('tags = %s' % tags, file=sys.stderr)
    print('word_sizes = %s' % options.word_sizes, file=sys.stderr)
    for tag in tags:
        name, version = tag.rsplit('-v', 1)
        #version = version[1:]
        print ('Tag: %s Name: %s Version: %s' % (tag, name, version), file=sys.stdout)
        for word_size in options.word_sizes.split(','):
            # Create build step for this package name.
            pkg_formats = options.package_formats.split(',')
            build_step = BuildStep(options.package_db, tag, name, version, options.os,
                                   word_size, pkg_formats, options.repository_url,
                                   options.make_args.split(), options, options.tmp_dir)
            # Check whether we need to build this.
            if not build_step.buildNeeded():
                continue  # Skip
            # Execute build step.
            build_step.execute()
    return 0


def workTrunk(options):
    """Run the individual steps for the trunk with fake tag name."""
    # Get the revisions and tag names.
    git = MinisculeGitWrapper()
    # Enumerate all package names that we could enumerate.
    print('fake tag = %s' % options.build_trunk_as, file=sys.stderr)
    print('word_sizes = %s' % options.word_sizes, file=sys.stderr)
    name, version = options.build_trunk_as.rsplit('-', 1)
    version = version[1:]
    for word_size in options.word_sizes.split(','):
        # Create build step for this package name.
        pkg_formats = options.package_formats.split(',')
        build_step = BuildStep(options.package_db, 'master', name, version, options.os,
                               word_size, pkg_formats, options.repository_url,
                               options.make_args.split(), options, options.tmp_dir)
        # Check whether we need to build this.
        if not build_step.buildNeeded():
            continue  # Skip
        # Execute build step.
        build_step.execute()
    return 0


def workSrcTar(options):
    """Build the source tarball."""
    # Get the revisions and tag names.
    git = MinisculeGitWrapper()
    revs_tags = [(rev, tag) for (rev, tag) in git.lsRemote(options.repository_url)
                 if re.match(TAG_RE, tag)]
    # Enumerate all package names that we could enumerate.
    for rev, tag in revs_tags:
        # Build URL.
        name, version = tag.rsplit('-', 1)
        version = version[1:]  # remove prefix "v"
        if name != 'seqan':
            continue  # only build source tarballs for seqan
        # Create destination file name.
        file_name = '%s-src-%s.tar.gz' % (name, version)
        dest = os.path.join(options.package_db, '%s-src' % name, file_name)
        # Check whether we need to rebuild.
        if os.path.exists(dest):
            print('Skipping %s; already exists.' % dest, file=sys.stderr)
            continue
        # Create temporary directory.
        if options.tmp_dir:
            if not os.path.exists(options.tmp_dir):
                os.makedirs(options.tmp_dir)
            tmp_dir = options.tmp_dir
        else:
            tmp_dir = tempfile.mkdtemp()
        print('Temporary directory is %s' % tmp_dir, file=sys.stderr)
        # Create git checkout in temporary directory.
        checkout_dir = os.path.join(tmp_dir, tag)
        print('Creating checkout in %s' % checkout_dir, file=sys.stderr)
        from_ = os.path.join(tmp_dir, file_name)
        git.clone(options.repository_url, tag, checkout_dir)
        # Create target directory if it does not exist yet.
        if not os.path.exists(os.path.dirname(dest)):  # Create directory if necessary.
            os.makedirs(os.path.dirname(dest))
        # Create tarball.
        git.archive(checkout_dir, tag, dest, prefix='%s-%s' % (name, version))
        # Remove temporary directory again.
        if tmp_dir and not options.keep_tmp_dir:
            # Only remove if not explicitely given and not forced to keep.
            print('Removing temporary directory %s' % (tmp_dir,), file=sys.stderr)
            shutil.rmtree(tmp_dir)
    return 0


def work(options):
    """Run the steps."""
    if options.src_tar:
        return workSrcTar(options)
    elif not options.build_trunk_as:
        return workTags(options)
    else:
        return workTrunk(options)


def main():
    """Program entry point."""
    # Parse Arguments.
    parser = optparse.OptionParser()

    parser.add_option('-u', '--repository-url', default=REPOSITORY_URL,
                      help='The git repository URL.', metavar='URL')
    parser.add_option('--package-db', dest='package_db', type='string',
                      default=DEFAULT_PACKAGE_DB,
                      help='Path the directory with the packages.')

    parser.add_option('--src-tar', dest='src_tar', action='store_true',
                      help='If specified then only the src tarball will be created')
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
    parser.add_option('--build-trunk-as', dest='build_trunk_as', type='string', default=None,
                      help='Build current trunk with this string as a tag name.')
    parser.add_option('--keep-build-dir', dest='keep_build_dir', default=False,
                      action='store_true', help='Keep build directory.')
    parser.add_option('--keep-tmp-dir', dest='keep_tmp_dir', default=False,
                      action='store_true', help='Keep temporary directory.')
    parser.add_option('--keep-co-dir', dest='keep_co_dir', default=False,
                      action='store_true', help='Keep checkout directory.')
    parser.epilog = ('The program will use the environment variable TMPDIR as '
                     'the directory for temporary files.')

    options, args = parser.parse_args()
    if args:
        parser.error('No arguments expected!')
        return 1

    options.package_db = os.path.abspath(options.package_db)

    # Fire up work.
    print('Running SeqAn Auto Builder', file=sys.stderr)
    return work(options)
