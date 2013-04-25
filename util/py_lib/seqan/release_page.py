#!/usr/bin/env python
"""Build the SeqAn Releases Website."""

import operator
import optparse
import os
import os.path
import re
import sys
import xml.sax.saxutils

import pyratemp

# Patterns matching seqan apps and library.
LIBRARY_PATTERN = (r'seqan-library-([0-9])\.([0-9])(?:\.([0-9]))?\.'
                   '(tar\.gz|tar\.bz2|zip)')
APPS_PATTERN = (r'seqan-apps-([0-9])\.([0-9])(?:\.([0-9]))?-'
                '(Linux|Mac|Windows)-(x86_64|i686)?'
                '\.(tar\.gz|tar\.bz2|zip|exe)')
# The regular expression to use for matching patterns.
PACKAGE_PATTERN = (r'(.*)-([0-9]+)\.([0-9]+)(?:\.([0-9]+))?-'
                   '(Linux|Mac|Windows)-(x86_64|i686)?'
                   '\.(tar\.gz|tar\.bz2|zip|exe)')
# The operating systems that we expect.
OPERATING_SYSTEMS = ['Linux', 'Mac', 'Windows', 'src']
# The architectures that we expect.
ARCHITECTURES = ['x86_64', 'i686', 'src']
# The file formats.
FORMATS = ['tar.gz', 'tar.bz2', 'zip', 'exe']
# Path to template.
TPL_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        'release_page.html')
PACKAGE_TPL_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                'one_package.html')
# Base URL for links.
BASE_URL='http://packages.seqan.de'


class Arch(object):
    def __init__(self, name):
        self.name = name
        self.files = {}


class Packages(object):
    def __init__(self, os_):
        self.os = os_
        self.archs = {}
        for arch in ARCHITECTURES:
            self.archs[arch] = Arch(arch)


class Version(object):
    def __init__(self, version):
        self.version = version
        self.packages = {}
        for os_ in OPERATING_SYSTEMS:
            self.packages[os_] = Packages(os_)


class Software(object):
    def __init__(self, name):
        self.name = name
        self.versions = {}


class PackageDatabase(object):
    def __init__(self, path):
        self.path = path
        self.seqan_apps = Software('SeqAn Apps')
        self.seqan_library = Software('SeqAn Library')
        self.softwares = {}

    def load(self):
        # Two craw directory structure by two levels.
        xs = []
        for x in os.listdir(self.path):
            if os.path.isdir(os.path.join(self.path, x)):
                for y in os.listdir(os.path.join(self.path, x)):
                    xs.append(y)
        for x in xs:
            if re.match(LIBRARY_PATTERN, x):
                major, minor, patch, suffix = re.match(LIBRARY_PATTERN, x).groups()
                if not patch:
                    patch = '0'
                major_minor_patch = '%s.%s.%s' % (major, minor, patch)
                software = self.seqan_library
                if not major_minor_patch in software.versions:
                    software.versions[major_minor_patch] = Version(major_minor_patch)
                version = software.versions[major_minor_patch]
                version.packages['src'].archs['src'].files[suffix] = x
            elif re.match(APPS_PATTERN, x):
                major, minor, patch, os_, arch, suffix = re.match(APPS_PATTERN, x).groups()
                if not patch:
                    patch = '0'
                major_minor_patch = '%s.%s.%s' % (major, minor, patch)
                software = self.seqan_apps
                if not major_minor_patch in software.versions:
                    software.versions[major_minor_patch] = Version(major_minor_patch)
                version = software.versions[major_minor_patch]
                version.packages[os_].archs[arch].files[suffix] = x
            elif re.match(PACKAGE_PATTERN, x):  # individual apps
                filename = x
                name, major, minor, patch, os_, arch, suffix = re.match(PACKAGE_PATTERN, x).groups()
                if not patch:
                    patch = '0'
                major_minor_patch = '%s.%s.%s' % (major, minor, patch)

                if not name in self.softwares:
                    self.softwares[name] = Software(name)
                software = self.softwares[name]

                if not major_minor_patch in software.versions:
                    software.versions[major_minor_patch] = Version(major_minor_patch)
                version = software.versions[major_minor_patch]

                version.packages[os_].archs[arch].files[suffix] = filename
            else:
                pass


class RssItem(object):
    """One RSS item."""
    def __init__(self, title, description, link):
        self.title = title
        self.description = description
        self.link = link

    def generate(self):
        tpl = ('<item>\n'
               '  <title>%s</title>\n'
               '  <summary>%s</summary>\n'
               '  <link>%s</link>\n'
               '</item>\n')
        return tpl % (self.title, self.description, self.link)


class RssFeed(object):
    """Feed with one channel."""

    def __init__(self, title, description, link):
        self.title = title
        self.description = description
        self.link = link
        self.items = []

    def generate(self):
        tpl = ('<?xml version="1.0" encoding="UTF-8" ?>\n'
               '<rss version="2.0">\n'
               '  <title>%s</title>\n'
               '  <description>%s</description>\n'
               '\n'
               '%s'
               '</rss>\n')
        items_s = '\n'.join([i.generate() for i in self.items])
        return tpl % (self.title, self.description, items_s)


class RssWriter(object):
    """Writing of RSS files for a PackageDB."""

    def __init__(self, out_dir, package_db, base_url):
        self.out_dir = out_dir
        self.package_db = package_db
        self.base_url = base_url

    def generate(self):
        """Create output RSS files."""
        for sname, software in self.package_db.softwares.items():
            feed = RssFeed(sname, '', '')
            for vname, version in software.versions.items():
                description = 'Version %s of %s.' % (vname, sname)
                link = '%s/%s#%s' % (self.base_url, sname, vname)
                item = RssItem('%s %s' % (sname, vname), description, link)
                feed.items.append(item)
            path = os.path.join(self.out_dir, sname, 'package.rss')
            print >>sys.stderr, 'Writing %s' % path
            with open(path, 'wb') as f:
                f.write(feed.generate())


def work(options):
    print >>sys.stderr, 'Generating Release Site.'
    print >>sys.stderr, 'Package Dir: %s' % (options.package_db,)
    print >>sys.stderr, 'Out file: %s' % (options.out_file,)
    db = PackageDatabase(options.package_db)
    db.load()
    # Load and render overview template.
    tpl = pyratemp.Template(filename=TPL_PATH)
    with open(options.out_file, 'wb') as f:
        f.write(tpl(FORMATS=FORMATS,
                    seqan_apps=db.seqan_apps,
                    seqan_library=db.seqan_library,
                    softwares=db.softwares,
                    sorted=sorted))
    # Load and render package template.
    tpl = pyratemp.Template(filename=PACKAGE_TPL_PATH)
    for sname, software in db.softwares.items():
        out_path = os.path.join(options.package_db, sname, 'index.html')
        print >>sys.stderr, 'Writing %s.' % out_path
        with open(out_path, 'wb') as f:
            f.write(tpl(FORMATS=FORMATS,
                        software=software,
                        sorted=sorted))
    # Write out RSS feeds for the packages.
    rss_writer = RssWriter(options.package_db, db, options.base_url)
    rss_writer.generate()


def main():
    parser = optparse.OptionParser()
    parser.add_option('-d', '--package-db', dest='package_db',
                      help='Path to directory with package files.')
    parser.add_option('-o', '--out-file', dest='out_file',
                      help='Path to the HTML file to generate.')
    parser.add_option('-b', '--base-url', dest='base_url',
                      help='Base URL.', default=BASE_URL)

    options, args = parser.parse_args()
    if args:
        parser.error('No arguments expected!')
        return 1
    if not options.package_db:
        parser.error('Option --package-db/-d is required!')
        return 1
    if not options.out_file:
        parser.error('Option --out-file/-o is required!')
        return 1

    return work(options)
