#!/usr/bin/env python2

import datetime
import optparse
import os
import sys

import core
import html


HEADER = """
----------------------------------------------------------------------
  Dot.Dot.Doc documentation generator.
----------------------------------------------------------------------
""".strip()

CMD_HELP = """
Usage: %s <base path> [ <module> ]*
""".strip()


# Names of dirs with .dddoc files.
# TODO(holtgrew): Actually, this should be given on the cmd line.
DOC_DIRS=['pages', 'concepts']

class DDDocRunner(object):
    """Runner object for dddoc.

    Attrs:
      index_only  Boolean, true iff only index pages are built.
      doc_dirs    List of strings.  Names of directories with dddoc files.
    """
    
    def __init__(self, index_only=False, doc_dirs=[], out_dir='html',
                 demos_dir='.', cache_only=False, include_dirs=[]):
        """Initialize, arguments correspond to attributes."""
        self.index_only = index_only
        self.doc_dirs = doc_dirs
        self.out_dir = out_dir
        self.demos_dir = demos_dir
        self.cache_only = cache_only
        self.include_dirs = include_dirs

    def run(self, base_paths):
        """Run dddoc on the modules below the given path.

        Args:
          base_paths Paths to build the documentation for.

        Returns:
          Return code of the application.  Is 0 for no problem, and 1 on
          errors and warnings.
        """
        print 'Scanning modules...'
        app = core.App()
        if self.cache_only:
            for fn in app.cache.content.iterkeys():
                core.parseFile(fn, app.cache)
        else:
            # Scan some/all modules.
            for path in base_paths:
                os.path.normpath(path)
                app.loadFiles(path)
                
            # Scan doc directories.
            for doc_dir in self.doc_dirs:
                print 'Scanning %s...' % doc_dir
                app.loadFiles(doc_dir)

        app.loadingComplete()

        # Actually build the HTML files.
        print 'Creating HTML Documentation...'
        tpl_path = os.path.abspath(os.path.join(os.path.dirname(__file__), 'tpl'))
        res = html.createDocs(app.error_logger, app.dddoc_tree, tpl_path, self.out_dir,
                              self.include_dirs)

        # Done, print end message.
        print 'Documentation created/updated.'
        return res
    
    
def main(argv):
    """Program entry point."""
    print '%s\n' % HEADER

    start_time = datetime.datetime.now()
    
    # Parse arguments.
    parser = optparse.OptionParser()
    parser.add_option('-d', '--doc-dir', dest='doc_dirs', action='append',
                      default=[],
                      help=('Read .dddoc files from this directory.  '
                            'Can be given multiple times.'))
    parser.add_option('-o', '--out-dir', dest='out_dir', default='html',
                      help='Name of output directory.  Default: "html".')
    parser.add_option('-e', '--demos-dir', dest='demos_dir',
                      default='../projects/library/demos',
                      help=('Directory to demos. Default: '
                            '"../projects/library/demos".'))
    parser.add_option('-I', '--include-dir', dest='include_dirs',
                      action='append', default=[],
                      help='Paths to the directories for files and snippets.')
    parser.add_option('-c', '--cache-only', dest='cache_only', default=False,
                      action='store_true',
                      help='Ignore files if cache file exists.')
    options, args = parser.parse_args(argv)
    print 'doc dirs: %s' % ', '.join(options.doc_dirs)
    print
    
    # Show help if no arguments are given.
    if len(args) < 2:
        print CMD_HELP % args[0]
        return 1
    # Create application object and run documentation generation.
    app = DDDocRunner(index_only=False, doc_dirs=options.doc_dirs,
                      out_dir=options.out_dir,
                      include_dirs=options.include_dirs,
                      demos_dir=options.demos_dir,
                      cache_only=options.cache_only)
    res = app.run(args)

    elapsed = datetime.datetime.now() - start_time
    print >>sys.stderr, 'Took %d s' % elapsed.seconds

    return res
    

if __name__ == '__main__':
    sys.exit(main(sys.argv))
