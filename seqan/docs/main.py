#!/usr/bin/env python
# import dddoc
# import dddoc_html
# import sys
# import os
# 
# # Command line: main.py <inpath> [ <module> ]*
# 
# if len(sys.argv) < 2: inpath = "../../projects/library"
# else: inpath =  sys.argv[1]
# 
# print("Welcome to Dot.Dot.Doc")
# 
# buildfull = (len(sys.argv) < 3)
# indexonly = not buildfull and (sys.argv[2] == 'indexonly')
# 
# if buildfull or indexonly:
#     print("Scanning " + inpath + "...")
#     os.path.normpath(inpath)
#     dddoc.loadFiles(inpath)
# else:
#     i = 2;
#     while (i < len(sys.argv)):
#         modulepath = inpath + "/seqan/" + sys.argv[i]
#         os.path.normpath(modulepath)
#         print("Scanning " + modulepath + "...")
#         dddoc.loadFiles(modulepath)
#         i += 1
# 
# 
# print("Scanning pages...")
# dddoc.loadFiles("pages")
# 
# print("Scanning concepts...")
# dddoc.loadFiles("concepts")
# 
# 
# dddoc.DATA.init()
# 
# print("Create HTML Documentation...")
# dddoc_html.createDocs("html", buildfull, indexonly)
# 
# if buildfull:
#     print("Documentation created.")
# else:
#     print("Documentation updated.")
# 
# #raw_input("press return")

#!/usr/bin/env python2.5

import optparse
import os
import sys

import dddoc
import dddoc_html


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
    
    def __init__(self, index_only=False, doc_dirs=[], out_dir='html', demos_dir='.'):
        """Initialize, arguments correspond to attributes."""
        self.index_only = index_only
        self.doc_dirs = doc_dirs
        self.out_dir = out_dir
        self.demos_dir = demos_dir

    def run(self, base_paths):
        """Run dddoc on the modules below the given path.

        Args:
          base_paths Paths to build the documentation for.

        Returns:
          Return code of the application.  Is 0 for no problem, and 1 on
          errors and warnings.
        """
        print 'Scanning modules...'
        dddoc_html.OUT_PATH = self.out_dir
        app = dddoc.App()
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
        # html_creator = dddoc_html.HtmlDocCreator(app, self.out_dir, not modules, self.index_only)
        # html_creator.run()
        dddoc_html.createDocs(self.out_dir, True, self.index_only)

        # Done, print end message.
        print 'Documentation created/updated.'
        return dddoc_html.WARNING_COUNT > 0
    
    
def main(argv):
    """Program entry point."""
    print '%s\n' % HEADER
    
    # Parse arguments.
    parser = optparse.OptionParser()
    parser.add_option('-d', '--doc-dir', dest='doc_dirs', action='append',
                      default=[],
                      help=('Read .dddoc files from this directory.  '
                            'Can be given multiple times.'))
    parser.add_option('-o', '--out-dir', dest='out_dir', default='html',
                      help='Name of output dirctory.  Default: "html".')
    parser.add_option('-e', '--demos-dir', dest='demos_dir',
                      default='../projects/library/demos',
                      help=('Directory to demos. Default: '
                            '"../projects/library/demos".'))
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
                      demos_dir=options.demos_dir)
    return app.run(args)
    

if __name__ == '__main__':
    sys.exit(main(sys.argv))
