#!/usr/bin/python
"""Implementation of the SeqAn Doxygen dialect.
"""

import argparse
import logging
import os
import re
import sys

import file_mgr
import lexer
import dox_tokens
import dox_parser
import proc_doc
import raw_doc
import write_html
import migration

class FileNameSource(object):
    def __init__(self, paths):
        self.paths = paths
        self.extensions = ['.h', '.cpp']
        self.ignore = ['.svn']

    def generate(self):
        for p in self.paths:
            if os.path.isfile(p):
                yield os.path.abspath(p)
            for root, subFolders, files in os.walk(p):
                for f in files:
                    if not any([f.endswith(s) for s in self.extensions]):
                        continue
                    if any([f.startswith(s) for s in self.ignore]):
                        continue
                    yield os.path.join(root, f)


def doMain(args):
    # Parse all legacy files.
    import seqan.dddoc.core as core
    app = core.App()
    for path in args.legacy_doc_dirs:
        print 'Scanning %s...' % path
        app.loadFiles(path)
    migrated_doc = raw_doc.RawDoc()
    if args.legacy_doc_dirs:
        app.loadingComplete()
        migrated_doc.entries = migration.migrate(app.dddoc_tree)
        print 'migrated_doc.entries', [e.name.text for e in migrated_doc.entries]
    # Parse all normal input files.
    fmgr = file_mgr.FileManager()
    master_doc = raw_doc.RawDoc()
    master_doc.merge(migrated_doc)
    fns = FileNameSource(args.inputs)
    for filename in fns.generate():
        if args.debug:
            print 'Processing %s' % filename
        the_file = fmgr.loadFile(filename)
        lex = lexer.Lexer(dox_tokens.LEXER_TOKENS, skip_whitespace=False)
        for comment in the_file.comments:
            # TODO(holtgrew): Also give offset.
            lex.input(comment.text, filename, comment.line, comment.col)
            parser = dox_parser.Parser()
            parser.parse(lex)
            master_doc.merge(parser.documentation)
    # Generate documentation.
    logging.basicConfig(format='%(message)s', level=logging.DEBUG)
    logger = logging.getLogger()
    processor = proc_doc.DocProcessor(logger=logger, include_dir=args.base_dir)
    doc_proc = processor.run(master_doc)
    html_writer = write_html.HtmlWriter(doc_proc)
    html_writer.generateFor()
    return 0


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--lex-only', dest='lex_only', help='Lex only.',
                        default=False, action='store_true')
    parser.add_argument('--debug', dest='debug', help='Debug.',
                        default=False, action='store_true')
    parser.add_argument('-i', dest='inputs', help='Path to input.',
                        action='append', default=[])
    parser.add_argument('-lx', dest='legacy_demo_dirs', help='Path to legacy demos.',
                        action='append', default=[])
    parser.add_argument('-ldd', dest='legacy_doc_dirs', help='Path to legacy doc dirs.',
                        action='append', default=[])
    parser.add_argument('-b', '--base-dir', help='Base directory for @include.',
                        default='.', dest='base_dir')
    args = parser.parse_args()
    #if not args.inputs:
    #    parser.error('Missing input.')
    return doMain(args)


if __name__ == '__main__':
    sys.exit(main())

