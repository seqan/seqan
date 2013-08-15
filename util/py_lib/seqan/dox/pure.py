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


# The expected HTML tags, useful for differentiating between F<T>::Type and real tags.
EXPECTED_TAGS = ['a', 'ul', 'ol', 'li', 'dl', 'dt', 'dd', 'em', 'i', 'b',
                 'strong', 'h1', 'h2', 'h3', 'h4', 'h5', 'h6', 'h7', 'tt',
                 'table', 'tbody', 'tr', 'th', 'td', 'caption', 'sup', 'img']


class FileNameSource(object):
    def __init__(self, paths):
        self.paths = paths
        self.extensions = ['.h', '.cpp', '.dox']
        self.ignore = ['.svn']

    def generate(self):
        for p in self.paths:
            if os.path.isfile(p):
                yield os.path.abspath(p)
            for root, subFolders, files in os.walk(p):
                for f in files:
                    if f.startswith('.'):
                        continue
                    if not any([f.endswith(s) for s in self.extensions]):
                        continue
                    if any([f.startswith(s) for s in self.ignore]):
                        continue
                    yield os.path.join(root, f)


def doMain(args):
    msg_printer = dox_parser.MessagePrinter(args.ignore_warnings_dirs)

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
            lex.input(comment.text, filename, comment.line + 1, comment.col, comment.offset_col)
            parser = dox_parser.Parser()
            try:
                parser.parse(lex)
            except dox_parser.ParserError, e:
                msg_printer.printParserError(e)
                return 1
            master_doc.merge(parser.documentation)
    # Generate documentation.
    logging.basicConfig(format='%(message)s', level=logging.DEBUG)
    logger = logging.getLogger()
    processor = proc_doc.DocProcessor(logger=logger, include_dirs=args.base_dirs,
                                      expected_tags=args.expected_tags,
                                      msg_printer=msg_printer)
    try:
        doc_proc = processor.run(master_doc)
    except dox_parser.ParserError, e:
        msg_printer.printParserError(e)
        return 1
    html_writer = write_html.HtmlWriter(doc_proc, args)
    html_writer.generateFor()


    msg_printer.printStats()
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
    parser.add_argument('--image-dir', dest='image_dirs', default=[],
                        action='append', help='Path to image directory.')
    parser.add_argument('-b', '--base-dir', help='Base directory for @include.',
                        default=['.'], dest='base_dirs', action='append')
    parser.add_argument('--expected-tags', help='Expected tags, warn about other tags.',
                        action='append', default=EXPECTED_TAGS)
    parser.add_argument('--ignore-warnings', help='Ignore warnings from directory.',
                        default=[], dest='ignore_warnings_dirs', action='append')
    args = parser.parse_args()
    #if not args.inputs:
    #    parser.error('Missing input.')
    return doMain(args)


if __name__ == '__main__':
    sys.exit(main())

