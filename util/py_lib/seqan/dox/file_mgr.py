#!/usr/bin/env python2
"""
Code for parsing comments from C++ files.
"""

import string
import re


class Comment(object):
    """Represents a comment in a file.

    @ivar file Back link to the File object.
    @ivar col, line, pos 0-based column, line, position.
    @ivar end_col, end_line, end_pos
          0-based end column, line, position.  Note that the end_line is the
          last line not the one after.
    @ivar offset_col 0-based column offset.
    @ivar text     str with the comment's content.
    @ivar raw_text star with comment's raw content.
    """

    def __init__(self, line, col, pos, end_line, end_col, end_pos, offset_line,
                 offset_col, text, raw_text):
        self.col = col
        self.line = line
        self.pos = pos
        self.end_col = end_col
        self.end_line = end_line
        self.end_pos = end_pos
        self.offset_line = offset_line
        self.offset_col = offset_col
        self.text = text
        self.raw_text = raw_text

    def __eq__(self, other):
        return (self.col == other.col and self.line == other.line and
                self.pos == other.pos and self.end_col == other.end_col and
                self.end_line == other.end_line and
                self.offset_col == other.offset_col and
                self.offset_line == other.offset_line and
                self.text == other.text and
                self.raw_text == other.raw_text)
        
    def __str__(self):
        return 'Comment(%d, %d, %d, %d, %d, %d, %d, %d, %s, %s)' % (
            self.line, self.col, self.pos, self.end_line, self.end_col,
            self.end_pos, self.offset_line, self.offset_col, repr(self.text),
            repr(self.raw_text))

    def __repr__(self):
        return str(self)
        

class File(object):
    """Represents one C++ files with simplified access to comments.
    """

    def __init__(self, path):
        """
        @ivar start_markers List of str with the start marker, e.g.
                            '**' for /** comments or '*!' for '/*$'.
        @ivar path          str with the path to the file.
        @ivar comments      List of Comment objects.
        """
        self.start_markers = ['*!']
        self.path = path
        self.comments = []

    def parse(self):
        with open(self.path, 'rb') as f:
            fcontents = f.read()
        self.comments = self._findComments(fcontents)

    def _findComments(self, text):
        """Parse C++ header as comments.

        @param text The text to parse.
        @returns List of Comment objects.
        """
        comments = []
        # Search for comments using regular expressions.
        def callback(match):
            if not any([match.group(0).startswith('/%s' % s) for s in self.start_markers]):
                return match.group(0)
            start, end = match.start(), match.end()
            start_line = string.count(text, '\n', 0, start)
            start_col = start - string.rfind(text, '\n', 0, start) - 1
            if start_col < 0:
                start_col = 0
            end_line = string.count(text, '\n', 0, end)
            end_col = end - string.rfind(text, '\n', 0, end) - 1
            offset_line = start_line + 1
            raw_text = match.group(0)
            # Compute offset column.
            idx = raw_text.find('\n') + 1
            for offset_col in xrange(idx, len(raw_text)):
                if raw_text[offset_col] not in [' ', '\t']:
                    break
            if raw_text[offset_col] == '*':
                offset_col += 1
                if raw_text[offset_col] in [' ', '\t']:
                    offset_col += 1
            offset_col -= idx

            proc_text = self._stripText(raw_text, offset_col)
            comments.append(Comment(start_line, start_col, start, end_line, end_col, end,
                                    offset_line, offset_col, proc_text, raw_text))
            return match.group(0)
        pattern = re.compile(
            r'//.*?$|/\*.*?\*/|\'(?:\\.|[^\\\'])*\'|"(?:\\.|[^\\"])*"',
            re.DOTALL | re.MULTILINE)
        re.sub(pattern, callback, text)
        return comments

    def _stripText(self, text, offset_col):
        """Strip leading stars etc."""
        res = []
        for line in text.splitlines()[1:-1]:
            res.append(line[offset_col:])
        return '\n'.join(res) + '\n'


class FileManager(object):
    """Handles loading of comments from files
    """
    
    def __init__(self, enable_cache=False, start_markers=['*!']):
        self.enable_cache = enable_cache
        self.start_markers = start_markers

    def loadFile(self, path):
        """Load a file from a path

        @returns File object with the parsed file.

        @raises IOError in case the file could not be opened for reading.
        """
        res = File(path)
        res.start_markers = self.start_markers
        res.parse()
        return res
