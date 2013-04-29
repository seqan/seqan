#!/usr/bin/env python
"""Management for including file and snippets.
"""

__author__ = 'Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>'

import os
import os.path
import sys


class IncludeException(Exception):
    """Raised when a file or snippet could not be loaded."""


class IncludeManager(object):
    """Manages inclusion of files.

    @ivar base_path: String, path to the include path.
    @ivar snippet_cache: Cache from (path, snippet name) to snippet text.
    @ivar file_cache: Cache form path to file content.
    """

    def __init__(self, base_path):
        self.base_path = base_path
        self.file_cache = {}
        self.snippet_cache = {}

    def loadSnippet(self, path, name):
        """
        @raises IncludeException When there is an error loading the file or snippet.
        """
        if not self.snippet_cache.get((path, name)):
            self._loadSnippets(path)
        return self.snippet_cache[(path, name)]

    def loadFile(self, path):
        """
        @raises IncludeException When there is an error loading the file.
        """
        self._loadFile(path)
        return self.file_cache[path]

    def _loadFile(self, path):
        """
        @raises IncludeException When there is an error loading the file.
        """
        if self.file_cache.get(path):
            return
        try:
            with open(os.path.join(self.base_path, path)) as f:
                self.file_cache[path] = f.read()
        except IOError, e:
            raise IncludeException('Could not load file %s: %s' % (path, e))
        
    def _loadSnippets(self, path):
        """Load snippet

        @raises IncludeException When there is an error loading the file or snippet.
        """
        if not self.file_cache.get(path):
            self._loadFile(path)
        current_key = None
        current_lines = []
        fcontents = self.file_cache[path]
        for line in fcontents.splitlines():
            line = line.rstrip()  # Strip line ending and trailing whitespace.
            if line.strip().startswith('//![') and line.strip().endswith(']'):
                key = line.strip()[4:-1].strip()
                if key == current_key:
                    self.snippet_cache[(path, key)] = '\n'.join(current_lines)
                    current_lines = []
                    current_key = None
                else:
                    current_key = key
            elif current_key:
                current_lines = [line]
        if current_lines and current_key:
            self.snippet_cache[(path, current_key)] = '\n'.join(current_lines)
        