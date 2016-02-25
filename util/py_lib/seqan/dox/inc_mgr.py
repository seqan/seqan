#!/usr/bin/env python2
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

    @ivar base_paths: List of strings, paths to the include paths.
    @ivar snippet_cache: Cache from (path, snippet name) to snippet text.
    @ivar file_cache: Cache form path to file content.
    """

    def __init__(self, base_paths):
        self.base_paths = base_paths
        self.file_cache = {}
        self.snippet_cache = {}

    def loadSnippet(self, path, name):
        """
        @raises IncludeException When there is an error loading the file or snippet.
        """
        if not self.snippet_cache.get((path, name)):
            self._loadSnippets(path)
        if not self.snippet_cache.get((path, name)):
            raise IncludeException("Could not include snippet %s from file %s." % (path, name))
        return self.snippet_cache[(path, name)]

    def loadFile(self, path):
        """
        @raises IncludeException When there is an error loading the file.
        """
        self._loadFile(path)
        return self.file_cache[path]

    def resolvePath(self, path):
        """Translate a path passed to @include to an absolute path.

        The path is determined by looking at self.base_paths.
        """
        for base_path in self.base_paths:
            if os.path.exists(os.path.join(base_path, path)):
                return os.path.join(base_path, path)
        raise IncludeException('Could not find file %s for inclusion.' % path)

    def _loadFile(self, path):
        """
        @raises IncludeException When there is an error loading the file.
        """
        if self.file_cache.get(path):
            return
        full_path = self.resolvePath(path)
        try:
            with open(full_path) as f:
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
                current_lines.append(line)
        if current_lines and current_key:
            self.snippet_cache[(path, current_key)] = '\n'.join(current_lines)
        
