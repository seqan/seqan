#!/usr/bin/env python2

import copy
import operator
import os
import os.path
import pickle
import string
import sys

# Constant for C++ files.
FILETYPE_CPP = 2
# Constant for DDDOC files.
FILETYPE_DDDOC = 1
# Constant for none of the above.
FILETYPE_OTHER = 0

SOURCE_ENCODING = 'iso8859-1'

# Extension of C++ files.
CPP_EXTS = ['c', 'C', 'cpp', 'CPP', 'c++', 'C++', 'h', 'H', 'hpp', 'HPP',
            'h++', 'H++']
# Extensions of DDDOC files.
DDDOC_EXTS = ['dddoc', 'DDDOC']

# List of ignored directory names.
IGNORED_DIRS = ['CSV', '.svn', 'seeds2', 'find2', 'cmake']

DATA = None
ID = 0

# Text attribute node keys.
TEXT_ATTRIBUTE_KEYS = set(['text', 'table', 'tableheader', 'code', 'console', 'section',
                           'subsection', 'image', 'contents', 'note', 'file', 'snippet',
                           'output'])

# Nodes having paths matching the following patterns are considered text
# container nodes.  Their children having only one more component which is in
# TEXT_ATTRIBUTE_KEYS are processed in a special way.  The last component is
# replaced with 'text' and their content is prefixed by "type=$key:" where $key
# is the original key.  The content of the text container nodes is prefixed with
# "type=$text:" and moved to a child with key 'text'.
TEXT_CONTAINER_PATHS = [
    'Indexpage.*.description',
    'Page.*.description',
    'Page.*.summary',
    'Page.*.glossary.*',
    'Function.*.example',
    'Function.*.summary',
    'Function.*.description',
    'Function.*.remarks',
    'Function.*.status',
    'Class.*.example',
    'Class.*.summary',
    'Class.*.description',
    'Class.*.remarks',
    'Class.*.status',
    'Metafunction.*.example',
    'Metafunction.*.summary',
    'Metafunction.*.description',
    'Metafunction.*.remarks',
    'Metafunction.*.status',
    'Memfunc.*.example',
    'Memfunc.*.summary',
    'Memfunc.*.description',
    'Memfunc.*.remarks',
    'Memfunc.*.status',
    'Memvar.*.example',
    'Memvar.*.summary',
    'Memvar.*.description',
    'Memvar.*.remarks',
    'Memvar.*.status',
    'Macro.*.example',
    'Macro.*.summary',
    'Macro.*.description',
    'Macro.*.remarks',
    'Macro.*.status',
    'Enum.*.example',
    'Enum.*.summary',
    'Enum.*.description',
    'Enum.*.remarks',
    'Enum.*.status',
    'Spec.*.example',
    'Spec.*.summary',
    'Spec.*.description',
    'Spec.*.remarks',
    'Spec.*.status',
    'Shortcut.*.example',
    'Shortcut.*.summary',
    'Shortcut.*.description',
    'Shortcut.*.remarks',
    'Shortcut.*.status',
    'Tag.*.example',
    'Tag.*.summary',
    'Tag.*.description',
    'Tag.*.remarks',
    'Tag.*.status',
    'Typedef.*.example',
    'Typedef.*.summary',
    'Typedef.*.description',
    'Typedef.*.remarks',
    'Typedef.*.status',
    'Demo.*.summary',
    'Demo.*.description',
    'Demo.*.remarks',
    'Demo.*.output',
    'Adaption.*.example',
    'Adaption.*.summary',
    'Adaption.*.description',
    'Adaption.*.remarks',
    'Adaption.*.status',
    'Concept.*.example',
    'Concept.*.summary',
    'Concept.*.description',
    'Concept.*.remarks',
    'Concept.*.status',
    ]

def _pathsMatch(path1, path2):
    """Compare two paths with wildcards."""
    if not type(path1) is list:
        path1 = splitKeys(path1[int(path1[0] == '.'):], '.')  # Strip leading '.', if any.
    if not type(path2) is list:
        path2 = splitKeys(path2[int(path2[0] == '.'):], '.')
    if len(path1) != len(path2):
        return False
    for i, p1 in enumerate(path1):
        p2 = path2[i]
        if not (p1 == '*' or p2 == '*' or p1 == p2):
            return False
    return True


def transformDddocEntry(entry):
    """Performs the text container node transformations.

    Returns list of entries to add if any.
    """
    for path in TEXT_CONTAINER_PATHS:
        if _pathsMatch(path, entry.path) and entry.content:  # Is text container.
            new_entry = copy.deepcopy(entry)
            new_entry.content = 'type=text:' + entry.content
            entry.content = ''
            return [new_entry]  # Done.
        if not entry.path[-1] in TEXT_ATTRIBUTE_KEYS:
            continue  # Skip if last component does not match.
        if not _pathsMatch(path, entry.path[:-1]):
            continue  # Skip non-matching path.
        # If we reach here, it is a text node.
        ## print 'TRANSFORMING ', entry
        last = entry.path[-1]
        entry.path = entry.path[:-1]
        entry.content = 'type=' + last + ':' + entry.content
        ## print '  to ', entry
        return []  # Done
    return []  # No updates.


class FileCache(object):
    """Simple file contents cache.

    Maps paths to (mtime, file contents) pairs.

    Attrs:

        path     Path to the cache file.
        content  Dict with cache content mapping file name to pair of mtime
                 and data associated with the cache.
    """
    
    def __init__(self, path):
        self.path = path
        self.content = {}
        self._tryLoad()

    def _tryLoad(self):
        try:
            with open(self.path, 'rb') as f:
                self.content = pickle.load(f)
        except:
            print >>sys.stderr, 'Could not load cache %s' % self.path
            return False
        print >>sys.stderr, 'Successfully loaded cache %s' % self.path
        return True

    def flush(self):
        """Store the cache to its file."""
        try:
            with open(self.path, 'wb') as f:
                pickle.dump(self.content, f)
        except:
            print >>sys.stderr, 'Could not store cache %s' % self.path
            return False
        print >>sys.stderr, 'Successfully stored cache %s' % self.path
        return True

    def has_key(self, key):
        """Returns True if the cache has data for this key."""
        return self.content.has_key(key)

    def isFresh(self, filename):
        """Returns True if the cache is fresh.

        The cache is fresh if the file at the given path is not newer than the
        data in the cache.
        """
        if not self.has_key(filename):
            return False
        mtime = os.stat(filename).st_mtime
        return mtime >= self.content[filename][0]

    def get(self, key, defaultValue=None):
        """Return content of the given entry."""
        return self.content.get(key, (None, defaultValue))[1]

    def set(self, filename, value):
        """Set cache content and mtime."""
        mtime = os.stat(filename).st_mtime
        self.content[filename] = (mtime, value)


class DddocEntry(object):
    def __init__(self, path, content, filename, line_no_begin, line_no_end):
        self.path = path
        self.content = content
        self.filename = filename
        self.line_no_begin = line_no_begin
        self.line_no_end = line_no_end

    def __str__(self):
        tpl = ('DddocEntry(path=%s, content=%s, filename=%s, line_no_begin=%s, '
               'line_no_end=%s)')
        values = (self.path, self.content, self.filename, self.line_no_begin,
                  self.line_no_end)
        return tpl % tuple(map(repr, values))

    def __repr__(self):
        return self.__str__()

    @classmethod
    def cmpPathLocation(klass, lhs, rhs):
        """Comparator, by entry path then filename and line number."""
        lhs_t = (lhs.path, lhs.filename, lhs.line_no_begin)
        rhs_t = (rhs.path, rhs.filename, rhs.line_no_begin)
        if lhs_t < rhs_t:
            return -1
        elif lhs_t > rhs_t:
            return 1
        else:
            return 0


def splitKeys(text, delimiters, limit=None, _cache={}):
    """Splitting that considers escaping of keys using quotes.

        >>> splitKeys('.Adaption.\'std::string\'.summary')
        ['', 'Adaption', '\'std::string\'', 'summary']
    """
    if '\u0001' in text:
        text = text.split('\u0001', 1)[0]  # Remove optional label, used in inheritance.
    if _cache.has_key((text, delimiters)):
        return _cache[(text, delimiters)]
    count = 0
    current = []
    result = []
    str_delimiter = None
    for i in range(0, len(text)):
        # Handle text in strings.
        if str_delimiter:
            if text[i] == str_delimiter:
                str_delimiter = None
            current.append(text[i])
            continue
        elif text[i] in '\'"':
            str_delimiter = text[i]
            current.append(text[i])
            continue
        # Handle non-in-string text.
        if text[i] in delimiters:
            result.append(''.join(current))
            current = []
            count += 1
            if limit and count >= limit:
                result.append(text[i+1:])
                _cache[(text, delimiters)] = result
                return result
        else:
            current.append(text[i])
    result.append(''.join(current))
    _cache[(text, delimiters)] = result
    return result


def cleanPath(path_arr):
    """Takes a list with a path and cleans its element.

    Cleaning its element currently only consists in removing singel and double
    quotes.
    """
    def _cleanPathElement(x):
        return x.strip().replace('\'', '').replace('"', '')
    return map(_cleanPathElement, path_arr)


class FileLoader(object):
    """File loader helper class.

    Attrs:

        cache    FileCache to use for caching.
        entries  List of DddocEntries objects.
    """
    
    def __init__(self, cache):
        self.cache = cache
        self.entries = []

    def _loadDDDOCFile(self, filename, cache):  # TODO(holtgrew): Make Top-Level Function?
        # Try to load from cache.
        if cache.isFresh(filename):
            return cache.get(filename)

        # Load file.
        with open(filename, 'rb') as f:
            text = [x.decode(SOURCE_ENCODING).encode("ascii", "xmlcharrefreplace") for x in f.readlines()]
        cache.set(filename, text)
        return text

    def _loadCPPFile(self, filename, cache):  # TODO(holtgrew): Make Top-Level Function?
        if cache.isFresh(filename):
            return cache.get(filename)

        # TODO(holtgrew): This looks overly complicated.
        f = open(filename)
        lines = [x.decode(SOURCE_ENCODING).encode("ascii", "xmlcharrefreplace") for x in f.readlines()]
        f.close()

        ret = []

        #test for SEQAN_NO_DDDOC
        for line in lines:
            if line.find("SEQAN_NO_DDDOC") >= 0:
                cache.set(filename, ret)
                return ret;


        incomment = False
        innextcomment = False
        inextract = False

        for line in lines:
            line = line.rstrip()
            str_line = ""
            if len(line) == 0:
                if not innextcomment and not incomment: 
                    str_line = "."
                else: 
                    str_line = " "

            while len(line) > 0 :
                if innextcomment:
                    if line[len(line)-1] == "\\" :
                        if inextract: str_line += line[: len(line)-1]
                    else:
                        if inextract: str_line += line
                        innextcomment = False
                    break

                elif incomment:
                    pos1 = line.find("*/")
                    if pos1 < 0:
                        if inextract: str_line += line;
                        break;
                    else:
                        if inextract: 
                            str_line += line[:pos1];
                            line = line[pos1 + 3:];
                        else:
                            line = line[pos1 + 2:];
                        incomment = False;

                else:
                    pos1 = line.find("/*")
                    pos2 = line.find("//")
                    pos3 = line.find('"')
                    if (pos1 >= 0) and ((pos2 < 0) or (pos1 < pos2)) and ((pos3 < 0) or (pos1 < pos3)):
                        pos9 = line.find("*/", pos1 + 2)
                        if (len(line) > pos1 + 2):
                            inextract = (line[pos1 + 2] == "/") or (line[pos1 + 2] == "*")
                        else:
                            inextract = False
                        if pos9 < 0 : 
                            if inextract: str_line += line[pos1 + 3:]
                            incomment = True
                            break
                        else: 
                            if inextract: 
                                str_line += line[pos1 + 3: pos3]
                                line = line[pos9 + 3:]
                            else:
                                line = line[pos9 + 2:]

                    elif (pos2 >= 0) and ((pos3 < 0) or (pos2 < pos3)):
                        pos2b = pos2 + 2;
                        while ((pos2b < len(line)) and ((line[pos2b] == "/") or (line[pos2b] == "*"))):
                            pos2b += 1
                        inextract = (pos2b > pos2 + 2)
                        if line[len(line)-1] == "\\" :
                            if inextract: str_line += line[pos2b: len(line)-1]
                            innextcomment = True
                        else:
                            if inextract: str_line += line[pos2b:]
                        break

                    elif pos3 >= 0:
                        pos9 = line.find('"', pos3 + 2)
                        if pos9 < 0:
                            line = line[pos9+1:]
                            break
                        else:
                            break

                    else:
                        break

            ret = ret + [str_line]

        cache.set(filename, ret)
        return ret

    def _getFileType(self, filename):  # TODO(holtgrew): Make Top-Level Function?
        """Determines file type from filename.

        Determines the file type from the extension of the given filename.

        >>> getFileType('test.cpp') == FILETYPE_CPP
        True
        >>> getFileType('path/file.h') == FILETYPE_CPP
        True
        >>> getFileType('test.dddoc') == FILETYPE_DDDOC
        True

        Args:

            filename  Filename to parse.

        Returns:

            One of {FILETYPE_CPP, FILETYPE_DDDOC, FILETYPE_OTHER}, depending
            on the extension of filename.
        """
        # Get file extension.
        base, ext = os.path.splitext(filename)
        if ext[1:] in CPP_EXTS:
            return FILETYPE_CPP
        elif ext[1:] in DDDOC_EXTS:
            return FILETYPE_DDDOC
        else:
            return FILETYPE_OTHER

    def _loadFile(self, filename):
        """Load the file with the given filename.

        The line is then split into DDDoc entries, unwrapping entries that span
        more than one line.  Finally, the keys are expanded, and surrounding
        whitespace is stripped.
        """
        ## print filename
        # Load file contents, through a cache.
        file_type = self._getFileType(filename)
        if file_type == FILETYPE_CPP:
            text = self._loadCPPFile(filename, self.cache)
        elif file_type == FILETYPE_DDDOC:
            text = self._loadDDDOCFile(filename, self.cache)
        else:
            raise Error("Unknown file type of file %s." % path)
        text.append('.')
        ## print 'LOADING', filename
        ## print '\n'.join(text)

        # Process all lines in the input, join lines that do not begin with a
        # dot with the previous ones.  This allows the wrapping of lines.
        str = False
        dddoc_entries = []  # [(path, filename, begin line no, end line no)]
        line_no_begin, line_no_end = 1, 1
        for line in text:
            ## if line and line != '.':
            ##     print 'LINE', line
            line_no_end += 1
            if not line:
                continue
            if line[0] == '.':
                if str is not False and str[0] == '.' and str != '.' and str.strip():  # Skip empty dummy lines.
                    dddoc_entries.append([str, filename, line_no_begin, line_no_end])
                    ## print dddoc_entries[-1]
                line_no_begin = line_no_end
                str = line
                if str == '.':
                    str = False
            elif str:
                if str[-1] != '\n':
                    str += '\n'
                str += line

        # Now, expand the keys of dddoc_entries, e.g. dddoc_entries[i][0].
        # TODO(holtgrew): Consider escaping of keys here.
        stack = []
        stack_len_sum = 0
        for entry in dddoc_entries:
            ## print 'ENTRY', entry
            ## print 'stack=%s' % (stack)
            # Split out $key:$value of the entry and $the.$path.$elements from $key.
            maybe_pair = splitKeys(entry[0].strip(), ':', 1)
            if len(maybe_pair) == 2:
                key, value = splitKeys(entry[0].strip(), ':', 1)
            else:
                key, value = entry[0].strip(), ''
            path = splitKeys(key, '.')[1:]
            # Count empty entries in the path.
            ## print ' ', path
            empty_count = reduce(operator.add, [1 for x in path if not x], 0)
            ## print '  empty_count', empty_count
            if empty_count <= len(stack):
                stack = stack[:empty_count]
                stack_len_sum = reduce(operator.add, map(len, stack), 0)
            stack.append(path[empty_count:])
            stack_len_sum += len(stack[-1])
            path = reduce(operator.add, stack, [])
            # Remove any leading and trailing whitespace from value and compute
            # updated begin and end line no.
            line_count = len(value.splitlines())
            value_no_leading = value.lstrip()
            line_count2 = len(value_no_leading.splitlines())
            line_no_begin = entry[2] + line_count - line_count2
            value_no_trailing = value_no_leading.rstrip()
            line_count3 = len(value_no_trailing.splitlines())
            line_no_end = entry[3] - line_count2 + line_count3
             
            # Store the DDDoc entry.
            if path:
                self.entries.append(DddocEntry(cleanPath(path), value_no_trailing, filename, line_no_begin, line_no_end))
                new_entries = transformDddocEntry(self.entries[-1])
                ## if new_entries:
                ##     print 'NEW ENTRIES', new_entries
                self.entries += new_entries
                ## print self.entries[-1]

    def run(self, search_path):
        """Call parseFile() on files.

        All files below search_path will be searched that have file type
        FILETYPE_CPP or FILETYPE_DOC as determined by getFileType().
        Directories with names of IGNORED_DIRS are skipped.

        Args:
            search_path  String, path to search files under.
        """
        for root, dirs, files in os.walk(search_path):
            # Parse all files.
            for file in files:
                if os.path.basename(file).startswith('.'):
                    continue  # Skipp hidden files.
                path = os.path.join(root, file)
                if self._getFileType(path) in [FILETYPE_CPP, FILETYPE_DDDOC]:
                    self._loadFile(path)
            # Exclude ignored diretories.
            for ignored in IGNORED_DIRS:
                if ignored in dirs:
                    dirs.remove(ignored)


class DddocTreeNode(object):
    """Represents one entry in the DddocTree.

    Attrs:

        tree      The DddocTree that the node belongs to.
        key       The key of this child, last element of path.
        path      The full path to the child.
        entry     Range [beg, end) of DddocEntry that this node represents.
        children  dict with the children as key/value pairs.
        texts     Array of strings with the texts.
    """
    
    def __init__(self, tree, key, path, entry, children={}):
        self.tree = tree
        self.key = key
        self.path = path
        self.entry = entry
        self.children = children
        self.texts = []

    def text(self, spacer=' '):
        return spacer.join(self.texts)

    def __str__(self):
        """Returns dump for the whole tree in a user-readable manner."""
        def _str(node, level=0, prefix=''):
            space = '  ' * level
            if prefix:
                prefix = prefix + ' --> '
            res = '%s %sDddocTreeNode(key=%s, texts=%s)' % (space, prefix, repr(node.key), repr(node.texts))
            for k, child in node.children.iteritems():
                res += '\n' + _str(child, level + 1, k)
            return res
        return _str(self)

    def dump(self, stream=sys.stdout):
        """Debug recursive dumping of a tree node."""
        print >>stream, self


class DddocTree(object):
    """Tree with the information from the DDDoc contents.

    Attrs:

        entries  The raw entries.
        root     The root DddocTreeNode.
        glossary_nodes  List of nodes that contain glossary entries.  Built
                        in finalize().
    """
    
    def __init__(self, entries):
        self.entries = entries
        #for e in self.entries:
        #    print e
        self.root = DddocTreeNode(self, 'ROOT', [], (0, 0), self._buildSubtree([], 0, len(entries), 0))
        self.cache = None
        self.glossary_nodes = []
        ## self.root.dump()
        ## for entry in self.entries:
        ##     print entry.path, entry.content

    def _enableFindCache(self):
        if self.cache is None:
            self.cache = {}

    def finalize(self):
        """Called after tree will not be modified any more.

        Enables caching and builds some indices.
        """
        self._enableFindCache()
        print >>sys.stderr, 'Indexing Glossary Pages'
        if 'Page' in self.root.children:
            for key, node in self.root.children['Page'].children.iteritems():
                if 'glossary' in node.children:
                    self.glossary_nodes.append(node.children['glossary'])
                    print >>sys.stderr, '  Found Page.%s' % node.key

    def _buildSubtree(self, path, begin_index, end_index, level):
        # First, identify the entries belonging to each node (entry.path[i] are
        # equal for i = level, inductively, also i <= level).
        prev_key = None
        prev_beg = None
        subseqs = []
        for i in range(begin_index, end_index):
            if prev_key != self.entries[i].path[level]:
                if prev_key != None:
                    subseqs.append((prev_beg, i))
                prev_key = self.entries[i].path[level]
                prev_beg = i
        if prev_key != None and prev_beg != end_index:  # Handle last.
            subseqs.append((prev_beg, end_index))
        # Now, subseqs contains a sequence of contiguous half-open intervals.
        # Each contains the data for one tree node.  There is a possibly empty
        # sequence of leading entries with paths of length level + 1 containing
        # the data for the current level node.  The rest is for the level below.
        result = {}
        for (b, c) in subseqs:
            assert b != c
            # Split into entries for this and for next level: [a, b); [b, c).
            a = b  # [a, b) will be for this vertex.
            while b < c and len(self.entries[b].path) == level + 1:
                b += 1
            # Compute node.
            path = self.entries[a].path[:(level + 1)]
            key = path[level]
            node = DddocTreeNode(self, key, path, (a, b))
            ## print 'new node', key
            for i in range(a, b):
                if self.entries[i].content:
                    node.texts.append(self.entries[i].content)
            # Compute subtree.
            node.children = self._buildSubtree(path, b, c, level + 1)
            result[key] = node
        return result

    def find(self, path):
        """Query tree for a DddocTreeNode.

        The argument path can either be a dot-separated string or a list with
        this information.  If path is a string then one optional leading dot is
        optional.  Returns None if nothing could be found.

            tree.find(['path', 'to', 'node'])
            tree.find('path.to.node')
            tree.find('.path.to.node')
        """
        ## print 'FIND(%s)' % repr(path)
        # Try to retrieve from cache if there is a cache.
        if not self.cache is None:
            if not type(path) is str:
                key = '.'.join(path)
            else:
                key = path
            if self.cache.has_key(key):
                return self.cache[key]
        # Split path if is string, ignore leading dot if any.
        if type(path) is str:
            path = splitKeys(path, '.')
            if path and path[0] == '':
                path = path[1:]
        # Now, query the tree.
        def findRecurse(node, path):
            """Helper function that searches for the node with given path."""
            if not path:
                return node
            if not node.children.has_key(path[0]):
                return None
            return findRecurse(node.children[path[0]], path[1:])
        res = findRecurse(self.root, path)
        if not self.cache is None:
            self.cache['.'.join(path)] = res
        return res

# Paths where the inline summary is moved into a .summary child.  See
# documentation of processInlineSummaries() for details.
SUMMARY_PATHS = [
    '*.*.param.*',
    '*.*.returns',
    '*.*.tag.*',
    '*.*.value.*',
    '*.*.returns.param.*',  # TODO(holtgrew): Used for metafunctions, could be improved.
    'Adaption.*',
    'Class.*',
    'Concept.*',
    'Demo.*',
    'Enum.*',
    'Function.*',
    'Macro.*',
    'Memfunc.*',
    'Metafunction.*',
    'Shortcut.*',
    'Spec.*',
    'Tag.*',
    ]

# TODO(holtgrew): Also use for generateAutomaticReferences()

def _matchTreesInNode(tree, node, path, func, block_paths=[['globals']], level=0):
    """Calls func on nodes matching path."""
    ## print '  ' * level, '_matchTreesInNode(tree', node.path, path, func, level, ')'
    if path:
        if path[0] == '*':
            for child in node.children.itervalues():
                _matchTreesInNode(tree, child, path[1:], func, block_paths, level+1)
        else:
            if node.children.has_key(path[0]):
                _matchTreesInNode(tree, node.children[path[0]], path[1:], func, block_paths, level+1)
    else:
        for block_path in block_paths:
            ## print node.path[:len(block_path)], '==', block_path
            if node.path[:len(block_path)] == block_path:
                return  # Path is blocked.
        func(node)


def processInlineSummaries(tree, paths):
    """Move inline documentation to .summary subpaths.

    The nodes matching the values in path are such that inline values are moved
    to .summary subnodes for greater consistency and lower variance.

    E.g. the following:

        .Function.f.param.x:This is param x.

    will be transformed into

        .Function.f.param.x
        ..summary:This is param x.
    """
    # First, collect nodes for the summary transfer.
    collected_nodes = set()
    def f(node):
        if node.texts:
            collected_nodes.add(node)
    for path in paths:
        _matchTreesInNode(tree, tree.root, splitKeys(path, '.'), f)
    # Then, move the inline summaries into a summary node.
    for node in collected_nodes:
        if not 'summary' in node.children:  # Create node if necessary.
            summaryNode = DddocTreeNode(tree, 'summary', node.path + ['summary'], (-2,-2))
            node.children['summary'] = summaryNode
        node.children['summary'].texts += node.texts
        node.texts = []


def generateAutomaticReferences(tree):
    """Interpret the globals.relations entries."""
    print >>sys.stderr, 'Generating Automatic References'
    relations_node = tree.find('globals.relations')
    if not relations_node:
        return  # Empty, do nothing.

    # We first collect all automatic links, scheduled to be added later.
    additions = []
    def appendToAdditions(node):
        for node_path in node.texts:
            node_path = splitKeys(node_path, '.')
            ## print '  ' * level, ' ', node_path
            res = tree.find(node_path)
            ## print '  ' * level, ' ', res is not None
            if not res:
                continue  # Not found, Skip  # TODO(holtgrew): Warning?
            additions.append((res.path + [key], '.'.join(node.path[:2])))
    for key, node in relations_node.children.iteritems():
        ## print 'GENERATE', key, node
        for txt in node.texts:
            path = splitKeys(txt, '.')
            _matchTreesInNode(tree, tree.root, splitKeys(txt, '.'), appendToAdditions)

    # Now, add these additions.  This circumvents problems leading to infinite
    # recursions.
    for path, text in additions:
        res = tree.find(path)
        if not res:
            parent = tree.find(path[:-1])
            assert parent
            res = DddocTreeNode(tree, path[-1], path, None)
            parent.children[path[-1]] = res
        if not text in res.texts:
            res.texts.append(text)


def generateInheritedElements(tree):
    """Push through inheritances."""
    print >>sys.stderr, 'Linking Inherited Entities'
    inherit_node = tree.find('globals.inherit')
    # Contains children: $TARGET_FIELD:$THROUGH_FIELD.$SOURCE_FIELD

    all_paths = set()
    depends_on = {}
    inheritance_rules = []

    # First build a dependency graph.
    for target_field, child in inherit_node.children.items():
        for txt in child.texts:
            arr = splitKeys(txt, '.')
            through_field = arr[0]
            if len(arr) > 1:
                source_field = arr[1]
            else:
                source_field = target_field
            inheritance_rules.append((target_field, through_field, source_field))
            def registerDependencies(node):
                all_paths.add('.'.join(node.path))
                if not through_field in node.children:
                    return
                for path in node.children[through_field].texts:
                    pth = '.'.join(node.path)
                    depends_on.setdefault(pth, set()).add(path)
            _matchTreesInNode(tree, tree.root, ['*', '*'], registerDependencies)
    ## print 'ALL PATHS', all_paths

    # Now, push through references by inheritance for all paths that are not
    # linked to and not completed yet.
    done = set()
    to_do = all_paths - done - set(depends_on.keys())
    while to_do:
        # Process all that are not completed and have no dependencies.
        if not to_do:
            raise Exception('Could not process all dependencies. Cyclic dependencies?')
        # Actually perform the preprocessing.
        for target_path in to_do:
            for target_field, through_field, source_field in inheritance_rules:
                target_node = tree.find(target_path)
                if not through_field in target_node.children:
                    continue  # Skip if no source children.
                ## print 'TRYING', target_path, through_field, source_field
                for source_path in target_node.children[through_field].texts:
                    source_node = tree.find(source_path)
                    if not source_field in source_node.children:
                        continue  # Skip if no source field.
                    for path in source_node.children[source_field].texts:
                        if not '\u0001' in path:  # We use this ugly hack to add the inheritance source here.
                            path = path + '\u0001' + '.'.join(source_node.path)
                        # If necessary then create child in target node.
                        if not target_field in target_node.children:
                            target_node.children[target_field] = DddocTreeNode(tree, target_field, target_node.path + [target_field], source_node.children[source_field].entry)
                        # Copy over path.
                        target_node.children[target_field].texts.append(path)
                        ## print '  appending', path
                
        # Clear out the stuff that we completed.
        to_delete = []
        for key in depends_on:  # Clear out all done.
            depends_on[key] -= to_do
            if not depends_on[key]:
                to_delete.append(key)
        for key in to_delete:
            del depends_on[key]
        done |= to_do  # Add done.
        to_do = all_paths - done - set(depends_on.keys())


def removeDuplicateTexts(tree):
    """Remove duplicates from texts members.

    Suffixes starting with '\u0001' are ignored for the comparisons
    and strings with these suffixes are preferred.
    """
    ##print 'remove duplicates'
    def recurse(node):
        in_cleaned = {}
        cleaned = []
        for txt in node.texts:
            clean = txt
            pos = txt.find('\u0001')
            if pos != -1:
                clean = txt[:pos]
            ##print cleaned, repr(clean)
            if clean in in_cleaned:
                if '\u0001' in clean and not '\u0001' in cleaned[in_cleaned[clean]]:
                    cleaned[in_cleaned[clean]] = txt
            else:
                in_cleaned[clean] = len(cleaned)
                cleaned.append(txt)
        node.texts = cleaned
        for child in node.children.itervalues():
            recurse(child)
    for child in tree.root.children.itervalues():
        recurse(child)


# TODO(holtgrew): If needed, this could easily be generalized.
def buildByTypeAndCatIndex(tree):
    """Build an index into the given DddocTree.

    The index will be a two-dimensional dict, mapping (first element of path,
    value of cat field) to a list of nodes in the DddocTree.
    """
    result = {}
    def recurse(result, path, node):
        ## print path, node.path
        if len(path) == 2:
            if node.children.has_key('cat'):
                for cat in node.children['cat'].texts:
                    result.setdefault(path[0], {}).setdefault(cat, []).append(node)
            else:
                result.setdefault(path[0], {})[path[1]] = node
        if len(path) < 2:
            for key, child in node.children.iteritems():
                recurse(result, path + [key], child)
    for key, child in tree.root.children.iteritems():
        recurse(result, [key], child)
    ## for k1, v1 in result.iteritems():
    ##     for k2, v2 in v1.iteritems():
    ##         print 'k1=%s\tk2=%s\tv=%s' % (k1, k2, [x.path for x in v2])
    return result


class ErrorLogger(object):
    def __init__(self):
        self.error_count = 0
    def invalidReference(self, txt, locations):
        self.error_count += 1
        if not locations:
            print >>sys.stderr, 'ERROR: Invalid Reference %s in unknown location (sorry).' % txt
        else:
            print >>sys.stderr, 'ERROR: Invalid Reference %s in one of the following locations:' % txt
            for filename, line in locations:
                print >>sys.stderr, '  %s:%s' % (filename, line)


class App(object):
    """Application object for DDDoc.

    Provides a facade to the functionality of the core module.

    Usage:
       app = App()
       app.loadFiles([<files>])
       app.loadFiles([<files>])
       app.loadingComplete()

    Attrs:
      data       The global state Data object.
    """

    def __init__(self):
        """Initialize object members."""
        self.cache = FileCache('dddoc_cache.bin')
        self.file_loader = FileLoader(self.cache)
        self.dddoc_tree = None
        self.error_logger = ErrorLogger()

    def loadFiles(self, path):
        """Load the files with the given file name."""
        self.file_loader.run(path)

    def loadingComplete(self):
        """Initialize data object.

        This method is called after all calls to loadFiles().
        """
        # Save the cache to disk again.
        self.cache.flush()
        # Sort Dddoc Entries and build tree.
        self.file_loader.entries.sort(cmp=DddocEntry.cmpPathLocation)
        self.dddoc_tree = DddocTree(self.file_loader.entries)
        # Generate automatic references.
        generateAutomaticReferences(self.dddoc_tree)
        # Perform inheritance as configured in global configuration.
        generateInheritedElements(self.dddoc_tree)
        # Clean duplicates from 'texts' members
        removeDuplicateTexts(self.dddoc_tree)
        # Move inline summaries into .summary children.
        processInlineSummaries(self.dddoc_tree, SUMMARY_PATHS)
        # Finally, after all modifications, enable caching and build indices in
        # tree.
        self.dddoc_tree.finalize()

    def getNextId(self):
        """Returns an identifier.

        Each id is only returned once.
        """
        assert False, "For future use."
        self.next_id += 1
        return self.next_id - 1
