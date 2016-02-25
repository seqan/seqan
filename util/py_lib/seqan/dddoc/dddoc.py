#!/usr/bin/env python2

import os
import os.path
import copy
import pickle
import string
import sys

# Constant for C++ files.
FILETYPE_CPP = 2
# Constant for DDDOC files.
FILETYPE_DDDOC = 1
# Constant for none of the above.
FILETYPE_OTHER = 0

# Extension of C++ files.
CPP_EXTS = ['c', 'C', 'cpp', 'CPP', 'c++', 'C++', 'h', 'H', 'hpp', 'HPP',
            'h++', 'H++']
# Extensions of DDDOC files.
DDDOC_EXTS = ['dddoc', 'DDDOC']

# List of ignored directory names.
IGNORED_DIRS = ['CSV', '.svn', 'seeds2', 'find2', 'cmake']

DATA = None
ID = 0

class DddocCache(object):
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
        try:
            with open(self.path, 'wb') as f:
                pickle.dump(self.content, f)
        except:
            print >>sys.stderr, 'Could not store cache %s' % self.path
            return False
        print >>sys.stderr, 'Successfully stored cache %s' % self.path
        return True

    def has_key(self, key):
        return self.content.has_key(key)

    def isFresh(self, filename):
        if not self.has_key(filename):
            return False
        mtime = os.stat(filename).st_mtime
        return mtime >= self.content[filename][0]

    def get(self, key, defaultValue=None):
        return self.content.get(key, (None, defaultValue))[1]

    def set(self, filename, value):
        mtime = os.stat(filename).st_mtime
        self.content[filename] = (mtime, value)
        

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
        global DATA
        DATA = Data([], 0)
        self.data = DATA
        self.next_id = ID
        self.cache = DddocCache('dddoc_cache.bin')

    def loadFiles(self, filenames):
        """Load the files with the given file name."""
        loadFiles(filenames, self.cache)

    def loadingComplete(self):
        """Initialize data object.

        This method is called after all calls to LoadFiles().
        """
        self.cache.flush()
        self.data.init()

    def getNextId(self):
        """Returns an identifier.

        Each id is only returned once.
        """
        assert False, "For future use."
        self.next_id += 1
        return self.next_id - 1
    


class Line:
    def __init__(self, _nodes, _text, file_name, line_no):
        global ID
        
        self.nodes = _nodes
        self._text = _text.rstrip()
        self.id = ID
        self.file_name = file_name
        self.line_no = line_no
        ID += 1
        
    def __repr__(self):
        return 'Line(%s, %s, %s, %s)' % (repr(self.nodes), repr(self._text), repr(self.file_name), repr(self.line_no))
        
    def name(self, index = 0):
        if len(self.nodes) > index:
            return self.nodes[index]
        else:
            return '(unknown)'
            
    def text(self):
        return self._text
        
            
class Data:
    def __init__ (self, _lines, _level):
        self.lines = _lines
        self.level = _level
        self.cache = {}
        
    def __repr__(self):
        return 'Data(\n  %s,\n%d)' % (',\n  '.join(map(repr, self.lines)), self.level)

    def init(self):
        self.lines.sort(sortLineCompare)
        
        relations = self["globals.relations"].at_level(1).lines
        for relation in relations:
            to = relation.name(2)
            arr = splitName(relation.text())
            if len(arr) > 0:
                findRelation(self, arr, to)
                
        self.lines.sort(sortLineCompare)
        
    def __getitem__(self, str):
        return self.find(str)
        
    def find(self, str):
        """Find all lines below the given path.

        Args:
          str  String with a dot separated path to the lines to find.

        Returns:
          Data object with the lines below the given path.
        """
        # If possible, return from cache.
        if self.cache.has_key(str):
            return self.cache[str]
        
        arr = splitName(str)
        lines = []

        # TODO(holtgrew): Keeping linear code for linear search to be able to fix things quickly.
        if True:
            # Use binary search for finding the first hit.
            def isHit(line, arr, level):
                """Return True iff arr matches the key of line, from the given level on."""
                i = 0
                while (i < len(arr)) and (i + level < len(line.nodes)) and (arr[i] == line.nodes[i + level]):
                    i += 1
                return i == len(arr)

            # Use binary search to find the first hit.
            query = arr
            lo = 0
            hi = len(self.lines)
            while lo < hi:
                mid = (lo + hi) // 2
                slice = self.lines[mid].nodes[self.level:self.level + len(query)]
                if slice < query:
                    lo = mid + 1
                else:
                    hi = mid
            result = lo

            # Output all consecutive hits, if any.
            if result < len(self.lines) and isHit(self.lines[result], arr, self.level):
                for line in self.lines[result:]:
                    if isHit(line, arr, self.level):
                        lines.append(line)
                    else:
                        break
        else:
            # Use linear search for finding the first hit.  Deactivated for now.
            maxi = 0
            for line in self.lines:
                i = 0
                while (i < len(arr)) and (i + self.level < len(line.nodes)) and (arr[i] == line.nodes[i + self.level]):
                    i += 1
                if i == len(arr):
                    lines.append(line)

                elif maxi > i:
                    break    
                maxi = i

            
        data = Data(lines, self.level + len(arr))
        # Cache result.
        self.cache[str] = data
        return data

    def at_level(self, level = 0):
        lines = []
        for line in self.lines:
            if len(line.nodes) == self.level + level:
                lines.append(line)
            
        data = Data(lines, self.level + level)
        return data
        
    def sub_level(self, level = 1):
        lines = []
        for line in self.lines:
            if len(line.nodes) >= self.level + level:
                lines.append(line)
            
        data = Data(lines, self.level + level)
        return data
        
    def by_occ(self):
        lines = copy.copy(self.lines)
        lines.sort(sortLinesByOcc)
        data = Data(lines, self.level)
        return data        

    def empty(self):
        return len(self.lines) == 0
        
    def name(self, index = 0):
        if len(self.lines) > 0:
            return self.lines[0].name(index)
        return '(empty)'
        
    def text(self):
        str = ''
        for line in self.at_level(0).lines:
            if (str != ''): str += '\n'
            str += line.text()
        return str
        
    def keys(self, level = 0):
        dict = {}
        for line in self.lines:
            if len(line.nodes) > self.level + level:
                dict[line.nodes[self.level + level]] = 1
                
        arr = dict.keys()
        arr.sort()
        return arr
            
    def keys_by_occ(self, level = 0):
        dict = {}
        for line in self.lines:
            if len(line.nodes) > self.level + level:
                key = line.nodes[self.level + level]
                if not dict.has_key(key) or (dict[key] > line.id):
                    dict[key] = line.id

        dict2 = {}
        for key in dict:
            dict2[dict[key]] = key
            
        arr2 = dict2.keys()
        arr2.sort()
        
        arr = []
        for i in arr2:
            arr.append(dict2[i])
        
        return arr
        
         
################################################################################
         
def sortLineCompare(left, right):
    l = left.nodes
    r = right.nodes
    
    i = 0
    while (i < len(l)) and (i < len(r)):
        ret = cmp(l[i], r[i])
        if ret != 0: 
            return ret
        i += 1
        
    if len(l) < len(r): return -1
    elif len(l) > len(r): return 1
    elif left.id < right.id: return -1
    else: return 1
       
################################################################################
         
def sortLinesByOcc(left, right):
    if left.id < right.id: return -1
    else: return 1
         
################################################################################

def findRelation(data, arr, to):
    global DATA
    
    if len(arr) > 0:
        if (arr[0] == '*'): 
            sub_data = data.sub_level(1)
        else: 
            lines = []
            for line in data.lines:
                if line.name(data.level) == arr[0]:
                    lines.append(line)
            sub_data = Data(lines, data.level + 1)

        findRelation(sub_data, arr[1:], to)
        
    else:
        for line in data.at_level(0).lines:
            text = line.name(0) + '.' + line.name(1)
            entry = splitName(line.text())
            entry = entry[:2] + [to]
            DATA.lines.append(Line(entry, text, '<through-relation>', 0))

################################################################################

def clearData():
    global DATA
    DATA = Data([], 0)
    global ID
    ID = 0

################################################################################

def loadDDDOCFile(filename, cache):
    if cache.isFresh(filename):
        return cache.get(filename)
    
    f = open(filename)
    text = f.readlines()
    f.close()

    cache.set(filename, text)
    return text

################################################################################

def loadCPPFile(filename, cache):
    if cache.isFresh(filename):
        return cache.get(filename)
    
    f = open(filename)
    lines = f.readlines()
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

################################################################################

def getFileType(filename):
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


################################################################################
    
def loadFile(filename, cache):
    file_type = getFileType(filename)
    if (file_type == 2): return loadCPPFile(filename, cache)
    elif (file_type == 1): return loadDDDOCFile(filename, cache)
    else: raise "unknown file type"


################################################################################

def parseFile(filename, cache):
    text = loadFile(filename, cache)
    text.append('.')
    
    context = [[]]
    str = False

    line_no = 0
    for line in text:
        line_no += 1
        if line != '': 
            if line[0] == '.':
                parseString(str, context, filename, line_no)
                str = line
            elif str:
                if str[len(str)-1] != '\n': str += '\n'
                str += line
                
################################################################################
                
def parseString(str, context, file_name, line_no):
    global DATA
    
    if not str or (str == '.'):
        return [[]]
        
    level = 0
    while (level < len(str)) and (str[level] == '.'): 
        level += 1
        
    str = str[level:]
        
    if (level < len(context)):
        del context[level:]
        
    if len(context) > 0:
        entry = copy.copy(context[len(context) - 1])
    else:
        entry = []

    key = ''
    text = ''
    
    pos = 0
    is_escaped = False
    c_quoted = ''
    while (pos < len(str)):
        c = str[pos]
        if c == "\x0d":
            pos += 1
            continue
        if c_quoted != "":
            if c_quoted == c: c_quoted = ""
            else: key += c                       
        elif is_escaped: 
            key += c
            is_escaped = False
        else:
            if c == '\\': is_escaped = True
            elif c in ['"', "'"]: c_quoted = c
            elif (c == ':'):
                key = str[0:pos]
                text = str[pos+1:]
                break
            else: key += c
                
        pos += 1    
        
    entry += splitName(key)
    DATA.lines.append(Line(entry, text, file_name, line_no))
    context.append(entry)


################################################################################

def splitName(line):
    pos = 0
    key = ""
    c_quoted = ""
    is_escaped = False
    li = []
    while (pos < len(line)):
        c = line[pos]
        if c_quoted != "":
            if c_quoted == c: c_quoted = ""
            else: key += c                       
        elif is_escaped: 
            key += c
            is_escaped = False
        else:
            if c == '\\': is_escaped = True
            elif c in ['"', "'"]: c_quoted = c
            elif c == '.':
                if key != "":
                    li.append(key)
                    key = ""
            elif c == '|':
                if key != "":
                    li.append(key)
                    key = ""
                rest = line[pos+1:]
                if len(rest)>0: li.append(rest)
                break;
            elif c != '\n':
                key += c
                
        pos += 1    
    
    if key != "": li.append(key)
    
    return li


def splitUrl(line):
    """Splits a tuple at separator characters.

    The separator character is '|'.  These characters can be escaped
    using the backslash sign '\', entries can also be quoted.

    >>> splitUrl('a|b|c')
    ['a', 'b', 'c']
    >>> splitUrl('a\|b|c')
    ['a|b', 'c']
    >>> splitUrl('"a|b"|c')
    ['a|b', 'c']

    Args:
      line  String to split.

    Returns
      List with strings, split at | symbols, excluding these symbols
      themselves.
    """
    pos = 0
    key = ""
    c_quoted = ""
    is_escaped = False
    li = []
    while (pos < len(line)):
        c = line[pos]
        if c_quoted != "":
            if c_quoted == c: c_quoted = ""
            else: key += c                       
        elif is_escaped: 
            key += c
            is_escaped = False
        else:
            if c == '\\': is_escaped = True
            elif c in ['"', "'"]: c_quoted = c
            elif c == '|':
                if key != "":
                    li.append(key)
                    key = ""
            elif c != '\n':
                key += c
                
        pos += 1    
    
    if key != "": li.append(key)

    return li


def loadFiles(search_path, cache):
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
            if getFileType(path) in [FILETYPE_CPP, FILETYPE_DDDOC]:
                parseFile(path, cache)
        # Exclude ignored diretories.
        for ignored in IGNORED_DIRS:
            if ignored in dirs:
                dirs.remove(ignored)
