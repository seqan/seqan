#!/usr/bin/env python2.5
"""SeqAn Automatic Forwards Generator.

Usage: build_forwards.py BASE_PATH [TARGET_PATH] [all]

The program is given the project base path.  All directories in this
base path contain one module which consists of all header files
(ending in ".h") in this module directory.

For each module, a module generated forwards header can be created.
The file name is the name of the module, followed by
"_generated_forwards.h".  If the module generated forward header does
not exist or is older than any file in the module, it is recreated.

Original Author: Andreas Gogol Doering <andreas.doering@mdc-berlin.de>
"""

__author__ = 'Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>'

import os
import os.path
import copy
import string
import sys


FUNCS = {}
CLASSES = {}
TYPEDEFS = {}

PROGRAM_USAGE = """
SeqAn Automatic Forwards Generator

USAGE: build_forwards.py BASE_PATH [all]

BASE_PATH is the path up to and including "seqan".  The option "all"
forces a rebuild.
""".strip()


def buildProject(project_path, target_path):
    if not os.path.exists(project_path):
        return

    print "create forwards for", project_path

    global FUNCS
    FUNCS = {}
    
    global CLASSES
    CLASSES = {}
    
    global TYPEDEFS
    TYPEDEFS = {}
    
    pos1 = project_path.rfind('/')
    if (pos1 < 0):
        exit ('ERROR: wrong argument "' + project_path + '"');
        
    project = project_path[pos1+1:]
    
    for root, dirs, files in os.walk(project_path):
        if 'CVS' in dirs:
            dirs.remove('CVS')
        if '.svn' in dirs:
            dirs.remove('.svn')
        for dir in dirs:
          if dir.startswith('.'):
            dirs.remove(dir)
        for file in files:
            if file.startswith('.'):
              continue  # Skip hidden files.
            if file == forwardFilename(project):
                continue
            path = os.path.join(root, file)
            if testFileType(path):
                parseFile(path)
                
#    if FUNCS != {}:        
    outAll(target_path, project)
        
    print


def forwardFilename(module):
    """Returns the generated forwards header filename for a module name."""
    return module + "_generated_forwards.h"


def testFileType(filename):
    pos = filename.rfind(".")
    if (pos >= 0): ext = filename[pos+1:]
    else: ext = ""

    return ext in ["c", "C", "cpp", "CPP", "c++", "C++", "h", "H", "hpp", "HPP", "h++", "H++"]
  
    
def parseFile(filename):
    f = open(filename)
    lines = f.readlines()
    f.close()
    
    for line in lines:
        if (line.find("SEQAN_NO_GENERATED_FORWARDS") >= 0):
            print "-",
            return
    
    print ".",

    sigs = preprocess(lines, filename);
    
    createEntries(sigs)


def isMultiLine(line):
    """Returns true iff line ends with '\'."""
    if (len(line) == 0): return False
    return line[len(line)-1] == '\\'
    

def preprocess(lines, filename):
    """Removes comments, linebreaks, and macro definitions returns a string."""
    ret = []
    
    inComment = False
    inString = False
    inDefine = False
    curlyCount = 0
    namespaces = []
    lineNumber = 0
    
    str = ""
    
    for line in lines:
        lineNumber += 1;

        line = line.strip()
    
        #remove some characters to make parsing simpler
        line = line.replace("'}'", "");
        line = line.replace("'{'", "");
        line = line.replace("';'", "");
        line = line.replace("'\"'", "");
        line = line.replace('\\"', "");
        

        #skip multiline defines
        if inDefine or ((not inComment) and (not inString) and (line.find("#define") >= 0)):
            inDefine = isMultiLine(line)
            continue
            
        if len(line) <= 0: 
            continue

        #skip preprocessor lines
        if (not inComment) and (not inString) and (line[0]=="#"):
            continue

        while len(line) > 0 :
            if inComment:
                pos1 = line.find("*/")
                if pos1 >= 0:
                    line = line[pos1 + 2:]
                    inComment = False
                else:
                    break
                    
            elif inString:
                pos1 = line.find('"')
                if pos1 >= 0:
                    line = line[pos1 + 1:]
                    inString = False
                else: 
                    break
                    
            else:
                pos1 = line.find("/*")
                pos2 = line.find("//")
                pos3 = line.find('"')
                pos4 = line.find('{')
                pos5 = line.find('}')
                pos6 = line.find(';')
                
                if (pos1 >= 0) and ((pos2 < 0) or (pos1 < pos2)) and ((pos3 < 0) or (pos1 < pos3)) and ((pos4 < 0) or (pos1 < pos4)) and ((pos5 < 0) or (pos1 < pos5)) and ((pos6 < 0) or (pos1 < pos6)):
                    if curlyCount == 0:
                        str += ' ' + line[0:pos1].strip()
                    line = line[pos1 + 2:]
                    inComment = True
                    
                elif (pos2 >= 0) and ((pos3 < 0) or (pos2 < pos3)) and ((pos4 < 0) or (pos2 < pos4)) and ((pos5 < 0) or (pos2 < pos5)) and ((pos6 < 0) or (pos2 < pos6)):
                    if curlyCount == 0:
                        str += ' ' + line[0:pos2].strip()
                    break
                    
                elif (pos3 >= 0) and ((pos4 < 0) or (pos3 < pos4)) and ((pos5 < 0) or (pos3 < pos5)) and ((pos6 < 0) or (pos3 < pos6)):
                    if curlyCount == 0:
                        str += ' ' + line[0:pos3].strip()
                    line = line[pos3 + 1:]
                    inString = True
                    
                elif (pos4 >= 0) and ((pos5 < 0) or (pos4 < pos5)) and ((pos6 < 0) or (pos4 < pos6)):
                    if curlyCount == 0:
                        entry = process2(str + ' ' + line[:pos4])
                        nam = isNamespace(entry)
                        if nam != "":
                            namespaces += [nam]
                        else:             
                            ret = ret + [[filename, lineNumber, entry, namespaces + []]]
                            curlyCount = 1
                        str = ""
                    else:
                        curlyCount += 1

                    line = line[pos4 + 1:]
                    
                elif (pos5 >= 0) and ((pos6 < 0) or (pos5 < pos6)):
                    line = line[pos5 + 1:]
                    if curlyCount > 0:
                        curlyCount -= 1
                    elif len(namespaces) > 0:
                        namespaces = namespaces[:len(namespaces)-1]
                    else:
                        print "ERROR in" , filename , "(", lineNumber, "): Too many }"
                        
                elif (pos6 >= 0):
                    if curlyCount == 0:
                        entry = process2(str + ' ' + line[:pos6])
                        if isStuctDeclaration(entry) or isTypedef(entry):
                            ret = ret + [[filename, lineNumber, entry, namespaces + []]]
                            
                    str = ""
                    line = line[pos6 + 1:]
                        
                else:
                    if curlyCount == 0:
                        str += ' ' + line
                    break
    
    return ret


def process2(str):
    """Shrinks to the interesing part of the signature."""
    str = str.replace('\t', ' ')
    str = str.replace('  ', ' ')
    
    pos1 = str.rfind(';')
    if pos1 >= 0:
        str = str[pos1+1:]
        
    str = removeBases(str)
        
    str = str.replace('template<', 'template <')
    str = str.replace('template < ', 'template <')
    str = str.replace(' operator ', ' operator')
        
    str = str.strip()
    return str;


def isNamespace(str):
    """Determines whether a signature is a namespace.

    Returns the namespace name or "" if it is no namespace
    """
    pos1 = str.rfind(' ')
    if pos1 < 0: return ""
    while (pos1 > 0) and (str[pos1] == ' '): pos1 -= 1
    if (pos1 >= 8) and (str[pos1-8:pos1+1] == 'namespace'):
        return str[pos1 + 1:].strip()
    else:
        return ""


def isStuctDeclaration(str):
    """Determines whether a signature is a struct declaration."""
    str = removeBases(str)
    
    pos1 = str.rfind(' ')
    if pos1 < 0: return False
    while (pos1 > 0) and (str[pos1] == ' '): pos1 -= 1
    return ((pos1 >= 5) and (str[pos1-5:pos1+1] == 'struct')) or ((pos1 >= 4) and (str[pos1-4:pos1+1] == 'class'))
   

def isTypedef(str):
    """Determines whether a signature is typedef."""
    str = str.strip()
    return (str[:7] == "typedef") and (str.find('::') < 0)


def removeBases(str):
    """Removes list of base class."""
    pos1 = -2
    while pos1 < len(str):
        pos1 = str.find(':', pos1 + 2)
        if pos1 < 0: return str
        if (pos1 == len(str)-1) or (str[pos1+1]!=':'): return str[:pos1]
        
    return str.strip()
    

def createEntries(sigs):
    """Gets a list of [filename, lineNumber, signature, namespaces] tupels.

    Analyse signature and adds entries in FUNCS and CLASSES.
    """
    for data in sigs:
        filename = data[0]
        lineNumber = data[1]
        sig = data[2]
        namespaces = data[3]
        
        entry = makeEntry(filename, lineNumber, sig)

        name = getTypedefName(sig)
        if name != '':
            addEntry(TYPEDEFS, name, entry, namespaces)
        else:
            name = getStructName(sig)
            if name != '':
                addEntry(CLASSES, name, deleteDefaultArguments(entry, '<'), namespaces)
            else:
                name = getFuncName(sig)
                if name != '':
                    addEntry(FUNCS, name, deleteDefaultArguments(entry, '('), namespaces)


def deleteDefaultArguments(str, start_delim):
    """Deletes all default arguments from argument lists.
    
    Use delim = '>' for template argument lists and ')' for function argument lists.
    """
    ret = ""

    start = str.find(start_delim);
    if start >= 0: 
        ret = str[:start+1]
        str = str[start+1:]
        
    str = str.replace("<<", "")
    str = str.replace(">>", "")

    while str != "":
        pos1 = findCharOutsideBrackets(str, 0, "=")
        if (pos1 < 0) or (str[pos1] != "="):
            ret += str
            break
        ret += str[:pos1]
        pos2 = findCharOutsideBrackets(str, pos1, ",")
        if pos2 >= 0:
            str = str[pos2:]
        else:
            print "ERROR while deleting default arguments"
            break
            
    return ret


def findCharOutsideBrackets(str, start_pos, char, verbose = False):
    """Returns position of the first occurence of "char", or the position of the
    first closing bracket, whatever comes first.
    
    Areas in brackets are ignored.
    """
    pos = start_pos
    edge_count = 0
    while pos < len(str):
        if verbose: 
            print pos, edge_count
            print str[pos:]
        
        p1 = str.find(char, pos)
        p2 = str.find("<", pos)
        p2a = str.find("(", pos)
        p3 = str.find(">", pos)
        p3a = str.find(")", pos)
        
        if (p2 < 0) or ((p2a >= 0) and (p2a < p2)): p2 = p2a
         
        if (p3 < 0) or ((p3a >= 0) and (p3a < p3)): p3 = p3a

        if (p1 >= 0) and ((p2 < 0) or (p1 < p2)) and ((p3 < 0) or (p1 < p3)):
            if edge_count == 0:
                return p1
            else:
                pos = p1 + 1
            
        elif (p2 >= 0) and ((p3 < 0) or (p2 < p3)):
            edge_count += 1
            pos = p2 + 1
            
        elif (p3 >= 0):
            if edge_count == 0:
                return p3
            else:
                edge_count -= 1
                pos = p3 + 1
            
        else:
            return -1


def makeEntry(filename, lineNumber, sig):
    """Returns the string that is inserted into the header."""
    text = sig + ";       \t// \"" + filename + "\"(" + str(lineNumber) + ")"
    return text


def getSortKey(name, namespaces):
    """Returns the key the functions and structs are sorted for."""
    return str(namespaces) + name


def addEntry(arr, name, entry, namespaces):
    """Adds a signature to FUNCS or CLASSES."""
    key = getSortKey(name, namespaces)
    if not arr.has_key(key):
        arr[key] = []
    arr[key] += [[name, entry, namespaces]]
    

def getFuncName(sig):
    """Get the function name from a signature or '', if signature is not a function."""
    sig = sig.strip()
    pos1 = sig.rfind('(')
    if pos1 < 0:
        return ""
    else:
        pos1 -= 1
        while (pos1 >= 0) and (sig[pos1] == ' '):
            pos1 -= 1
            
        pos2 = sig.rfind(' ', 0, pos1)
        return sig[pos2 + 1: pos1 + 1]

def getStructName(sig):
    """Get the class name from a signature or '', if signature is not a class."""
    sig = sig.strip()

    pos1 = sig.rfind(' ')
    if pos1 < 0: return ""
    
    while (pos1 > 0) and (sig[pos1] == ' '): pos1 -= 1
    
    if ((pos1 >= 5) and (sig[pos1-5:pos1+1] == 'struct')) or ((pos1 >= 4) and (sig[pos1-4:pos1+1] == 'class')):
        name = sig[pos1 + 1:].strip()
        if (name.find('<') >= 0): return ""
        else: return name
       
    return ""


def getTypedefName(sig):
    """Get the typedef name from a signature or '', if signature is not a typedef."""
    sig = sig.strip()
    if not isTypedef(sig): return ""
    pos1 = sig.rfind(' ')
    if pos1 < 0: return ""
    return sig[pos1+1:].strip()



def outAll(path, project):
    """Main Function for output of forward header."""
    
    header_switch = "SEQAN_HEADER_" + project + "_GENERATED_FORWARDS_H"

    str = ""
    str += '// ==========================================================================\n'
    str += '//                 SeqAn - The Library for Sequence Analysis\n'
    str += '// ==========================================================================\n'
    str += '// Copyright (c) 2006-2010, Knut Reinert, FU Berlin\n'
    str += '// All rights reserved.\n'
    str += '//\n'
    str += '// Redistribution and use in source and binary forms, with or without\n'
    str += '// modification, are permitted provided that the following conditions are met:\n'
    str += '//\n'
    str += '//     * Redistributions of source code must retain the above copyright\n'
    str += '//       notice, this list of conditions and the following disclaimer.\n'
    str += '//     * Redistributions in binary form must reproduce the above copyright\n'
    str += '//       notice, this list of conditions and the following disclaimer in the\n'
    str += '//       documentation and/or other materials provided with the distribution.\n'
    str += '//     * Neither the name of Knut Reinert or the FU Berlin nor the names of\n'
    str += '//       its contributors may be used to endorse or promote products derived\n'
    str += '//       from this software without specific prior written permission.\n'
    str += '//\n'
    str += '// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"\n'
    str += '// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE\n'
    str += '// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE\n'
    str += '// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE\n'
    str += '// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL\n'
    str += '// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR\n'
    str += '// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER\n'
    str += '// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT\n'
    str += '// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY\n'
    str += '// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH\n'
    str += '// DAMAGE.\n'
    str += '//\n'
    str += '// ==========================================================================\n'
    str += '\n\n'

    str += "#ifndef " + header_switch.upper() + " \n"
    str += "#define " + header_switch.upper() + " \n\n"
    
    str += "//////////////////////////////////////////////////////////////////////////////\n"
    str += "// NOTE: This file is automatically generated by build_forwards.py\n"
    str += "//       Do not edit this file manually!\n"
    str += "//////////////////////////////////////////////////////////////////////////////\n\n\n"

    str += "//////////////////////////////////////////////////////////////////////////////\n"
    str += "// CLASSES\n"
    str += outList(CLASSES)
    
    str += "\n//////////////////////////////////////////////////////////////////////////////\n"
    str += "// TYPEDEFS\n"
    str += outList(TYPEDEFS)

    str += "\n//////////////////////////////////////////////////////////////////////////////\n"
    str += "// FUNCTIONS\n"
    str += outList(FUNCS)
    
    str += "#endif\n"

    if not os.path.exists(path):
        os.makedirs(path)

    filename = os.path.join(path, forwardFilename(project))
    fl = file(filename, "w")
    fl.write(str + "\n")
    fl.close()

def outList(lst):
    keys = lst.keys()
    keys.sort()
    
    namespaces = []
    str = ""
    
    for key in keys:
        first_entry = lst[key][0]
        new_namespaces = first_entry[2]
        str += outChangeNamespaces(namespaces, new_namespaces)
        namespaces = new_namespaces
            
        str += "//____________________________________________________________________________\n"
        str += "// " + first_entry[0] + "\n\n"
        
        for entry in lst[key]:
            str += entry[1] + "\n"
            
    str += outChangeNamespaces(namespaces, [])
    
    return str
    

def outChangeNamespaces(old_namespaces, new_namespaces):
    """If old_namespaces != new_namespaces, close old namespace and open new one."""
    str = ""
    if old_namespaces != new_namespaces:
        if len(old_namespaces) > 0: 
            str += "\n"
        
        while len(old_namespaces) > 0:
            str += "} //namespace " + old_namespaces[len(old_namespaces)-1] + "\n"
            old_namespaces = old_namespaces[:len(old_namespaces)-1]
            
        if len(new_namespaces) > 0:
            str += "//////////////////////////////////////////////////////////////////////////////\n\n"
            
        while len(new_namespaces) > 0:
            str += "namespace " + new_namespaces[0] + " {\n"
            new_namespaces = new_namespaces[1:]
            
    return str + "\n"


def forwardsNeedsRebuild(module_path, forwards_filename):
    """Determine whether the given forwards header needs rebuild.

    Args:
      module_path        String, path to the module.
      forwards_filename  String, path to the forwards header.
    """
    forwards_mtime = os.path.getmtime(forwards_filename)
    for f in os.listdir(module_path):
        if f.startswith('.'):
            continue  # Ignore hidden files.
        f_mtime = os.path.getmtime(os.path.join(module_path,f))
        if f_mtime > forwards_mtime:
            return True  # There is a newer file!
    return False


def buildAllForwards(project_path, target_path, force_recreate=False):
    """Build forwards for all modules below project_path.

    Searches for projects without generated forward and build forward
    file if force_recreate is True, all forward files are rebuilt
    returns True if one or more files were built.
    """
    # Strip trailing slash from project_path.
    if project_path[-1] == '/':
        project_path = project_path[:-1]
    if target_path[-1] == '/':
        target_path = target_path[:-1]
    
    ret = False
    
    for f in os.listdir(project_path):
        if f in ['CVS', '.svn']:
            continue  # Skip CVS and .svn directories.
        if f.startswith('.'):
            continue  # Skip hidden files.
        # All directories below the project are module directories.
        module_name = f
        module_path = project_path + "/" + module_name
        if os.path.isdir(module_path):
            file = target_path + "/" + module_path + "/" + forwardFilename(module_name)
            if force_recreate or not os.path.exists(file) or \
                    forwardsNeedsRebuild(module_path, file):
                ret = True
                buildProject(module_path, target_path + "/" + module_name)
    return ret


def main():
  """Main entry point for the forwards generator."""
  if len(sys.argv) < 2: 
    print >>sys.stderr, 'ERROR: Too few arguments.'
    print >>sys.stderr, PROGRAM_USAGE
    return 1

  force_rebuild = sys.argv[-1] == 'all'
  if force_rebuild:
      target_path = sys.argv[-2]
  else:
      target_path = sys.argv[-1]

  buildAllForwards(sys.argv[1], target_path, force_rebuild)

  return 0

if __name__ == '__main__':
  sys.exit(main())
