#!/usr/bin/env python2.5
"""SeqAn IO-Revision code tagger step1

Usage: build_forwards.py BASE_PATH [all]

This program parses seqan for io-related functions, classes etc
and generates in each folder a file that lists all files in the folder
and corresponding lines of class definition etc

Step2 is awn script that parses thhis programs output and modifies the files

Author: Hannes Hauswedell <hauswedell@mi.fu-berlin.de>

Based on the Automatic Forwards Generator by
Andreas Gogol Doering <andreas.doering@mdc-berlin.de>
"""

__author__ = 'Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>'

import os
import os.path
import copy
import string
import sys
import re
from collections import defaultdict
import shutil
#from functools import cmp_to_key


#FUNCS = {}
#CLASSES = {}
#TYPEDEFS = {}


FILES_COMPLETELY_IORELATED = [ r'.*/seqan/file/.*\.h', \
                               ".*parsing.h" ] # todo finish
#FILES_COMPLETELY_IORELATED = [ ] # todo finish

FUNCTION_KEYWORDS = [ ".*read.*", \
                      ".*write.*", \
                      "_stream.*", \
                      "_parse.*", \
                      ".*file.*", \
                      ".*File.*", \
                      "_is.*"] # todo finish
#FUNCTION_KEYWORDS = [ ".*open\(.*" ] # todo finish

#FILES_LINES = {}

FILTERED_SIGS = []

PROGRAM_USAGE = """
SeqAn IO-Revision code tagger step1

USAGE: iorev_tagger.py BASE_PATH

BASE_PATH is the path up to and including "seqan".
""".strip()


def buildProject(project_path):
    if not os.path.exists(project_path):
        return

    #print "parsing: ", project_path

    #global FUNCS
    #FUNCS = {}
    
    #global CLASSES
    #CLASSES = {}
    
    #global TYPEDEFS
    #TYPEDEFS = {}

    #global FILE_LINES
    #FILE_LINES = {}
    
    pos1 = project_path.rfind('/')
    if (pos1 < 0):
        exit ('ERROR: wrong argument "' + project_path + '"');
        
    project = project_path[pos1+1:]
    
    
    for root, dirs, files in os.walk(project_path):
        if 'CVS' in dirs:
            dirs.remove('CVS')
        if '.svn' in dirs:
            dirs.remove('.svn')
        for file in files:
            if file.startswith('.'):
              continue  # Skip hidden files.
            if file == forwardFilename(project):
                continue
            path = os.path.join(root, file)
            if testFileType(path):
                parseFile(path)

#    if FUNCS != {}:        
    #outAll(project_path, project)
        



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
    
    #for line in lines:
        #if (line.find("SEQAN_NO_GENERATED_FORWARDS") >= 0):
            #print "-",
            #return
    
#    print ".",

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
    inClass = 0
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

                ## attempt to traverse classes
                #elif (pos4 >= 0) and ((pos5 < 0) or (pos4 < pos5)) and ((pos6 < 0) or (pos4 < pos6)):
                    #if curlyCount == 0:
                        #entry = process2(str + ' ' + line[:pos4])
                        #nam = isNamespace(entry)
                        #if nam != "":
                            #namespaces += [nam]
                        #else:             
                            #ret = ret + [[filename, lineNumber, entry, namespaces + []]]
                            #if isStuctDeclaration(entry):
                                #inClass += 1
                            #else:
                                #curlyCount = 1
                                
                        #str = ""
                    #else:
                        #curlyCount += 1

                    #line = line[pos4 + 1:]
                    
                #elif (pos5 >= 0) and ((pos6 < 0) or (pos5 < pos6)):
                    #line = line[pos5 + 1:]
                    #if curlyCount > 0:
                        #curlyCount -= 1
                    #elif len(namespaces) > 0:
                        #namespaces = namespaces[:len(namespaces)-1]
                    #elif inClass > 0:
                        #inClass -= 1
                    #else:
                        #print "ERROR in" , filename , "(", lineNumber, "): Too many }"

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
    

def comp_str(s1, s2):
    if (s1 == s2):
        return 0
    elif (s1 < s2 ):
        return -1
    else:
        return 1

#def filename_line_tupel_isBigger(tup1, tup2):
    
    #filediff  = comp_str(tup1[0], tup2[0])
    #if (filediff != 0 ):
        #return filediff
    #return int(tup1[1]) - int(tup2[1])
    
def matchesAny(str1, list_of_re):
    for reg in list_of_re:
        if reg.search(str1):
            return True
    return False

def compile_res(strings):
    res = []
    for s in strings:
        res.append(re.compile(s))
    return res
    
def createEntries(sigs):
    """Gets a list of [filename, lineNumber, signature, namespaces] tupels.

    Analyse signature and adds entries in FUNCS and CLASSES.
    """
    
    completefile_lines = defaultdict(list)
    noncompletefile_lines = defaultdict(list)

    #print FILES_COMPLETELY_IORELATED
    #print FUNCTION_KEYWORDS
    re1 = compile_res(FILES_COMPLETELY_IORELATED)
    re2 = compile_res(FUNCTION_KEYWORDS)

    
    for data in sigs: #sorted(sigs, key=lambda x:(x[0],int(x[1]))):
        filename = data[0]
        lineNumber = int(data[1])
        sig = data[2]
        namespaces = data[3]
        t = ''

        name = getTypedefName(sig)
        if name != '':
            t= "t"
        else:
            name = getStructName(sig)
            if name != '':
                t="s"
            else:
                name = getFuncName(sig)
                if name != '':
                    t="f"
                else:
                    continue
                    #print sig

        if matchesAny(filename, re1) or matchesAny(name, re2):
            FILTERED_SIGS.append((filename, lineNumber, t, name, namespaces))
            
        

        #print filename
        #print re1
        
        #if filename in completefile_lines.keys() or matchesAny(filename, re1):
            #completefile_lines[filename].append(int(lineNumber))
        #elif matchesAny(name, re2):
            #noncompletefile_lines[filename].append(int(lineNumber))

    ##print completefile_lines
    #FILES_LINES.update(completefile_lines)
    #FILES_LINES.update(noncompletefile_lines)
    

    

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

    #pos1 = sig.rfind(' ')
    #if pos1 < 0: return ""
    
    #while (pos1 > 0) and (sig[pos1] == ' '): pos1 -= 1
    
    #if ((pos1 >= 5) and (sig[pos1-5:pos1+1] == 'struct')) or ((pos1 >= 4) and (sig[pos1-4:pos1+1] == 'class')):
        #name = sig[pos1 + 1:].strip()
        #if (name.find('<') >= 0): return ""
        #else: return name
        
        
    pos1 = sig.find("struct ")
    while pos1 > 0 and sig[pos1-1] != ' ' and pos1 + len("struct ") < len(sig):
        pos1 = sig[pos1:].find("struct ")
    if pos1 < 0:
        pos1 = sig.find("class ")
        while pos1 > 0 and sig[pos1-1] != ' ' and pos1 + len("class ") < len(sig):
            pos1 = sig[pos1:].find("class ")
        if pos1 < 0:
            return ""
        pos1 = pos1 + len("class ")
    else:
        pos1 = pos1 + len("struct ")

    r = ""
    for c in sig[pos1:]:
        if c == ' ' or c == '<' or c == ';':
            break
        r += c
    return r


def getTypedefName(sig):
    """Get the typedef name from a signature or '', if signature is not a typedef."""
    sig = sig.strip()
    if not isTypedef(sig): return ""
    pos1 = sig.rfind(' ')
    if pos1 < 0: return ""
    return sig[pos1+1:].strip()


def printOutFilesAndLines():
    print FILES_LINES["projects/library/seqan/statistics/statistics_markov_model.h"]
    #TODO further processing




def addTags():
    lastfilename = ""
    lines_func_class = ""
    lines_typedef = ""
    for i in sorted(FILTERED_SIGS, key=lambda x:(x[0],x[1])):
        if i[0] != lastfilename:
            if (lastfilename != ""):
                cmd = 'awk -f ' + sys.path[0] + '/iorev_tagger_level2.awk'
                if (lines_func_class != ""):
                    cmd += " -v lines_fc=" + lines_func_class
                if (lines_typedef != ""):
                    cmd += " -v lines_t=" + lines_typedef
                cmd += " " + lastfilename + " > " + lastfilename + ".new"
                print cmd
                #if (lastfilename.find("file") >= 0):
                os.system(cmd)
                shutil.move(lastfilename, lastfilename + ".old")
                shutil.move(lastfilename+".new", lastfilename)
            lastfilename = i[0]
            lines_func_class = ""
            lines_typedef = ""
        if i[2] != "t":
            lines_func_class = lines_func_class + "_" + str(i[1])
        else:
            lines_typedef = lines_typedef + "_" + str(i[1])

    ## run for the last file aswell
    cmd = 'awk -f ' + sys.path[0] + '/iorev_tagger_level2.awk'
    if (lines_func_class != ""):
        cmd += " -v lines_fc=" + lines_func_class
    if (lines_typedef != ""):
        cmd += " -v lines_t=" + lines_typedef
    cmd += " " + lastfilename + " > " + lastfilename + ".new"
    print cmd
    #if (lastfilename.find("file") >= 0):
    os.system(cmd)
    shutil.move(lastfilename, lastfilename + ".old")
    shutil.move(lastfilename+".new", lastfilename)


def display(order="ntfl"):
    #TODO implement proper parsing of display and sort order

    if order == "fltn": # file line type name
        for i in sorted(FILTERED_SIGS, key=lambda x:(x[0],x[1],x[2],x[3])):
            print i[0], "\t", i[1], "\t", i[2], "\t", i[3]
    elif order == "ntfl": # name type file line
        for i in sorted(FILTERED_SIGS, key=lambda x:(x[3],x[2],x[0],x[1])):
            print i[3], "\t", i[2], "\t", i[0], "\t", i[1]
    elif order == "fl": # file and lines in one line
        lastfilename = ""
        for i in sorted(FILTERED_SIGS, key=lambda x:(x[0],x[1])):
                if i[0] != lastfilename:
                    lastfilename = i[0]
                    print "\n", i[0],
                print  " ", i[1],

        print

    

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


def buildAllForwards(project_path, force_recreate=False):
    """Build forwards for all modules below project_path.

    Searches for projects without generated forward and build forward
    file if force_recreate is True, all forward files are rebuilt
    returns True if one or more files were built.
    """
    # Strip trailing slash from project_path.
    if project_path[-1] == '/':
        project_path = project_path[:-1]
    
    ret = False
    
    for f in os.listdir(project_path):
        if (f == 'CVS') or (f == '.svn'):
            continue  # Skip CVS and .svn directories.
        if f.startswith('.'):
            continue  # Skip hidden files.
        # All directories below the project are module directories.
        module_name = f
        module_path = project_path + "/" + module_name
        if os.path.isdir(module_path):
            file = module_path + "/" + forwardFilename(module_name)
            if force_recreate or not os.path.exists(file) or \
                    forwardsNeedsRebuild(module_path, file):
                ret = True
                buildProject(module_path)
    return ret


def main():
  """Main entry point for the forwards generator."""
  if len(sys.argv) < 2: 
    print >>sys.stderr, 'ERROR: Too few arguments.'
    print >>sys.stderr, PROGRAM_USAGE
    return 1

  force_rebuild = 1

  buildAllForwards(sys.argv[1], force_rebuild)

  #display("fl")
  display()
  #addTags()

  return 0

if __name__ == '__main__':
  sys.exit(main())
