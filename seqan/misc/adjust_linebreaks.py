# Dieses Tool durchsucht alle Unterverzeichnisse und baut alle linebreaks zu \r\n um
# Vorsicht! Katastrophen nicht ausgeschlossen! Backups machen!

import os
import copy
import string


def main(search_path):
    counter = 0
    for root, dirs, files in os.walk(search_path):
        for file in files:
            path = os.path.join(root, file)
            if testFileType(path):
                if parseFile(path): counter += 1
        if 'CVS' in dirs:
            dirs.remove('CVS')
    print counter, "lf corrupted files corrected"

def testFileType(filename):
    pos = filename.rfind(".")
    if (pos >= 0): ext = filename[pos+1:]
    else: ext = ""

    return ext in ["c", "C", "cpp", "CPP", "c++", "C++", "h", "H", "hpp", "HPP", "h++", "H++", "dddoc", "DDDOC", "txt"]

  
    
def parseFile(path):
    f = open(path, "rb")
    text = f.read()
    f.close()
    
    s = ""
    i = 0;
    corrupt = False
    found_d = False
    found_first_a = False
    
    while i < len(text):
        c = text[i]
        i += 1  
        
        if ord(c) == 0x0a:
            if not found_d:
                corrupt = True
                if not found_first_a:
                    s += '\r\n'
                    found_first_a = True
            found_d = False      
            
        elif ord(c) == 0x0d:
            if found_d:
                corrupt = True
            else:
                found_d = True
                if not found_first_a:
                    s += '\r\n'
            found_first_a = False
        else:
            s += c
            found_d = False
            found_first_a = False
            
    if corrupt:
        f = open(path, "wb")
        f.write(s)
        f.close()
        print '\n', path, 
    else:
        print '.',

    return corrupt
        
main ('.')
