import dddoc
import dddoc_tex
import sys

# Command line: main.py <inpath> [ <module> ]*

if len(sys.argv) < 2: inpath = "..\\..\\projects\\library"
else: inpath =  sys.argv[1]

print("Welcome to Dot.Dot.Doc")

if len(sys.argv) < 3:
    print("Scanning " + inpath + "...")
    dddoc.loadFiles(inpath)
else:
    i = 2;
    while (i < len(sys.argv)):
        modulepath = inpath + "\\seqan\\" + sys.argv[i]
        print("Scanning " + modulepath + "...")
        dddoc.loadFiles(modulepath)
        i += 1


print("Scanning .\\pages...")
dddoc.loadFiles("pages")

print("Scanning .\\concepts...")
dddoc.loadFiles("concepts")


dddoc.DATA.init()

print("Create LaTeX Documentation...")
dddoc_tex.createDocs("tex")

print("Documentation created.")

#raw_input("press return")