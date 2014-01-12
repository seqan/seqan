from distutils.core import setup
import sys
import py2exe
from glob import glob
import os

# ISSUE: import classes (via includes) instead of copying it (via data_files)

sys.path.append("C:\\Program Files\\CMake 2.8\\bin")
sys.path.append(os.path.normpath(os.getcwd() + "/../bin"))
setup(options={
    "py2exe": { "includes": ["platform", "glob", "shutil", "webbrowser", "cgi", "netrc", "urllib2", "Queue", "decimal", "os", "hashlib"], "skip_archive": True, "bundle_files": 3 } },
    data_files=[
        ("Microsoft.VC90.CRT", glob(r'.\Microsoft.VC90.CRT\*.*')),
        ("classes", glob(r'..\bin\classes\*.*')),
        ("classes/requests", glob(r'..\bin\classes\requests\*.*')),
        ("classes/requests/packages", glob(r'..\bin\classes\requests\packages\*.*')),
        ("classes/requests/packages/oreos", glob(r'..\bin\classes\requests\packages\oreos\*.*')),
        ("classes/requests/packages/urllib3", glob(r'..\bin\classes\requests\packages\urllib3\*.*')),
        ("classes/requests/packages/urllib3/packages", glob(r'..\bin\classes\requests\packages\urllib3\packages\*.*')),
        ("classes/requests/packages/urllib3/packages/mimetools_choose_boundary", glob(r'..\bin\classes\requests\packages\urllib3\packages\mimetools_choose_boundary\*.*')),
        ("classes/requests/packages/urllib3/packages/ssl_match_hostname", glob(r'..\bin\classes\requests\packages\urllib3\packages\ssl_match_hostname\*.*')),        
        ("classes/simplejson", glob(r'..\bin\classes\simplejson\*.*'))],
    console=["../bin/seqan_instrumentation.py"])
