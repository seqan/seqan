import os
import subprocess
import glob
from datetime import datetime
import socket
import re # regular expressions
import zipfile
import hashlib

class Diff(object):
    def __init__(self, bin_dir, results_dir, root_dir, from_dir_local, to_dir_local, id, excluded_resources):
        self.results_dir = results_dir
        self.root_dir = root_dir
        self.from_dir_local = from_dir_local
        self.to_dir_local = to_dir_local
        self.id = id
        self.excluded_resources = excluded_resources
        self.bin_dir = bin_dir
        
    def get_timezone_string(self):
        import time
        offset = -time.timezone
        if(time.localtime().tm_isdst):
            offset = -time.altzone
        hours = abs(offset)/3600
        minutes = (abs(offset)%3600)/60
        return ("+" if offset >= 0 else "-") + str(hours).rjust(2, "0") + str(minutes).rjust(2, "0")

    def get_diff_filename(self, datetime):
        root_dir_hash = hashlib.md5(self.root_dir).hexdigest()[:4]
        return self.results_dir + "/" + self.id.get() + "_" + root_dir_hash + "_" + datetime.strftime("%Y-%m-%dT%H-%M-%S.%f") + self.get_timezone_string() + ".diff.zip"
    
    def save_diff(self, filename):
        # old parameters (to create a huge textual diff file: "-a", "-u", "-r", "-N"
        # change the language to C to avoid unparseable non-english output 
        # of localized shell tools
        my_env = os.environ.copy()
        my_env["LANG"] = "C"
        if(os.name == "nt"):
            p = subprocess.Popen([self.bin_dir + "/diffutils/bin/diff.exe", "-q", "-r", "-N"] + self.excluded_resources.get_exclude_list_diff() + ["." + self.from_dir_local, self.to_dir_local], stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=self.root_dir, env=my_env)
        else:
            p = subprocess.Popen(["diff", "-q", "-r", "-N"] + self.excluded_resources.get_exclude_list_diff() + ["." + self.from_dir_local, self.to_dir_local], stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=self.root_dir, env=my_env)
        output = p.communicate()[0]
        
        zipFile = zipfile.ZipFile(filename, "w")
        
        for line in output.split("\n"):
            match = re.search(" and \\.(.*) differ", line)
            if(match and match.groups > 0):
                currentFile = self.root_dir + match.group(1);
                oldFile = self.root_dir + self.from_dir_local + match.group(1);
                if(os.path.isfile(currentFile)):
                    zipFile.write(currentFile, match.group(1), zipfile.ZIP_STORED)

                #else:
                    # currentFile does not exist anymore
                    # ignore, since the old revision was already saved
                    
        zipFile.close()

    def sync(self, url):
        import requests
        
        # get remote file list
        remote_files = []
        
        readFileListUrl = url + "/" + self.id.get()
        try:
            r = requests.get(readFileListUrl, timeout=5.0)
            if(r.status_code != 200):
                print("Connection to " + readFileListUrl + " failed with status code " + r.status_code + "! Retrying next time.")
                return
            remote_files = {}
            for line in r.text.split("\n"):
                file = line.strip().split("\t")
                if(len(file) < 2): continue
                remote_files[file[0]] = file[1]
        except requests.exceptions.ConnectionError as e:
            print("\n!!! Connection to " + readFileListUrl + " failed! Retrying next time.") 
            return       
        except requests.exceptions.Timeout as e:
            print("\n!!! Connection to " + readFileListUrl + " timed out! Retrying next time.")
            return
        except Exception as e:
            print(type(e))
            print("\n !!!Connection to " + readFileListUrl + " failed for an unknow reason! Please inform bjoern.kahlert@fu-berlin.de.")
            print(e)
            return

        # upload all files that are not present on the ftp server
        for filename in glob.glob(self.results_dir + "/" + self.id.get() + "_*"):
            copy = True
            if(os.path.basename(filename) in remote_files.keys()):
                remote_size = int(remote_files[os.path.basename(filename)])
                local_size = os.path.getsize(filename)
                if(remote_size == local_size): copy = False
            
            if(copy):
                files = {self.id.get(): open(filename, 'rb')}
                r = requests.post(url, files=files)
