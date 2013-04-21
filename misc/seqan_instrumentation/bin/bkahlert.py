#!/usr/bin/env python

import os
import sys
import threading
import binascii
import shutil
from datetime import datetime
from classes import diff
from classes import dirs
from classes import flushfile
from classes import id
from classes import stats
from classes import sync

sys.stdout = flushfile.Flushfile(sys.stdout)

class DiffCollector(object):
    def __init__(self, cMake_binary_dir, src_dir):
        self.cMake_binary_dir = cMake_binary_dir
        self.src_dir = src_dir        
        
        self.dirs = dirs.Dirs(self.src_dir, [ ".git", ".svn", "bin", "build", "util", "misc", "docs", "docs2", "extras", "core" ], [ "*.o", "Thumbs.db", ".DS_Store", "CMakeCache.txt" ])
        self.bin_dir = self.dirs.get_abs_dir_path("misc/seqan_instrumentation/bin")
        try:
            os.mkdir(self.dirs.get_abs_dir_path("misc/seqan_instrumentation/last_revision_copy"))
        except Exception:
            pass
        self.last_revision_dir = self.dirs.get_abs_dir_path("misc/seqan_instrumentation/last_revision_copy")
        try:
            os.mkdir(self.dirs.get_abs_dir_path("misc/seqan_instrumentation/userdata"))
        except Exception:
            pass
        self.userdata_dir = self.dirs.get_abs_dir_path("misc/seqan_instrumentation/userdata")

        # if possible save user ID in the user's home directory
        # otherwise use sub folder of SeqAn installation
        old_id_file = self.dirs.get_abs_file_path("misc/seqan_instrumentation/userdata/id.txt")
        self.id_file = os.getenv("HOME", os.getenv("USERPROFILE"))
        if(self.id_file is None):
            self.id_file = old_id_file
        else:
            self.id_file += "/.SUA.ID"
            # migrate old ID if no new one exists
            if(os.path.isfile(old_id_file) and not os.path.isfile(self.id_file)):
                shutil.move(old_id_file, self.id_file)
                
        self.id = id.ID(self.id_file)

        self.stats_file = self.dirs.get_abs_file_path("misc/seqan_instrumentation/userdata/" + self.id.get() + "_stats.txt")
        self.stats = stats.Stats(cMake_binary_dir, src_dir, self.stats_file)
        self.stats.save("id", self.id.get())
        
    def getID(self):
        return self.id

    def getStatsFile(self):
        return self.stats_file

    def getStats(self):
        return self.stats

    def create_diff(self, differ):
        diff_filename = differ.get_diff_filename(datetime.today())
        shortform = diff_filename[len(self.src_dir):]
        print("create..."),    
        if not os.path.isdir(self.last_revision_dir):
            os.makedirs(self.last_revision_dir)
            
        differ.save_diff(diff_filename)
        print("create finished;"),

    def upload_diff(self, differ):
        print("upload..."),
        differ.sync("http://dalak.imp.fu-berlin.de/SUAsrv/diff")
        print("upload finished;"),

    def update_comparision_copy(self, syncer):
        print("sync..."),
        syncer.make_comparison_copy(self.src_dir, self.last_revision_dir, self.dirs)
        print("sync finished;"),
    
    def prepare(self):
        if(not os.path.exists(self.last_revision_dir)):
            os.makedirs(self.last_revision_dir)
        if(not os.path.exists(self.userdata_dir)):
            os.makedirs(self.userdata_dir)
    
    def build(self):
        differ = diff.Diff(bin_dir=self.bin_dir, results_dir=self.userdata_dir, root_dir=self.src_dir, from_dir_local=self.last_revision_dir[len(self.src_dir):], to_dir_local=".", id=self.id, excluded_resources=self.dirs)
        syncer = sync.Sync(bin_dir=self.bin_dir)

        print("[SUA]"),
        t1 = threading.Thread(target=self.create_diff, args=(differ,))
        t1.start()
        t1.join()

        t2 = threading.Thread(target=self.upload_diff, args=(differ,))
        t2.start()

        t3 = threading.Thread(target=self.update_comparision_copy, args=(syncer,))
        t3.start()
    
        t2.join()
        t3.join()
        print("finished."),
