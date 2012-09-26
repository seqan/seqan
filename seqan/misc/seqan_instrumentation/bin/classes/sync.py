import os
import subprocess

class Sync(object):
	MAX_REVISION_LENGTH = 16
        
	def __init__(self, bin_dir):
		self.bin_dir = bin_dir
	
	def make_comparison_copy(self, from_dir, to_dir, excluded_resources):
		if(os.name == "nt"):
			p = subprocess.Popen([self.bin_dir + "/robocopy.exe", from_dir, to_dir, "/w:1", "/r:1", "/MIR"] + excluded_resources.get_exclude_list_win32_robocopy(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		else:
			p = subprocess.Popen(["rsync", "-auv", "--delete", "--delete-excluded"] + excluded_resources.get_exclude_list_other_rsync() + [from_dir + "/", to_dir], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		p.communicate()
	
