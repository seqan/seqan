import os
import subprocess
import glob
from datetime import datetime
import socket
import requests

class Diff(object):
	MAX_REVISION_LENGTH = 8
        
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

	def get_diff_filename(self, revision, datetime):
		return self.results_dir + "/" + self.id.get() + "_r" + str.rjust(str(revision), self.MAX_REVISION_LENGTH, "0") + "_" + datetime.strftime("%Y-%m-%dT%H-%M-%S") + self.get_timezone_string() + ".diff"
	
	def get_diff_filenames(self):
		return glob.glob(self.results_dir + "/" + self.id.get() + "_r*_*.diff")

	def get_next_diff_filename(self):
		r = len(self.get_diff_filenames())
		if(r >= 10**self.MAX_REVISION_LENGTH):
			raise Exception("Maximum number of " + (10**self.MAX_REVISION_LENGTH) + " exceeded!")
		return self.get_diff_filename(r, datetime.today())

	def save_diff_reset(self, filename):
		f = open(filename, "w")
		f.write("RESET\n")
		f.close()

	def save_diff(self, filename):
		if(os.name == "nt"):
			p = subprocess.Popen([self.bin_dir + "/diffutils/bin/diff.exe", "-u", "-r", "-N"] + self.excluded_resources.get_exclude_list_diff() + ["." + self.from_dir_local, self.to_dir_local], stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=self.root_dir)
		else:
			p = subprocess.Popen(["diff", "-u", "-r", "-N"] + self.excluded_resources.get_exclude_list_diff() + ["." + self.from_dir_local, self.to_dir_local], stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=self.root_dir)	
		f = open(filename, "wb")
		output = p.communicate()[0]
		f.write(output)
		f.close()

	def sync(self, url):
		# get remote file list
		remote_files = []
		
		try:
			r = requests.get(url + "/" + self.id.get())
			if(r.status_code != 200):
				raise Exception;
			remote_files = {}
			for line in r.text.split("\n"):
				file = line.strip().split("\t")
				if(len(file) < 2): continue
				remote_files[file[0]] = file[1]
		except Exception:
			print("Error connecting to " + url + "! Retrying next time.")
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
