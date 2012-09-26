import os
from datetime import datetime
import random
import string

"""Generates a permanent ID."""
class ID(object):
	"""
	@param cache_file is the file the generated ID is stored in
	"""
	def __init__(self, cache_file):
		self.minlength = 8
		self.length = 16
		self.cache_file = cache_file
		self.is_linked_file = cache_file + ".LINKED"
		
	def isLinked(self):
		return os.path.isfile(self.is_linked_file)
	
	def link(self):
		f = open(self.is_linked_file, "w")
		f.write("")
		f.close
		return
	
	def getLinkedDate(self):
		return datetime.fromtimestamp(os.path.getmtime(self.is_linked_file))

	def get(self):
		if(os.path.isfile(self.cache_file)):
			f = open(self.cache_file, "r")
			id = f.readlines()[0].strip()
			f.close
			if(len(id) >= self.minlength): return id
		
		id = ''.join(random.choice(string.ascii_lowercase + string.digits) for x in range(self.length))

		self.create_dir()
		
		f = open(self.cache_file, "w")
		f.write(id)
		f.close
		return id

	def __repr__(self):
		return get()

	def create_dir(self):
		try:
			os.makedirs(os.path.dirname(self.cache_file))
		except:
			pass
