class Dirs(object):
	def __init__(self, root, excluded_folders, excluded_files):
		self.root = root
		self.excluded_folders = excluded_folders
		self.excluded_files = excluded_files

	def get_root(self):
                return self.root

	def get_abs_dir_path(self, rel_dir_path):
		abs_dir_path = self.root + "/" + rel_dir_path
		self.excluded_folders.append(rel_dir_path);
		self.excluded_folders.append(abs_dir_path);
		return abs_dir_path
	
	def get_abs_file_path(self, rel_file_path):
		abs_file_path = self.root + "/" + rel_file_path
		self.excluded_files.append(rel_file_path);
		self.excluded_files.append(abs_file_path);
		return abs_file_path

	

	def get_exclude_list_diff(self):
		exclude_list = []
		for excluded_file in self.excluded_files + self.excluded_folders:
			exclude_list.append("-x")
			exclude_list.append(excluded_file)
		return exclude_list

	def get_exclude_list_win32_robocopy(self):
		exclude_list = []
		exclude_list.append("/XF")
		for excluded_file in self.excluded_files:
			exclude_list.append(excluded_file)
		exclude_list.append("/XD")
		for excluded_folder in self.excluded_folders:
			exclude_list.append(excluded_folder)
		return exclude_list

	def get_exclude_list_other_rsync(self):
		exclude_list = []
		for exclude in self.excluded_files + self.excluded_folders:
			exclude_list.append("--exclude")
			exclude_list.append(exclude)
		return exclude_list
