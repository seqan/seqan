#!/usr/bin/env python

import sys
import os
import webbrowser
from bkahlert import DiffCollector

def main():
	# sys.argv[1]: event name
	# sys.argv[2]: cmake binary dir (e.g. [SEQAN]/build/Debug)
	# sys.argv[3]: seqan src directory
	diffCollector = DiffCollector(sys.argv[2], sys.argv[3])
	id = diffCollector.getID().get()
	url = "https://dalak.imp.fu-berlin.de/SUAsrv/static/register.html?id=" + id
        
	# cmake run
	if(sys.argv[1] == "cmake"):
		diffCollector.prepare()
		print("[NOTE] Your ID is " + id)
		print("[NOTE] Your documented data are saved in file\n       "
			+ diffCollector.getStatsFile() + ".")
		
		# register ID if never happened before
		if(diffCollector.getID().isLinked()):
			#print("[NOTE] Your ID has been linked on "
			#	+ diffCollector.getID().getLinkedDate().strftime("%Y-%m-%d %H:%M:%S")
			#	+ ". You can always refresh your link by opening\n       " + url)
			sys.path
		else:
			if(not webbrowser.open(url=url, new=2)):
				print("\n\n")
				print("  Please register your ID:")
				print("  ------------------------")
				print("  1) Open your favorite browser")
				print("  2) Enter the following website:\n")
				print("     >>>   " + url + "   <<<\n")
				print("  3) Press ENTER when you have finished to continue ...")
				print("\n\n")
				raw_input("")
			diffCollector.getID().link()
		
		diffCollector.build()
		return 0

        # build (only called no matter how many target need to be rebuild)
	if(sys.argv[1] == "build"):
		diffCollector.prepare()
		diffCollector.build()
		return 0

	# pre build called before each target build
	# sys.argv[4]: target name
	if(sys.argv[1] == "pre_build"):
		return 0

	# post build called after each target build
	# sys.argv[4]: target name
	if(sys.argv[1] == "post_build"):
		return 0

	return 0

if __name__ == '__main__':
	sys.exit(main())

