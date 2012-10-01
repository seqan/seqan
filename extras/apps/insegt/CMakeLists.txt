cmake_minimum_required (VERSION 2.6)
project (extras_apps_insegt)

# Update the list of file names below if you add source files to your application.
seqan_add_executable(insegt insegt.cpp
			    base.h
			    create_gff.h
			    fusion.h
			    overlap_module.h
			    read_gff.h
			    read_gtf.h)
