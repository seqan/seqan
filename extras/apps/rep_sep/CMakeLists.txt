cmake_minimum_required (VERSION 2.6)
project (extras_apps_rep_sep)

# Update the list of file names below if you add source files to your application.
seqan_add_executable(rep_sep rep_sep.cpp
			     assembly_parser.h
			     column_scanner.h
			     rgraph.h
			     rgraph_base.h
			     rgraph_construction.h
			     rgraph_heuristics.h
			     rgraph_score.h
			     utils.h)
