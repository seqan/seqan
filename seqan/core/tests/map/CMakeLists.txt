cmake_minimum_required (VERSION 2.6)
project (core_tests_map)

# Update the list of file names below if you add source files to your test.
seqan_add_test_executable(test_map test_map.cpp
				   test_map_map.h
				   test_map_sumlist.h)
