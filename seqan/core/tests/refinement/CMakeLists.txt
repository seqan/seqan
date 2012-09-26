cmake_minimum_required (VERSION 2.6)
project (core_tests_refinement)

# Update the list of file names below if you add source files to your test.
seqan_add_test_executable(test_refinement test_refinement.cpp
					  test_graph_impl_align.h
					  test_graph_match_refinement.h)
