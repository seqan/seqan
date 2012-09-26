cmake_minimum_required (VERSION 2.6)
project (core_tests_graph_types)



# Update the list of file names below if you add source files to your test.
seqan_add_test_executable(test_graph_types test_graph_types.cpp
					   test_graph_derived.h
					   test_graph_iterators.h
					   test_graph_properties.h
					   test_graph_types.h
					   test_graph_utils.h
					   test_graph_basic.h)


