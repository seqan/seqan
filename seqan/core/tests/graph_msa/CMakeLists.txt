cmake_minimum_required (VERSION 2.6)
project (core_tests_graph_msa)

# Update the list of file names below if you add source files to your test.
seqan_add_test_executable(test_graph_msa test_graph_msa.cpp
                                         test_graph_msa_guide_tree.h
                                         test_graph_tcoffee.h)
