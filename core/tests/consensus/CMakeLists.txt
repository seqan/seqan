cmake_minimum_required (VERSION 2.6)
project (core_tests_consensus)

# Update the list of file names below if you add source files to your test.
seqan_add_test_executable(test_consensus test_consensus.cpp)
