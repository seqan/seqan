cmake_minimum_required (VERSION 2.6)
project (extras_tests_blast)

# Update the list of file names below if you add source files to your test.
seqan_add_test_executable(test_blast test_blast.cpp test_blast_parsing.h)
