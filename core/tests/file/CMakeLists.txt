cmake_minimum_required (VERSION 2.6)
project (core_tests_file)

# Update the list of file names below if you add source files to your test.
seqan_add_test_executable(test_file test_embl.h
                                    test_file.cpp
                                    test_file.h)
