cmake_minimum_required (VERSION 2.6)
project (core_tests_find)

# Update the list of file names below if you add source files to your test.
seqan_add_test_executable(test_find test_find.cpp
                                    test_find_hamming.h
                                    test_find_myers_banded.h)
