cmake_minimum_required (VERSION 2.6)
project (core_tests_sequence_journaled)

# Update the list of file names below if you add source files to your test.
seqan_add_test_executable(test_sequence_journaled test_sequence_journaled.cpp
						  test_sequence_journaled.h
						  test_sequence_journaled_iterator.h)
