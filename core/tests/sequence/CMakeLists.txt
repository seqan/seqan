cmake_minimum_required (VERSION 2.6)
project (core_tests_sequence)

# Update the list of file names below if you add source files to your test.
seqan_add_test_executable(test_sequence test_sequence.cpp
					test_segment.h
					test_sequence_std_adaptions.h
					test_string.h
					test_stringset.h)
seqan_add_test_executable(test_sequence_v2 test_sequence_v2.cpp
					test_sequence.h
					test_segment_beta.h)
