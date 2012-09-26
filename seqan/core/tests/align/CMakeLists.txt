cmake_minimum_required (VERSION 2.6)
project (core_tests_align)

# Update the list of file names below if you add source files to your test.
seqan_add_test_executable(test_align test_align.cpp
                                     test_align_stream.h
                                     test_align_gaps.h
                                     test_align_gaps_iterator.h
                                     test_align_global_alignment.h
                                     test_align_global_alignment_banded.h
                                     test_align_global_alignment_specialized.h)
