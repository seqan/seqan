cmake_minimum_required (VERSION 2.6)
project (extras_tests_statistics)

# Update the list of file names below if you add source files to your test.
seqan_add_test_executable(test_statistics test_statistics.cpp test_statistics_base.h test_statistics_markov_model.h)
