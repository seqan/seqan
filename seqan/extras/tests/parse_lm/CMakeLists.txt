cmake_minimum_required (VERSION 2.6)
project (sandbox_holtgrew_tests_parse_lm)

# Update the list of file names below if you add source files to your test.
seqan_add_test_executable(test_parse_lm test_parse_lm.cpp test_parse_lm.h)
