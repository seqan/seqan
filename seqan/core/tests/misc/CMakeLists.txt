cmake_minimum_required (VERSION 2.6)
project (core_tests_misc)

# Update the list of file names below if you add source files to your test.
seqan_add_test_executable(test_misc test_misc.cpp
                          test_misc_accumulators.h
                          test_misc_interval_tree.h
                          test_misc_bit_twiddling.h
                          test_misc_long_word.h
                          test_misc_edit_environment.h)

seqan_add_test_executable(test_misc_cmdparser 
                          test_misc_cmdparser.cpp)