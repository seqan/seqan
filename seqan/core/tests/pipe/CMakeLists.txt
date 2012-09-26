cmake_minimum_required (VERSION 2.6)
project (core_tests_pipe)

# Update the list of file names below if you add source files to your test.
seqan_add_test_executable(test_pipe test_pipe.cpp
				    test_pipe.h)
