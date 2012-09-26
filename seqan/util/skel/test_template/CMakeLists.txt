cmake_minimum_required (VERSION 2.6)
project (%(CMAKE_PROJECT_NAME)s)

# Update the list of file names below if you add source files to your test.
seqan_add_test_executable(test_%(NAME)s test_%(NAME)s.cpp test_%(NAME)s.h)
