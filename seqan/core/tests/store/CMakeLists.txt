cmake_minimum_required (VERSION 2.6)
project (core_tests_store)

# Update the list of file names below if you add source files to your test.
seqan_add_test_executable(test_store test_store.cpp
				     test_store_io_bam.h
				     test_store_io.h)
