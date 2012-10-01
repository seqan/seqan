cmake_minimum_required (VERSION 2.6)
project (seqan_extras_tests)

seqan_setup_tests(extras_tests)  # Setup flags etc.

seqan_add_all_subdirectories()  # include_directory(*)
