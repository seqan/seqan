cmake_minimum_required (VERSION 2.6)
project (seqan_core_tests)

seqan_setup_tests(core_tests)  # Setup flags etc.

seqan_add_all_subdirectories()  # include_directory(*)
