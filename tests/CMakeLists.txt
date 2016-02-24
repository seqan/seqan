# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /tests/CMakeLists.txt
#
# CMakeLists.txt file for core tests.
# ===========================================================================

cmake_minimum_required (VERSION 3.0.0)
project (seqan_tests CXX)
message (STATUS "Configuring tests")

# Setup flags and register all subdirectories.
seqan_register_tests ()
