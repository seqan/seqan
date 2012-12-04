# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /core/tests/misc/CMakeLists.txt
#
# CMakeLists.txt file for the misc module tests.
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (seqan_core_tests_misc)
message (STATUS "Configuring core/tests/misc")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search SeqAn and select dependencies.
set (SEQAN_FIND_DEPENDENCIES NONE)
find_package (SeqAn REQUIRED)

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Update the list of file names below if you add source files to your test.
add_executable (test_misc
               test_misc.cpp
               test_misc_accumulators.h
               test_misc_interval_tree.h
               test_misc_bit_twiddling.h
               test_misc_long_word.h
               test_misc_edit_environment.h)
target_link_libraries (test_misc ${SEQAN_LIBRARIES})

add_executable (test_misc_cmdparser 
               test_misc_cmdparser.cpp)
target_link_libraries (test_misc_cmdparser ${SEQAN_LIBRARIES})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Register with CTest
# ----------------------------------------------------------------------------

add_test (NAME test_test_misc COMMAND $<TARGET_FILE:test_misc>)
add_test (NAME test_test_misc_cmdparser COMMAND $<TARGET_FILE:test_misc_cmdparser>)
