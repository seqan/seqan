# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /core/tests/sequence_journaled/CMakeLists.txt
#
# CMakeLists.txt file for the sequence_journaled module tests.
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (seqan_core_tests_sequence_journaled)
message (STATUS "Configuring core/tests/sequence_journaled")

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
add_executable (test_sequence_journaled
               test_sequence_journaled.cpp
               test_sequence_journaled.h
               test_sequence_journaled_iterator.h)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (test_sequence_journaled ${SEQAN_LIBRARIES})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS})

# ----------------------------------------------------------------------------
# Register with CTest
# ----------------------------------------------------------------------------

add_test (NAME test_test_sequence_journaled COMMAND $<TARGET_FILE:test_sequence_journaled>)
