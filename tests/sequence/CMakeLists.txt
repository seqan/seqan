# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /tests/sequence/CMakeLists.txt
#
# CMakeLists.txt file for the sequence module tests.
# ===========================================================================

cmake_minimum_required (VERSION 3.0.0)
project (seqan_tests_sequence CXX)
message (STATUS "Configuring tests/sequence")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search SeqAn and select dependencies.
if (NOT "${SEQAN_BUILD_SYSTEM}" STREQUAL "DEVELOP")
    find_package (SeqAn REQUIRED)
endif ()
# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Update the list of file names below if you add source files to your test.
add_executable (test_sequence
               test_sequence.cpp
               test_segment.h
               test_sequence_std_adaptions.h
               test_string.h
               test_stringset.h
               test_string_packed_extension.h)
target_link_libraries (test_sequence ${SEQAN_LIBRARIES})

add_executable (test_sequence_v2 test_sequence_v2.cpp
               test_sequence.h
               test_segment_beta.h)
target_link_libraries (test_sequence_v2 ${SEQAN_LIBRARIES})

add_executable (test_stringset_v2 test_stringset_v2.cpp
               test_string_set.h)
target_link_libraries (test_stringset_v2 ${SEQAN_LIBRARIES})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Register with CTest
# ----------------------------------------------------------------------------

add_test (NAME test_test_sequence COMMAND $<TARGET_FILE:test_sequence>)
add_test (NAME test_test_sequence_v2 COMMAND $<TARGET_FILE:test_sequence_v2>)
add_test (NAME test_test_stringset_v2 COMMAND $<TARGET_FILE:test_stringset_v2>)
