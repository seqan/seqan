# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /core/tests/sequence/CMakeLists.txt
#
# CMakeLists.txt file for the sequence module tests.
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (seqan_core_tests_sequence)
message (STATUS "Configuring core/tests/sequence")

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
add_executable (test_sequence
               test_sequence.cpp
               test_segment.h
               test_sequence_std_adaptions.h
               test_string.h
               test_stringset.h)
target_link_libraries (test_sequence ${SEQAN_LIBRARIES})

# Only compile sequence tests if explicitely enabled.
option (SEQAN_COMPILE_SEQUENCE_TESTS_V2 "Enable new sequence tests." FALSE)
if (SEQAN_COMPILE_SEQUENCE_TESTS_V2)
add_executable (test_sequence_v2 test_sequence_v2.cpp
               test_sequence.h
               test_segment_beta.h)
target_link_libraries (test_sequence_v2 ${SEQAN_LIBRARIES})
endif ()

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Register with CTest
# ----------------------------------------------------------------------------

add_test (NAME test_test_sequence COMMAND $<TARGET_FILE:test_sequence>)
if (SEQAN_COMPILE_SEQUENCE_TESTS_V2)
add_test (NAME test_test_sequence_v2 COMMAND $<TARGET_FILE:test_sequence_v2>)
endif ()
