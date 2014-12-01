# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /extras/tests/align_profile/CMakeLists.txt
#
# CMakeLists.txt file for the align_profile module tests.
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (seqan_extras_tests_align_profile)
message (STATUS "Configuring extras/tests/align_profile")

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

# The test program.
add_executable (test_align_profile test_align_profile.cpp)
target_link_libraries (test_align_profile ${SEQAN_LIBRARIES})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Register with CTest
# ----------------------------------------------------------------------------

add_test (NAME test_test_align_profile COMMAND $<TARGET_FILE:test_align_profile>)
