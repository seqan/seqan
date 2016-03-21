# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /extras/tests/align_extend/CMakeLists.txt
#
# CMakeLists.txt file for the align_extend module tests.
# ===========================================================================

cmake_minimum_required (VERSION 3.0.0)
project (seqan_extras_tests_blast CXX)
message (STATUS "Configuring tests/blast")

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
add_executable (test_blast test_blast.cpp
                           test_blast_statistics.h
                           test_blast_output.h
                           test_blast_input.h
                           test_blast_misc.h)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (test_blast ${SEQAN_LIBRARIES})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Register with CTest
# ----------------------------------------------------------------------------

add_test (NAME test_blast COMMAND $<TARGET_FILE:test_blast>)
