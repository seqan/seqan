# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /%(CMAKE_PROJECT_PATH)s/CMakeLists.txt
#
# CMakeLists.txt file for the %(NAME)s module tests.
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (seqan_%(CMAKE_PROJECT_NAME)s)
message (STATUS "Configuring %(CMAKE_PROJECT_PATH)s")

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
add_executable (test_%(NAME)s
                test_%(NAME)s.cpp
                test_%(NAME)s.h)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (test_%(NAME)s ${SEQAN_LIBRARIES})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS})

# ----------------------------------------------------------------------------
# Register with CTest
# ----------------------------------------------------------------------------

add_test (NAME test_test_%(NAME)s COMMAND $<TARGET_FILE:test_%(NAME)s>)
