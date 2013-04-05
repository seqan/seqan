# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /extras/tests/multiple_translation_units/CMakeLists.txt
#
# CMakeLists.txt file for the multiple_translation_units module tests.
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (seqan_extras_tests_multiple_translation_units)
message (STATUS "Configuring extras/tests/multiple_translation_units")

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
add_executable (test_multiple_translation_units
                test_multiple_translation_units.cpp
                test_multiple_translation_units_2.cpp)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (test_multiple_translation_units ${SEQAN_LIBRARIES})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Register with CTest
# ----------------------------------------------------------------------------

add_test (NAME test_test_multiple_translation_units COMMAND $<TARGET_FILE:test_multiple_translation_units>)
