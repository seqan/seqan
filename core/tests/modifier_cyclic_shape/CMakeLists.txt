# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /core/tests/modifier_cyclic_shape/CMakeLists.txt
#
# CMakeLists.txt file for the modifier_cyclic_shape module tests.
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (seqan_core_tests_modifier_cyclic_shape)
message (STATUS "Configuring core/tests/modifier_cyclic_shape")

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
add_executable (test_modifier_cyclic_shape
                test_modifier_cyclic_shape.cpp
                test_modifier_cyclic_shape.h
		test_modifier_cylcic_shape_iterator.h)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (test_modifier_cyclic_shape ${SEQAN_LIBRARIES})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Register with CTest
# ----------------------------------------------------------------------------

add_test (NAME test_test_modifier_cyclic_shape COMMAND $<TARGET_FILE:test_modifier_cyclic_shape>)
