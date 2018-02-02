# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /tests/basic/CMakeLists.txt
#
# CMakeLists.txt file for the basic module tests.
# ===========================================================================

cmake_minimum_required (VERSION 3.0.0)
project (seqan_tests_simd CXX)
message (STATUS "Configuring tests/simd")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search SeqAn and select dependencies.
if (NOT "${SEQAN_BUILD_SYSTEM}" STREQUAL "DEVELOP")
    find_package (SeqAn REQUIRED)
endif ()

include (SeqAnSimdUtility)

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Update the list of file names below if you add source files to your test.
add_executable (
  test_simd_vector
  test_simd_vector.cpp
  test_simd_vector.h)
target_link_libraries (test_simd_vector ${SEQAN_LIBRARIES})
if (COMPILER_CLANG)
  target_compile_options(test_simd_vector PRIVATE -ftemplate-depth=1024)
endif()

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Register with CTest
# ----------------------------------------------------------------------------

add_simd_platform_tests(test_simd_vector)
