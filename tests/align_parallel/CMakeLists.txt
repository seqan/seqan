# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /tests/align_parallel/CMakeLists.txt
#
# CMakeLists.txt file for the align_parallel module tests.
# ===========================================================================

cmake_minimum_required (VERSION 3.0.0)
project (seqan_tests_align_parallel CXX)
message (STATUS "Configuring tests/align_parallel")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

set (ALIGN_PARALLEL_SIMD_TEST TRUE CACHE INTERNAL "Whether to build test_align_parallel.")
# workaround a bug in llvm35 on FreeBSD
if ((${CMAKE_SYSTEM_NAME} STREQUAL "FreeBSD") AND
    (COMPILER_CLANG) AND
    (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 3.6.0))
    set (ALIGN_PARALLEL_SIMD_TEST FALSE CACHE INTERNAL "Whether to build test_align_parallel.")
endif ()

# Increase recursive template instantiation depth for clang compiler.
if (COMPILER_CLANG)
    set(SEQAN_CXX_FLAGS "${SEQAN_CXX_FLAGS} -ftemplate-depth=1024")
endif ()

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Update the list of file names below if you add source files to your test.
add_executable (test_align_parallel_data_structures
                test_align_parallel_data_structures.cpp
                test_align_wavefront_task_scheduler.h
                test_align_wavefront_alignment_scheduler.h
                test_align_wavefront_intermediate_dp_result.h
                test_align_wavefront_alignment_thread_local.h)

add_executable (test_align_parallel_algorithm
                test_align_parallel_algorithm.cpp
                test_align_parallel_wavefront_alignment.h)

add_executable (test_align_parallel_interface
                test_align_parallel_interface.cpp
                test_align_parallel_interface.h
                ../align/test_mock.h)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (test_align_parallel_data_structures ${SEQAN_LIBRARIES})

# Add dependencies found by find_package (SeqAn).
target_link_libraries (test_align_parallel_algorithm ${SEQAN_LIBRARIES})

# Add dependencies found by find_package (SeqAn).
target_link_libraries (test_align_parallel_interface ${SEQAN_LIBRARIES})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Register with CTest
# ----------------------------------------------------------------------------

add_test (NAME test_align_parallel_data_structures COMMAND $<TARGET_FILE:test_align_parallel_data_structures>)
if (ALIGN_PARALLEL_SIMD_TEST)
    include (SeqAnSimdUtility)
    add_simd_platform_tests(test_align_parallel_interface)
    add_simd_platform_tests(test_align_parallel_algorithm)
endif ()
