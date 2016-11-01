# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /tests/align/CMakeLists.txt
#
# CMakeLists.txt file for the align module tests.
# ===========================================================================

cmake_minimum_required (VERSION 3.0.0)
project (seqan_tests_align CXX)
message (STATUS "Configuring tests/align")

set (ALIGN_SIMD_TEST TRUE CACHE INTERNAL "Whether to build test_align_simd.")
# workaround a bug in llvm35 on FreeBSD
if ((${CMAKE_SYSTEM_NAME} STREQUAL "FreeBSD") AND
    (COMPILER_CLANG) AND
    (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 3.6.0))
    set (ALIGN_SIMD_TEST FALSE CACHE INTERNAL "Whether to build test_align_simd.")
endif ()

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
add_executable (test_align
                test_align.cpp
                test_align_stream.h
                test_align_fragment.h
                test_align_gaps.h
                test_align_gaps_iterator.h
                test_alignment_dp_band.h
                test_alignment_dp_formula.h
                test_alignment_dp_cell.h
                test_alignment_dp_profile.h
                test_alignment_dp_matrix.h
                test_alignment_dp_matrix_navigator.h
                test_alignment_dp_trace_segment.h
                test_alignment_dp_adapt_tracesegments.h
                test_alignment_dp_traceback.h
                test_alignment_algorithms_band_position.h
                test_alignment_algorithms_global.h
                test_alignment_algorithms_local.h
                test_alignment_algorithms_global_banded.h
                test_alignment_algorithms_local_banded.h
                test_align_global_alignment_specialized.h
                test_evaluate_alignment.h)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (test_align ${SEQAN_LIBRARIES})

if (ALIGN_SIMD_TEST)
    # Add executable for simd tests.
    add_executable (test_align_simd
                    test_align_simd.cpp
                    test_align_simd.h)

    # Add dependencies found by find_package (SeqAn).
    target_link_libraries (test_align_simd ${SEQAN_LIBRARIES})
endif()

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Register with CTest
# ----------------------------------------------------------------------------

add_test (NAME test_test_align COMMAND $<TARGET_FILE:test_align>)
if (ALIGN_SIMD_TEST)
    add_test (NAME test_test_align_simd COMMAND $<TARGET_FILE:test_align_simd>)
endif ()
