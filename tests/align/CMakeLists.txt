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
if ((SEQAN_TRAVIS_BUILD) OR
    ((${CMAKE_SYSTEM_NAME} STREQUAL "FreeBSD") AND
    (COMPILER_CLANG) AND
    (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 3.6.0)))
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
    add_executable (test_align_simd_global_equal_length
                    test_align_simd_global_equal_length.cpp
                    test_align_simd_base.h
                    test_align_simd_global.h
                    test_mock.h)

    add_executable (test_align_simd_global_variable_length
                    test_align_simd_global_variable_length.cpp
                    test_align_simd_base.h
                    test_align_simd_global.h
                    test_mock.h)

    add_executable (test_align_simd_local_equal_length
                    test_align_simd_local_equal_length.cpp
                    test_align_simd_base.h
                    test_align_simd_local.h
                    test_mock.h)

    add_executable (test_align_simd_local_variable_length
                    test_align_simd_local_variable_length.cpp
                    test_align_simd_base.h
                    test_align_simd_local.h
                    test_mock.h)

    # Add dependencies found by find_package (SeqAn).
    target_link_libraries (test_align_simd_global_equal_length ${SEQAN_LIBRARIES})
    target_link_libraries (test_align_simd_global_variable_length ${SEQAN_LIBRARIES})
    target_link_libraries (test_align_simd_local_equal_length ${SEQAN_LIBRARIES})
    target_link_libraries (test_align_simd_local_variable_length ${SEQAN_LIBRARIES})
    # note(marehr): there is a bug when using <=clang3.8 with gcc4.9's stdlib,
    # where the default -ftemplate-depth=256 of clang is insufficient.
    # test_align_simd_avx2 needs a depth of at least 266.
    if (COMPILER_CLANG)
      target_compile_options(test_align_simd_global_equal_length PRIVATE -ftemplate-depth=1024)
      target_compile_options(test_align_simd_global_variable_length PRIVATE -ftemplate-depth=1024)
      target_compile_options(test_align_simd_local_equal_length PRIVATE -ftemplate-depth=1024)
      target_compile_options(test_align_simd_local_variable_length PRIVATE -ftemplate-depth=1024)
    endif()
endif()

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Register with CTest
# ----------------------------------------------------------------------------

add_test (NAME test_test_align COMMAND $<TARGET_FILE:test_align>)
if (ALIGN_SIMD_TEST)
    include (SeqAnSimdUtility)
    add_simd_platform_tests(test_align_simd_global_equal_length)
    add_simd_platform_tests(test_align_simd_global_variable_length)
    add_simd_platform_tests(test_align_simd_local_equal_length)
    add_simd_platform_tests(test_align_simd_local_variable_length)
endif ()
