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

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search SeqAn and select dependencies.
find_package (SeqAn REQUIRED)

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

# Add executable for simd tests.
add_executable(test_align_simd
               test_align_simd.cpp
               test_align_simd.h)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (test_align_simd ${SEQAN_LIBRARIES})

#TODO(rrahn) Enable for windows platforms.
if (MSVC)
  # Disable vectorization support on windows for now.
  # SET (_FLAGS "/arch:AVX")
else (MSVC)
  SET (_FLAGS "-march=native")
# Check why we would check for the abi version?
#if (CMAKE_COMPILER_IS_GNUCXX)
#    SET (_FLAGS "${_FLAGS} -fabi-version=6")
#    if (APPLE)
      #exec_program("/usr/sbin/sysctl -n hw.optional.avx2_0" OUTPUT_VARIABLE _HAS_AVX2)
      #if (NOT "${_HAS_AVX2}" EQUAL "1")
#      SET (_FLAGS "${_FLAGS} -mno-avx")
      #endif ()
#    endif (APPLE)
#  endif (CMAKE_COMPILER_IS_GNUCXX)
  SET_TARGET_PROPERTIES(test_align_simd PROPERTIES COMPILE_FLAGS ${_FLAGS})
endif (MSVC)

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Register with CTest
# ----------------------------------------------------------------------------

add_test (NAME test_test_align COMMAND $<TARGET_FILE:test_align>)
add_test (NAME test_test_align_simd COMMAND $<TARGET_FILE:test_align_simd>)
