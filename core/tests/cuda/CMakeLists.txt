# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /core/tests/cuda/CMakeLists.txt
#
# CMakeLists.txt file for the cuda module tests.
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (seqan_core_tests_cuda)
message (STATUS "Configuring core/tests/cuda")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search SeqAn and select dependencies.
set (SEQAN_FIND_DEPENDENCIES CUDA)
find_package (SeqAn REQUIRED)

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Setup CUDA flags.
seqan_setup_cuda_vars (ARCH sm_20 DEBUG_DEVICE DISABLE_WARNINGS)

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Update the list of file names below if you add source files to your test.
if (SEQAN_HAS_CUDA)
  cuda_add_executable (test_cuda_basic test_cuda_basic.cu)
  cuda_add_executable (test_cuda_sequence test_cuda_sequence.cu test_cuda_common.h)
  cuda_add_executable (test_cuda_index test_cuda_index.cu test_cuda_common.h)
endif ()

# Add dependencies found by find_package (SeqAn).
if (SEQAN_HAS_CUDA)
  target_link_libraries (test_cuda_basic ${SEQAN_LIBRARIES})
  target_link_libraries (test_cuda_sequence ${SEQAN_LIBRARIES})
  target_link_libraries (test_cuda_index ${SEQAN_LIBRARIES})
endif ()

# ----------------------------------------------------------------------------
# Register with CTest
# ----------------------------------------------------------------------------

if (SEQAN_HAS_CUDA)
    add_test (NAME test_test_cuda_basic COMMAND $<TARGET_FILE:test_cuda_basic>)
    add_test (NAME test_test_cuda_sequence COMMAND $<TARGET_FILE:test_cuda_sequence>)
    add_test (NAME test_test_cuda_index COMMAND $<TARGET_FILE:test_cuda_index>)
endif ()
