# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /tests/seeds/CMakeLists.txt
#
# CMakeLists.txt file for the seeds module tests.
# ===========================================================================

cmake_minimum_required (VERSION 3.0.0)
project (seqan_tests_seeds CXX)
message (STATUS "Configuring tests/seeds")

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

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# Create executables for all seeds module tests.
add_executable (test_seeds_combination
                test_seeds_combination.cpp)
add_executable (test_seeds_extension
                test_seeds_extension.cpp)
add_executable (test_seeds_global_chaining
                test_seeds_global_chaining.cpp)
add_executable (test_seeds_seed_base
                test_seeds_seed_base.cpp)
add_executable (test_seeds_seed_chained
                test_seeds_seed_chained.cpp)
add_executable (test_seeds_seed_diagonal
                test_seeds_seed_diagonal.cpp)
add_executable (test_seeds_seed_set_base
                test_seeds_seed_set_base.cpp)
add_executable (test_seeds_seed_set_unordered
                test_seeds_seed_set_unordered.cpp)
add_executable (test_seeds_align_banded_chain_impl
                test_align_banded_chain_impl.cpp)
add_executable (test_seeds_banded_chain_alignment_interface
                test_banded_chain_alignment_interface.cpp)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (test_seeds_combination ${SEQAN_LIBRARIES})
target_link_libraries (test_seeds_extension ${SEQAN_LIBRARIES})
target_link_libraries (test_seeds_global_chaining ${SEQAN_LIBRARIES})
target_link_libraries (test_seeds_seed_base ${SEQAN_LIBRARIES})
target_link_libraries (test_seeds_seed_chained ${SEQAN_LIBRARIES})
target_link_libraries (test_seeds_seed_diagonal ${SEQAN_LIBRARIES})
target_link_libraries (test_seeds_seed_set_base ${SEQAN_LIBRARIES})
target_link_libraries (test_seeds_seed_set_unordered ${SEQAN_LIBRARIES})
target_link_libraries (test_seeds_align_banded_chain_impl ${SEQAN_LIBRARIES})
target_link_libraries (test_seeds_banded_chain_alignment_interface ${SEQAN_LIBRARIES})

# ----------------------------------------------------------------------------
# Register with CTest
# ----------------------------------------------------------------------------

add_test (NAME test_test_seeds_combination COMMAND $<TARGET_FILE:test_seeds_combination>)
add_test (NAME test_test_seeds_extension COMMAND $<TARGET_FILE:test_seeds_extension>)
add_test (NAME test_test_seeds_global_chaining COMMAND $<TARGET_FILE:test_seeds_global_chaining>)
add_test (NAME test_test_seeds_seed_base COMMAND $<TARGET_FILE:test_seeds_seed_base>)
add_test (NAME test_test_seeds_seed_chained COMMAND $<TARGET_FILE:test_seeds_seed_chained>)
add_test (NAME test_test_seeds_seed_diagonal COMMAND $<TARGET_FILE:test_seeds_seed_diagonal>)
add_test (NAME test_test_seeds_seed_set_base COMMAND $<TARGET_FILE:test_seeds_seed_set_base>)
add_test (NAME test_test_seeds_seed_set_unordered COMMAND $<TARGET_FILE:test_seeds_seed_set_unordered>)
add_test (NAME test_test_seeds_align_banded_chain_impl COMMAND $<TARGET_FILE:test_seeds_align_banded_chain_impl>)
add_test (NAME test_test_seeds_banded_chain_alignment_interface COMMAND $<TARGET_FILE:test_seeds_banded_chain_alignment_interface>)
