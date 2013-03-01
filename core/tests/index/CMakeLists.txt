# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /core/tests/index/CMakeLists.txt
#
# CMakeLists.txt file for the index module tests.
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (seqan_core_tests_index)
message (STATUS "Configuring core/tests/index")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search SeqAn and select dependencies.
set (SEQAN_FIND_DEPENDENCIES OpenMP)
find_package (SeqAn REQUIRED)

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Update the list of file names below if you add source files to your test.
add_executable (test_index_creation
                test_index_creation.cpp)
target_link_libraries (test_index_creation ${SEQAN_LIBRARIES})

add_executable (test_index_crosscompare_char
                test_index_crosscompare_char.cpp
                test_crosscompare.h)
target_link_libraries (test_index_crosscompare_char ${SEQAN_LIBRARIES})

add_executable (test_index_crosscompare_dna
                test_index_crosscompare_dna.cpp
                test_crosscompare.h)
target_link_libraries (test_index_crosscompare_dna ${SEQAN_LIBRARIES})

add_executable (test_index_qgram
                test_index_qgram.cpp
                test_qgram_index.h)
target_link_libraries (test_index_qgram ${SEQAN_LIBRARIES})

add_executable (test_index_sa_bwtwalk
                test_index_sa_bwtwalk.cpp
                test_sa_bwtwalk.h)
target_link_libraries (test_index_sa_bwtwalk ${SEQAN_LIBRARIES})

add_executable (test_index_shapes
                test_index_shapes.cpp
                test_shapes.h)
target_link_libraries (test_index_shapes ${SEQAN_LIBRARIES})

add_executable (test_index_drawing
                test_index_drawing.cpp)
target_link_libraries (test_index_drawing ${SEQAN_LIBRARIES})

add_executable (test_index_fm
                test_index_fm.cpp
			    test_index_fm.h
			    test_index_fm_compressed_sa_iterator.h
			    test_index_fm_compressed_sa.h
			    test_index_fm_prefix_sum_table.h
			    test_index_fm_rank_support_bit_string_iterator.h
			    test_index_fm_rank_support_bit_string.h
			    test_index_fm_right_array_binary_tree_iterator.h
			    test_index_fm_right_array_binary_tree.h
			    test_index_fm_sparse_string.h
			    test_index_fm_stree.h)
target_link_libraries (test_index_fm ${SEQAN_LIBRARIES})

add_executable (test_index_fm_rank_dictionary
                test_index_fm_rank_dictionary.cpp
			    test_index_fm_rank_dictionary.h
			    test_index_fm_sentinel_rank_dictionary.h)
target_link_libraries (test_index_fm_rank_dictionary ${SEQAN_LIBRARIES})


# GCC < 4.5 segfaults when compiling the test_index_stree_iterators test.
# Enable if the compiler is not GCC or the version is >= 4.5.
if (NOT CMAKE_COMPILER_IS_GNUCXX OR (450 LESS _GCC_VERSION))
    add_executable (
        test_index_stree_iterators
        test_index_stree_iterators.cpp
        test_stree_iterators.h)
    target_link_libraries (test_index_stree_iterators ${SEQAN_LIBRARIES})
endif ()

add_executable (test_index_repeats
               test_index_repeats.cpp
               test_index_repeats.h)
target_link_libraries (test_index_repeats ${SEQAN_LIBRARIES})


# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Register with CTest
# ----------------------------------------------------------------------------

add_test (NAME test_test_index_creation COMMAND $<TARGET_FILE:test_index_creation>)
add_test (NAME test_test_index_crosscompare_char COMMAND $<TARGET_FILE:test_index_crosscompare_char>)
add_test (NAME test_test_index_crosscompare_dna COMMAND $<TARGET_FILE:test_index_crosscompare_dna>)
add_test (NAME test_test_index_qgram COMMAND $<TARGET_FILE:test_index_qgram>)
add_test (NAME test_test_index_sa_bwtwalk COMMAND $<TARGET_FILE:test_index_sa_bwtwalk>)
add_test (NAME test_test_index_shapes COMMAND $<TARGET_FILE:test_index_shapes>)
add_test (NAME test_test_index_drawing COMMAND $<TARGET_FILE:test_index_drawing>)
add_test (NAME test_test_index_fm COMMAND $<TARGET_FILE:test_index_fm>)
if (NOT CMAKE_COMPILER_IS_GNUCXX OR (450 LESS _GCC_VERSION))
    add_test (NAME test_test_index_stree_iterators COMMAND $<TARGET_FILE:test_index_stree_iterators>)
endif ()
add_test (NAME test_test_index_repeats COMMAND $<TARGET_FILE:test_index_repeats>)

