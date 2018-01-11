# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /tests/index/CMakeLists.txt
#
# CMakeLists.txt file for the index module tests.
# ===========================================================================

cmake_minimum_required (VERSION 3.0.0)
project (seqan_tests_index CXX)
message (STATUS "Configuring tests/index")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search SeqAn and select dependencies.
if (NOT "${SEQAN_BUILD_SYSTEM}" STREQUAL "DEVELOP")
    find_package (OpenMP)
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
add_executable (test_index_creation
                test_index_creation.cpp
                test_index_creation.h
                test_index_helpers.h)
target_link_libraries (test_index_creation ${SEQAN_LIBRARIES})

add_executable (test_index_crosscompare_char
                test_index_crosscompare_char.cpp
                test_crosscompare.h)
target_link_libraries (test_index_crosscompare_char ${SEQAN_LIBRARIES})

# We need an additional TU for the DFI test otherwise GCC 4.7 warns in 32 bit mode.
add_executable (test_index_crosscompare_char_dfi
                test_index_crosscompare_char_dfi.cpp
                test_crosscompare.h)
target_link_libraries (test_index_crosscompare_char_dfi ${SEQAN_LIBRARIES})

add_executable (test_index_crosscompare_dna
                test_index_crosscompare_dna.cpp
                test_crosscompare.h)
target_link_libraries (test_index_crosscompare_dna ${SEQAN_LIBRARIES})

# We need an additional TU for the DFI test otherwise GCC 4.7 warns in 32 bit mode.
add_executable (test_index_crosscompare_dna_dfi
                test_index_crosscompare_dna_dfi.cpp
                test_crosscompare.h)
target_link_libraries (test_index_crosscompare_dna_dfi ${SEQAN_LIBRARIES})

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

#add_executable (test_index_fm_right_array_binary_tree
#                test_index_fm_right_array_binary_tree.cpp
#                test_index_fm_right_array_binary_tree.h
#                test_index_fm_right_array_binary_tree_iterator.h
#                test_index_helpers.h)
#target_link_libraries (test_index_fm_right_array_binary_tree ${SEQAN_LIBRARIES})

add_executable (test_index_fm_rank_dictionary
                test_index_fm_rank_dictionary.cpp
                test_index_helpers.h)
target_link_libraries (test_index_fm_rank_dictionary ${SEQAN_LIBRARIES})

add_executable (test_index_fm_sparse_string
                test_index_fm_sparse_string.cpp)
target_link_libraries (test_index_fm_sparse_string ${SEQAN_LIBRARIES})

add_executable (test_index_base
                test_index_base.cpp
                test_index_helpers.h)
target_link_libraries (test_index_base ${SEQAN_LIBRARIES})

add_executable (test_index_fm
                test_index_fm.cpp
                test_index_helpers.h)
target_link_libraries (test_index_fm ${SEQAN_LIBRARIES})

add_executable (test_index_bifm
                test_index_bifm.cpp
                test_index_helpers.h)
target_link_libraries (test_index_bifm ${SEQAN_LIBRARIES})

add_executable (test_index_vstree
                test_index_vstree.cpp
                test_index_fm_stree.h
                test_index_helpers.h)
target_link_libraries (test_index_vstree ${SEQAN_LIBRARIES})

# GCC < 4.5 segfaults when compiling the test_index_stree_iterators test.
# Enable if the compiler is not GCC or the version is >= 4.5.
if (NOT CMAKE_COMPILER_IS_GNUCXX OR (450 LESS _GCC_VERSION))
    add_executable (
        test_index_stree_iterators
        test_index_stree_iterators.cpp
        test_stree_iterators.h)
    target_link_libraries (test_index_stree_iterators ${SEQAN_LIBRARIES})
endif ()

add_executable (test_index_view
                test_index_view.cpp
                test_index_view.h)
target_link_libraries (test_index_view ${SEQAN_LIBRARIES})

add_executable (test_index_finder
                test_index_finder.cpp
                test_index_helpers.h)
target_link_libraries (test_index_finder ${SEQAN_LIBRARIES})

add_executable (test_find_backtracking
                test_find_backtracking.cpp
                test_find_backtracking.h
                test_find_base.h)
target_link_libraries (test_find_backtracking ${SEQAN_LIBRARIES})

add_executable (test_find2_index_approx
                test_find2_index_approx.cpp
                test_find2_index_approx.h)
target_link_libraries (test_find2_index_approx ${SEQAN_LIBRARIES})

add_executable (test_index_repeats
               test_index_repeats.cpp
               test_index_repeats.h)
target_link_libraries (test_index_repeats ${SEQAN_LIBRARIES})

add_executable (test_index_swift
                test_index_swift.cpp)
target_link_libraries (test_index_swift ${SEQAN_LIBRARIES})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Register with CTest
# ----------------------------------------------------------------------------

add_test (NAME test_test_index_creation COMMAND $<TARGET_FILE:test_index_creation>)
add_test (NAME test_test_index_crosscompare_char COMMAND $<TARGET_FILE:test_index_crosscompare_char>)
add_test (NAME test_test_index_crosscompare_char_dfi COMMAND $<TARGET_FILE:test_index_crosscompare_char_dfi>)
add_test (NAME test_test_index_crosscompare_dna COMMAND $<TARGET_FILE:test_index_crosscompare_dna>)
add_test (NAME test_test_index_crosscompare_dna_dfi COMMAND $<TARGET_FILE:test_index_crosscompare_dna_dfi>)
add_test (NAME test_test_index_qgram COMMAND $<TARGET_FILE:test_index_qgram>)
add_test (NAME test_test_index_sa_bwtwalk COMMAND $<TARGET_FILE:test_index_sa_bwtwalk>)
add_test (NAME test_test_index_shapes COMMAND $<TARGET_FILE:test_index_shapes>)
add_test (NAME test_test_index_drawing COMMAND $<TARGET_FILE:test_index_drawing>)
#add_test (NAME test_test_index_fm_right_array_binary_tree COMMAND $<TARGET_FILE:test_index_fm_right_array_binary_tree>)
add_test (NAME test_test_index_fm_rank_dictionary COMMAND $<TARGET_FILE:test_index_fm_rank_dictionary>)
add_test (NAME test_test_index_fm_sparse_string COMMAND $<TARGET_FILE:test_index_fm_sparse_string>)
add_test (NAME test_test_index_base COMMAND $<TARGET_FILE:test_index_base>)
add_test (NAME test_test_index_fm COMMAND $<TARGET_FILE:test_index_fm>)
add_test (NAME test_test_index_bifm COMMAND $<TARGET_FILE:test_index_bifm>)
add_test (NAME test_test_index_vstree COMMAND $<TARGET_FILE:test_index_vstree>)
if (NOT CMAKE_COMPILER_IS_GNUCXX OR (450 LESS _GCC_VERSION))
    add_test (NAME test_test_index_stree_iterators COMMAND $<TARGET_FILE:test_index_stree_iterators>)
endif ()
add_test (NAME test_test_index_view COMMAND $<TARGET_FILE:test_index_view>)
add_test (NAME test_test_index_finder COMMAND $<TARGET_FILE:test_index_finder>)
add_test (NAME test_test_find_backtracking COMMAND $<TARGET_FILE:test_find_backtracking>)
add_test (NAME test_test_index_repeats COMMAND $<TARGET_FILE:test_index_repeats>)
