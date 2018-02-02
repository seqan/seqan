# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /apps/rabema/CMakeLists.txt
#
# CMakeLists.txt file for Rabema.
# ===========================================================================

cmake_minimum_required (VERSION 3.0.0)
project (seqan_apps_rabema CXX)
message (STATUS "Configuring apps/rabema")

set (SEQAN_APP_VERSION "1.2.10")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search SeqAn and select dependencies.
if (NOT "${SEQAN_BUILD_SYSTEM}" STREQUAL "DEVELOP")
    find_package (ZLIB)
    find_package (SeqAn CONFIG REQUIRED)
endif ()

if (NOT ZLIB_FOUND)
  message (STATUS "  ZLIB not found: RABEMA will not be built.")
  return ()
endif (NOT ZLIB_FOUND)

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Add definitions set by the build system.
add_definitions (-DSEQAN_APP_VERSION="${SEQAN_APP_VERSION}")
add_definitions (-DSEQAN_REVISION="${SEQAN_REVISION}")
add_definitions (-DSEQAN_DATE="${SEQAN_DATE}")

# Preparation of SAM files, i.e. filling SEQ field of each record.
add_executable (rabema_prepare_sam prepare_sam.cpp sorting.h)
target_link_libraries (rabema_prepare_sam ${SEQAN_LIBRARIES})

# Building of GSI file from "perfect" SAM file and reference sequence.
add_executable (rabema_build_gold_standard
               rabema_build_gold_standard.cpp
               curve_smoothing.h
               find_hamming_simple_ext.h
               find_myers_ukkonen_ext.h
               find_myers_ukkonen_reads.h
               io_gsi.h
               rabema_stats.h
               ref_id_mapping.h
               sorting.h)
target_link_libraries (rabema_build_gold_standard ${SEQAN_LIBRARIES})

# Comparison of arbitrary read mapper SAM against GSI file.
add_executable (rabema_evaluate rabema_evaluate.cpp
               curve_smoothing.h
               find_hamming_simple_ext.h
               find_myers_ukkonen_ext.h
               find_myers_ukkonen_reads.h
               io_gsi.h
               rabema_stats.h
               ref_id_mapping.h
               sorting.h)
target_link_libraries (rabema_evaluate ${SEQAN_LIBRARIES})

# Small tool for helping in debugging.
add_executable (rabema_do_search do_search.cpp find_myers_ukkonen_ext.h)
target_link_libraries (rabema_do_search ${SEQAN_LIBRARIES})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------

# Set variables for installing, depending on the selected build type.
if (NOT SEQAN_PREFIX_SHARE_DOC)
  seqan_setup_install_vars (rabema)
endif (NOT SEQAN_PREFIX_SHARE_DOC)

# Install rabema programs in ${PREFIX}/bin directory
install (TARGETS rabema_prepare_sam
                 rabema_build_gold_standard
                 rabema_evaluate
         DESTINATION ${CMAKE_INSTALL_BINDIR})

# Install rabema in /bin directory
install (TARGETS rabema_prepare_sam rabema_build_gold_standard rabema_evaluate
        DESTINATION ${CMAKE_INSTALL_BINDIR})

# Install non-binary files for the package to "." for app builds and
# ${PREFIX}/share/doc/pair_align for SeqAn release builds.
install (FILES LICENSE
               README
         DESTINATION ${SEQAN_PREFIX_SHARE_DOC})

# ----------------------------------------------------------------------------
# App Test
# ----------------------------------------------------------------------------

seqan_add_app_test (rabema)

# ----------------------------------------------------------------------------
# Setup Common Tool Description for Generic Workflow Nodes
# ----------------------------------------------------------------------------

# Include executable razers in CTD structure.
set (SEQAN_CTD_EXECUTABLES ${SEQAN_CTD_EXECUTABLES} rabema_prepare_sam rabema_build_gold_standard rabema_evaluate CACHE INTERNAL "")

# ----------------------------------------------------------------------------
# CPack Install
# ----------------------------------------------------------------------------

if (SEQAN_BUILD_SYSTEM STREQUAL "APP:rabema")
    set (CPACK_PACKAGE_NAME "rabema")
    set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "Rabema - Read Alignment Benchmark")
    set (CPACK_DEBIAN_PACKAGE_MAINTAINER "Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>")
    set (CPACK_PACKAGE_VENDOR "SeqAn Team, FU Berlin")

    seqan_configure_cpack_app (rabema "Rabema")
endif (SEQAN_BUILD_SYSTEM STREQUAL "APP:rabema")

