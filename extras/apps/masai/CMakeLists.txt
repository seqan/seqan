# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /extras/apps/masai/CMakeLists.txt
#
# CMakeLists.txt file for Masai.
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (extras_apps_masai)
message (STATUS "Configuring extras/apps/masai")

# ----------------------------------------------------------------------------
# App-Level Configuration
# ----------------------------------------------------------------------------

option (SEQAN_EXTRAS_MASAI_DISABLE_MMAP "Set to ON to disable memory mapping of the reference genome index." OFF)

if (SEQAN_EXTRAS_MASAI_DISABLE_MMAP)
    add_definitions (-DSEQAN_EXTRAS_MASAI_DISABLE_MMAP=1)
endif (SEQAN_EXTRAS_MASAI_DISABLE_MMAP)

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search SeqAn and select dependencies.
set (SEQAN_FIND_DEPENDENCIES NONE)
find_package (SeqAn REQUIRED)

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Update the list of file names below if you add source files to your application.
add_executable (masai_indexer indexer.cpp
                             store.h
                             index.h
                             options.h)

add_executable (masai_mapper mapper.cpp
                            mapper.h
                            tags.h
                            store.h
                            index.h
                            seeder.h
                            extender.h
                            manager.h
                            matches.h
                            stream.h
                            writer.h
                            options.h
                            index/index_qgram_stretched.h
                            index/find_backtracking_stretched.h
                            store/store_io.h
                            store/store_io_sam.h)

add_executable (masai_output_se sorter.cpp
                               sorter.h
                               tags.h
                               store.h
                               matches.h
                               stream.h
                               writer.h
                               options.h
                               store/store_io.h
                               store/store_io_sam.h)

add_executable (masai_output_pe pairer.cpp
                               pairer.h
                               tags.h
                               store.h
                               matches.h
                               stream.h
                               writer.h
                               options.h
                               store/store_io.h
                               store/store_io_sam.h)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (masai_indexer ${SEQAN_LIBRARIES})
target_link_libraries (masai_mapper ${SEQAN_LIBRARIES})
target_link_libraries (masai_output_se ${SEQAN_LIBRARIES})
target_link_libraries (masai_output_pe ${SEQAN_LIBRARIES})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------

# Set variables for installing, depending on the selected build type.
if (NOT SEQAN_PREFIX_SHARE_DOC)
  seqan_setup_install_vars (masai)
endif (NOT SEQAN_PREFIX_SHARE_DOC)

# Install masai in ${PREFIX}/bin directory
install (TARGETS masai_indexer masai_mapper masai_output_se masai_output_pe
         DESTINATION bin)

# Install non-binary files for the package to "." for app builds and
# ${PREFIX}/share/doc/masai for SeqAn release builds.
install (FILES LICENSE
               README
         DESTINATION ${SEQAN_PREFIX_SHARE_DOC})
#install (FILES
#         DESTINATION ${SEQAN_PREFIX_SHARE_DOC}/example)

# ----------------------------------------------------------------------------
# App Test
# ----------------------------------------------------------------------------

seqan_add_app_test (masai)

# ----------------------------------------------------------------------------
# Setup Common Tool Description for Generic Workflow Nodes
# ----------------------------------------------------------------------------

# Include executable masai_mapper in CTD structure.
set (SEQAN_CTD_EXECUTABLES ${SEQAN_CTD_EXECUTABLES} masai_mapper CACHE INTERNAL "")

# ----------------------------------------------------------------------------
# CPack Install
# ----------------------------------------------------------------------------

if (SEQAN_BUILD_SYSTEM STREQUAL "APP:masai")
  set (CPACK_PACKAGE_NAME "masai")
  set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "Masai - Fast and accurate read mapping with approximate seeds and multiple backtracking")
  set (CPACK_DEBIAN_PACKAGE_MAINTAINER "Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>")
  set (CPACK_PACKAGE_VENDOR "SeqAn Team, FU Berlin")

  seqan_configure_cpack_app (masai "masai")
endif (SEQAN_BUILD_SYSTEM STREQUAL "APP:masai")
