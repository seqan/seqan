# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /extras/apps/fx_tools/CMakeLists.txt
#
# CMakeLists.txt file for fx_tools.
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (seqan_extras_apps_fx_tools)
message (STATUS "Configuring extras/apps/fx_tools")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search SeqAn and select dependencies.
#set (SEQAN_FIND_DEPENDENCIES NONE)
find_package (SeqAn REQUIRED)

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Setup SeqAn requirements.
include_directories (${SEQAN_INCLUDE_DIRS})
add_definitions (${SEQAN_DEFINITIONS})
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# Computation of coverage and C+G content from genome and mapping.
add_executable (fx_bam_coverage fx_bam_coverage.cpp)
target_link_libraries (fx_bam_coverage ${SEQAN_LIBRARIES})

# Computation of statistics on FASTQ files.
add_executable (fx_fastq_stats fx_fastq_stats.cpp)
target_link_libraries (fx_fastq_stats ${SEQAN_LIBRARIES})

# ----------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------

# Set variables for installing, depending on the selected build type.
if (NOT SEQAN_PREFIX_SHARE_DOC)
  seqan_setup_install_vars (fx_tools)
endif (NOT SEQAN_PREFIX_SHARE_DOC)

# Install fx_tools in ${PREFIX}/bin directory
install (TARGETS fx_bam_coverage
         DESTINATION bin)

# Install non-binary files for the package to "." for app builds and
# ${PREFIX}/share/doc/fx_tools for SeqAn release builds.
install (FILES LICENSE
               README
         DESTINATION ${SEQAN_PREFIX_SHARE_DOC})
#install (FILES example/fasta1.fa
#               example/fasta2.fa
#         DESTINATION ${SEQAN_PREFIX_SHARE_DOC}/example)

# ----------------------------------------------------------------------------
# App Test
# ----------------------------------------------------------------------------

#seqan_add_app_test (fx_tools)

# ----------------------------------------------------------------------------
# Setup Common Tool Description for Generic Workflow Nodes
# ----------------------------------------------------------------------------

# Include executable razers in CTD structure.
set (SEQAN_CTD_EXECUTABLES ${SEQAN_CTD_EXECUTABLES} fx_fastq_stats fx_bam_coverage CACHE INTERNAL "")

# ----------------------------------------------------------------------------
# CPack Install
# ----------------------------------------------------------------------------

if (SEQAN_BUILD_SYSTEM STREQUAL "APP:fx_tools")
  set (CPACK_PACKAGE_NAME "fx_tools")
  set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "fx_tools")
  set (CPACK_DEBIAN_PACKAGE_MAINTAINER "Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>")
  set (CPACK_PACKAGE_VENDOR "Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>")

  seqan_configure_cpack_app (fx_tools "fx_tools")
endif (SEQAN_BUILD_SYSTEM STREQUAL "APP:fx_tools")

