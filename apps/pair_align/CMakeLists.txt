# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /core/apps/pair_align/CMakeLists.txt
#
# CMakeLists.txt file for pair_align.
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (core_apps_pair_align)
message (STATUS "Configuring core/apps/pair_align")

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
add_executable (pair_align pair_align.cpp)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (pair_align ${SEQAN_LIBRARIES})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------

# Set variables for installing, depending on the selected build type.
if (NOT SEQAN_PREFIX_SHARE_DOC)
  seqan_setup_install_vars (pair_align)
endif (NOT SEQAN_PREFIX_SHARE_DOC)

# Install pair_align in ${PREFIX}/bin directory
install (TARGETS pair_align
         DESTINATION bin)

# Install non-binary files for the package to "." for app builds and
# ${PREFIX}/share/doc/pair_align for SeqAn release builds.
install (FILES LICENSE
               README
         DESTINATION ${SEQAN_PREFIX_SHARE_DOC})

# ----------------------------------------------------------------------------
# App Test
# ----------------------------------------------------------------------------

seqan_add_app_test (pair_align)

# ----------------------------------------------------------------------------
# Setup Command Tool Description for Generic Workflow Nodes
# ----------------------------------------------------------------------------

# Currently, pair_align uses a lot of memory to compile so we do not add it
# to the list of programs to build for workflow engines.

#Include executable pair_align in CTD structure.
#set (SEQAN_CTD_EXECUTABLES ${SEQAN_CTD_EXECUTABLES} pair_align CACHE INTERNAL "")

# ----------------------------------------------------------------------------
# CPack Install
# ----------------------------------------------------------------------------

if (SEQAN_BUILD_SYSTEM STREQUAL "APP:pair_align")
  set (CPACK_PACKAGE_NAME "pair_align")
  set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "pair_align - Pairwise Sequence Alignment")
  set (CPACK_DEBIAN_PACKAGE_MAINTAINER "Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>")
  set (CPACK_PACKAGE_VENDOR "SeqAn Team, FU Berlin")

  seqan_configure_cpack_app (pair_align "pair_align")
endif (SEQAN_BUILD_SYSTEM STREQUAL "APP:pair_align")

