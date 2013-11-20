# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /core/apps/razers2/CMakeLists.txt
#
# CMakeLists.txt file for RazerS 2.
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (core_apps_razers2)
message (STATUS "Configuring core/apps/razers2")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search SeqAn and select dependencies.
# TODO(holtgrew): Specifying ZLIB here is for demo purposes only. Remove again.
set (SEQAN_FIND_DEPENDENCIES ZLIB)
find_package (SeqAn REQUIRED)

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# TODO(holtgrew): Add dependency on SeqAn library headers.

# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Add definitions set by the build system.
add_definitions (-DSEQAN_REVISION="${SEQAN_REVISION}")
add_definitions (-DSEQAN_DATE="${SEQAN_DATE}")

# Update the list of file names below if you add source files to your application.
add_executable (razers2 razers.cpp
                       razers.h
                       param_tabs.h
                       param_tabs.cpp
                       param_tabs.inc
                       outputFormat.h
                       paramChooser.h
                       razers_matepairs.h
                       readSimulator.h)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (razers2 ${SEQAN_LIBRARIES})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------

# Set variables for installing, depending on the selected build type.
if (NOT SEQAN_PREFIX_SHARE_DOC)
  seqan_setup_install_vars (razers2)
endif (NOT SEQAN_PREFIX_SHARE_DOC)

# Install razers2 in ${PREFIX}/bin directory
install (TARGETS razers2
         DESTINATION bin)

# Install non-binary files for the package to "." for app builds and
# ${PREFIX}/share/doc/razers2 for SeqAn release builds.
install (FILES LICENSE
               README
         DESTINATION ${SEQAN_PREFIX_SHARE_DOC})
install (FILES example/genome.fa
               example/reads.fa
               example/reads2.fa
         DESTINATION ${SEQAN_PREFIX_SHARE_DOC}/example)

# ----------------------------------------------------------------------------
# App Test
# ----------------------------------------------------------------------------

seqan_add_app_test (razers2)

# ----------------------------------------------------------------------------
# Setup Common Tool Description for Generic Workflow Nodes
# ----------------------------------------------------------------------------

# Include executable razers in CTD structure.
set (SEQAN_CTD_EXECUTABLES ${SEQAN_CTD_EXECUTABLES} razers2 CACHE INTERNAL "")

# ----------------------------------------------------------------------------
# CPack Install
# ----------------------------------------------------------------------------

if (SEQAN_BUILD_SYSTEM STREQUAL "APP:razers2")
  set (CPACK_PACKAGE_NAME "razers2")
  set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "RazerS 2 - Read Mapping With Controllable Sensitivity")
  set (CPACK_DEBIAN_PACKAGE_MAINTAINER "David Weese <david.weese@fu-berlin.de>")
  set (CPACK_PACKAGE_VENDOR "SeqAn Team, FU Berlin")

  seqan_configure_cpack_app (razers2 "RazerS 2")
endif (SEQAN_BUILD_SYSTEM STREQUAL "APP:razers2")
