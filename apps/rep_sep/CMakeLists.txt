# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /apps/rep_sep/CMakeLists.txt
#
# CMakeLists.txt file for rep_sep.
# ===========================================================================

cmake_minimum_required (VERSION 3.0.0)
project (seqan_apps_rep_sep CXX)
message (STATUS "Configuring apps/rep_sep")

set (SEQAN_APP_VERSION "0.1.11")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search SeqAn and select dependencies.
if (NOT "${SEQAN_BUILD_SYSTEM}" STREQUAL "DEVELOP")
    find_package (SeqAn CONFIG REQUIRED)
endif ()

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

# Update the list of file names below if you add source files to your application.
add_executable (rep_sep rep_sep.cpp
                       assembly_parser.h
                       column_scanner.h
                       rgraph.h
                       rgraph_base.h
                       rgraph_construction.h
                       rgraph_heuristics.h
                       rgraph_score.h
                       utils.h)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (rep_sep ${SEQAN_LIBRARIES})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------

# Set variables for installing, depending on the selected build type.
if (NOT SEQAN_PREFIX_SHARE_DOC)
  seqan_setup_install_vars (rep_sep)
endif (NOT SEQAN_PREFIX_SHARE_DOC)

# Install rep_sep in ${PREFIX}/bin directory
install (TARGETS rep_sep
         DESTINATION ${CMAKE_INSTALL_BINDIR})

# Install non-binary files for the package to "." for app builds and
# ${PREFIX}/share/doc/rep_sep for SeqAn release builds.
install (FILES LICENSE
               README
         DESTINATION ${SEQAN_PREFIX_SHARE_DOC})
#install (FILES
#         DESTINATION ${SEQAN_PREFIX_SHARE_DOC}/example)

# ----------------------------------------------------------------------------
# App Test
# ----------------------------------------------------------------------------

#seqan_add_app_test (rep_sep)

# ----------------------------------------------------------------------------
# CPack Install
# ----------------------------------------------------------------------------

if (SEQAN_BUILD_SYSTEM STREQUAL "APP:rep_sep")
  set (CPACK_PACKAGE_NAME "rep_sep")
  set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "rep_sep - Repeat Separation")
  set (CPACK_DEBIAN_PACKAGE_MAINTAINER "Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>")
  set (CPACK_PACKAGE_VENDOR "SeqAn Team, FU Berlin")

  seqan_configure_cpack_app (rep_sep "rep_sep")
endif (SEQAN_BUILD_SYSTEM STREQUAL "APP:rep_sep")

