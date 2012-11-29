# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /extras/apps/indel_simulator/CMakeLists.txt
#
# CMakeLists.txt file for indel_simulator.
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (extras_apps_indel_simulator)
message (STATUS "Configuring extras/apps/indel_simulator")

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
add_executable (indel_simulator
                indel_simulator.cpp
                indel_simulator.h)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (indel_simulator ${SEQAN_LIBRARIES})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------

# Set variables for installing, depending on the selected build type.
if (NOT SEQAN_PREFIX_SHARE_DOC)
  seqan_setup_install_vars (indel_simulator)
endif (NOT SEQAN_PREFIX_SHARE_DOC)

# Install indel_simulator in ${PREFIX}/bin directory
install (TARGETS indel_simulator
         DESTINATION bin)

# Install non-binary files for the package to "." for app builds and
# ${PREFIX}/share/doc/indel_simulator for SeqAn release builds.
install (FILES LICENSE
               README
         DESTINATION ${SEQAN_PREFIX_SHARE_DOC})
#install (FILES
#         DESTINATION ${SEQAN_PREFIX_SHARE_DOC}/example)

# ----------------------------------------------------------------------------
# App Test
# ----------------------------------------------------------------------------

#seqan_add_app_test (indel_simulator)

# ----------------------------------------------------------------------------
# CPack Install
# ----------------------------------------------------------------------------

if (SEQAN_BUILD_SYSTEM STREQUAL "APP:indel_simulator")
  set (CPACK_PACKAGE_NAME "indel_simulator")
  set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "indel_simulator - Indel Simulation")
  set (CPACK_DEBIAN_PACKAGE_MAINTAINER "Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>")
  set (CPACK_PACKAGE_VENDOR "SeqAn Team, FU Berlin")

  seqan_configure_cpack_app (indel_simulator "indel_simulator")
endif (SEQAN_BUILD_SYSTEM STREQUAL "APP:indel_simulator")

