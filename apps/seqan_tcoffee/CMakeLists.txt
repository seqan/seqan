# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /core/apps/seqan_tcoffee/CMakeLists.txt
#
# CMakeLists.txt file for SeqAn::T-Coffee.
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (core_apps_seqan_tcoffee)
message (STATUS "Configuring core/apps/seqan_tcoffee")

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
add_executable (seqan_tcoffee seqan_tcoffee.cpp)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (seqan_tcoffee ${SEQAN_LIBRARIES})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------

# Set variables for installing, depending on the selected build type.
if (NOT SEQAN_PREFIX_SHARE_DOC)
  seqan_setup_install_vars (seqan_tcoffee)
endif (NOT SEQAN_PREFIX_SHARE_DOC)

# Install seqan_tcoffee in ${PREFIX}/bin directory
install (TARGETS seqan_tcoffee
         DESTINATION bin)

# Install non-binary files for the package to "." for app builds and
# ${PREFIX}/share/doc/seqan_tcoffee for SeqAn release builds.
install (FILES LICENSE
               README
         DESTINATION ${SEQAN_PREFIX_SHARE_DOC})
install (FILES example/seq.fa
         DESTINATION ${SEQAN_PREFIX_SHARE_DOC}/example)

# ----------------------------------------------------------------------------
# App Test
# ----------------------------------------------------------------------------

seqan_add_app_test (seqan_tcoffee)

# ----------------------------------------------------------------------------
# Setup Common Tool Description for Generic Workflow Nodes
# ----------------------------------------------------------------------------
	
# Include executable seqan_tcoffee in CTD structure.
set (SEQAN_CTD_EXECUTABLES ${SEQAN_CTD_EXECUTABLES} seqan_tcoffee CACHE INTERNAL "")

# ----------------------------------------------------------------------------
# CPack Install
# ----------------------------------------------------------------------------

if (SEQAN_BUILD_SYSTEM STREQUAL "APP:seqan_tcoffee")
  set (CPACK_PACKAGE_NAME "seqan_tcoffee")
  set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "SeqAn::T-Coffee - Multiple Sequence Alignment")
  set (CPACK_DEBIAN_PACKAGE_MAINTAINER "David Weese <david.weese@fu-berlin.de>")
  set (CPACK_PACKAGE_VENDOR "SeqAn Team, FU Berlin")

  seqan_configure_cpack_app (seqan_tcoffee "seqan_tcoffee")
endif (SEQAN_BUILD_SYSTEM STREQUAL "APP:seqan_tcoffee")

