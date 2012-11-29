# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /core/apps/mason/CMakeLists.txt
#
# CMakeLists.txt file for Mason.
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (core_apps_mason)
message (STATUS "Configuring core/apps/mason")

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
add_executable (mason
                mason.cpp
                mason.h
                simulate_454.h
                simulate_454_base_calling.h
                simulate_illumina.h
                simulate_sanger.h
                store_config.h
                util.h)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (mason ${SEQAN_LIBRARIES})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------

# Set variables for installing, depending on the selected build type.
if (NOT SEQAN_PREFIX_SHARE_DOC)
  seqan_setup_install_vars (mason)
endif (NOT SEQAN_PREFIX_SHARE_DOC)

# Install mason in ${PREFIX}/bin directory
install (TARGETS mason
         DESTINATION bin)

# Install non-binary files for the package to "." for app builds and
# ${PREFIX}/share/doc/mason for SeqAn release builds.
install (FILES LICENSE
               README
         DESTINATION ${SEQAN_PREFIX_SHARE_DOC})
install (FILES example/SRR018294.1.1M.error_dist
               example/SRR018311.1.1M.error_dist
               example/SRR026674.1.1M.error_dist
               example/SRR049254.1.1M.error_dist
         DESTINATION ${SEQAN_PREFIX_SHARE_DOC}/example)

# ----------------------------------------------------------------------------
# App Test
# ----------------------------------------------------------------------------

seqan_add_app_test (mason)

# ----------------------------------------------------------------------------
# CPack Install
# ----------------------------------------------------------------------------

if (SEQAN_BUILD_SYSTEM STREQUAL "APP:mason")
  set (CPACK_PACKAGE_NAME "mason")
  set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "Mason - A Read Simulator")
  set (CPACK_DEBIAN_PACKAGE_MAINTAINER "Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>")
  set (CPACK_PACKAGE_VENDOR "SeqAn Team, FU Berlin")

  seqan_configure_cpack_app (mason "Mason")
endif (SEQAN_BUILD_SYSTEM STREQUAL "APP:mason")

