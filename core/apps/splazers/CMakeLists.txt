# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /core/apps/splazers/CMakeLists.txt
#
# CMakeLists.txt file for SplazerS.
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (core_apps_splazers)
message (STATUS "Configuring core/apps/splazers")

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
add_executable (splazers splazers.cpp
                        outputFormat.h
                        paramChooser.h
                        razers.h
                        razers_matepairs.h
                        razers_parallel.h
                        razers_spliced.h
                        readSimulator.h)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (splazers ${SEQAN_LIBRARIES})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------

# Set variables for installing, depending on the selected build type.
if (NOT SEQAN_PREFIX_SHARE_DOC)
  seqan_setup_install_vars (splazers)
endif (NOT SEQAN_PREFIX_SHARE_DOC)

# Install splazers in ${PREFIX}/bin directory
install (TARGETS splazers
         DESTINATION bin)

# Install non-binary files for the package to "." for app builds and
# ${PREFIX}/share/doc/splazers for SeqAn release builds.
install (FILES LICENSE
               README
         DESTINATION ${SEQAN_PREFIX_SHARE_DOC})
install (FILES example/genome.fa
               example/reads.fa
               example/reads.fa.result
         DESTINATION ${SEQAN_PREFIX_SHARE_DOC}/example)

# ----------------------------------------------------------------------------
# App Test
# ----------------------------------------------------------------------------

seqan_add_app_test (splazers)

# ----------------------------------------------------------------------------
# CPack Install
# ----------------------------------------------------------------------------

if (SEQAN_BUILD_SYSTEM STREQUAL "APP:splazers")
  set (CPACK_PACKAGE_NAME "splazers")
  set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "SplazerS - Split Read Mapping")
  set (CPACK_DEBIAN_PACKAGE_MAINTAINER "David Weese <david.weese@fu-berlin.de>")
  set (CPACK_PACKAGE_VENDOR "SeqAn Team, FU Berlin")

  seqan_configure_cpack_app (splazers "SplazerS")
endif (SEQAN_BUILD_SYSTEM STREQUAL "APP:splazers")

