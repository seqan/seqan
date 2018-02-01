# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /apps/razers/CMakeLists.txt
#
# CMakeLists.txt file for RazerS.
# ===========================================================================

cmake_minimum_required (VERSION 3.0.0)
project (seqan_apps_razers CXX)
message (STATUS "Configuring apps/razers")

set (SEQAN_APP_VERSION "1.5.8")

if (SEQAN_TRAVIS_BUILD)
    message (STATUS "  Skipping build and test of razers on Travis CI.")
    return ()
endif (SEQAN_TRAVIS_BUILD)

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
add_executable (razers razers.cpp
                      razers.h
                      param_tabs.h
                      param_tabs.cpp
                      param_tabs.inc
                      outputFormat.h
                      paramChooser.h
                      razers_matepairs.h
                      razers_spliced.h
                      readSimulator.h
                      param_tabs.h)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (razers ${SEQAN_LIBRARIES})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------

# Set variables for installing, depending on the selected build type.
if (NOT SEQAN_PREFIX_SHARE_DOC)
  seqan_setup_install_vars (razers)
endif (NOT SEQAN_PREFIX_SHARE_DOC)

# Install razers in ${PREFIX}/bin directory
install (TARGETS razers
         DESTINATION ${CMAKE_INSTALL_BINDIR})

# Install non-binary files for the package to "." for app builds and
# ${PREFIX}/share/doc/razers for SeqAn release builds.
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

seqan_add_app_test (razers)

# ----------------------------------------------------------------------------
# Setup Common Tool Description for Generic Workflow Nodes
# ----------------------------------------------------------------------------

# Include executable razers in CTD structure.
set (SEQAN_CTD_EXECUTABLES ${SEQAN_CTD_EXECUTABLES} razers CACHE INTERNAL "")

# ----------------------------------------------------------------------------
# CPack Install
# ----------------------------------------------------------------------------

if (SEQAN_BUILD_SYSTEM STREQUAL "APP:razers")
  set (CPACK_PACKAGE_NAME "razers")
  set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "RazerS - Read Mapping With Controllable Sensitivity")
  set (CPACK_DEBIAN_PACKAGE_MAINTAINER "David Weese <david.weese@fu-berlin.de>")
  set (CPACK_PACKAGE_VENDOR "SeqAn Team, FU Berlin")

  seqan_configure_cpack_app (razers "RazerS")
endif (SEQAN_BUILD_SYSTEM STREQUAL "APP:razers")

