# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /apps/micro_razers/CMakeLists.txt
#
# CMakeLists.txt file for MicroRazerS.
# ===========================================================================

cmake_minimum_required (VERSION 3.0.0)
project (seqan_apps_micro_razers CXX)
message (STATUS "Configuring apps/micro_razers")

set (SEQAN_APP_VERSION "1.0.11")

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
add_executable (micro_razers
               micro_razers.cpp
               outputFormat.h
               ../razers/paramChooser.h
               ../razers/param_tabs.h
               ../razers/param_tabs.cpp
               ../razers/param_tabs.inc)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (micro_razers ${SEQAN_LIBRARIES})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------

# Set variables for installing, depending on the selected build type.
if (NOT SEQAN_PREFIX_SHARE_DOC)
  seqan_setup_install_vars (micro_razers)
endif (NOT SEQAN_PREFIX_SHARE_DOC)

# Install micro_razers in ${PREFIX}/bin directory
install (TARGETS micro_razers
         DESTINATION ${CMAKE_INSTALL_BINDIR})

# Install non-binary files for the package to "." for app builds and
# ${PREFIX}/share/doc/micro_razers for SeqAn release builds.
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

seqan_add_app_test (micro_razers)

# ----------------------------------------------------------------------------
# CPack Install
# ----------------------------------------------------------------------------

if (SEQAN_BUILD_SYSTEM STREQUAL "APP:micro_razers")
  set (CPACK_PACKAGE_NAME "micro_razers")
  set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "MicroRazerS - Read Mapping With Controllable Sensitivity")
  set (CPACK_DEBIAN_PACKAGE_MAINTAINER "David Weese <david.weese@fu-berlin.de>")
  set (CPACK_PACKAGE_VENDOR "SeqAn Team, FU Berlin")

  seqan_configure_cpack_app (micro_razers "MicroRazerS")
endif (SEQAN_BUILD_SYSTEM STREQUAL "APP:micro_razers")

